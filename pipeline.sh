#### Set of commands to process Salmonella GWAS data ####

## 1. Convert the output of Axiom Analysis Suite to BED file.
# 1a. many samples are missing sex info - thus use the --allow-no-sex flag.
# 1b. plink defaults to 1 for controls, 2 for cases. --1 tells plink that controls are 0 and cases are 1.
# 1c. --chr-set sets the number of autosomes to 18, and the chrX to 18+1, chrY = 18+2, pseudoautosomal region of X to 18+3, and chrM to 18+4
plink2 --file complete --make-bed --allow-no-sex --1 --chr-set 18 --out complete


## 2. Subset the data to include only the specific desired samples
# 2a. Include MAF filtering here (could be included above)
# 2b. Keep-file should be two column tab-separated, Family_ID and Individual_ID (include a header)
# 2c. Can also exclude some SNPs at this step: e.g. in initial testing, there were 105 SNPs mapped to chr0,21,and 22, but the porcine karyotype is 18,XY... can filter those out if desired, but **note** 21 may be pseudoautosomal X and 22 may be chrM. Not sure about 0.
# awk '$1==0 || $1==21 || $1 ==22' complete.map >> exclude.snps
# It might be worth excluding these along with the X chromosome. 
# 2d. Check-sex. Plink2 has a build in tool for this, but I haven't quite figured out i) how to properly code the pseudoautosomal region and ii) what the coordinates of the pseudoautosomal region actually are. 
# https://doi.org/10.1159/000351310 has identified one border at  X:6743567, but the second (tail) border is still unknown... Could just use the end of X? X:144288218 - tried this, didn't help at all.
# Alternatively, AAS outputs computed sex, which can be compared to phenotype if desired using the script "check_gender.py" - need to do a little file finaggling (eg cut -f 2,5 complete.ped > known.sex, cut -f 1,8 sample_info.txt | sed 's/_.*CEL//' | sed 's/female/2/' sed 's/male/1/' > computed.sex)
plink2 --bfile complete --chr-set 18 --make-bed --allow-no-sex --keep salmonella_list.txt --exclude exclude.snps --out salmonella --maf 0.05


## 3. Calculate missingness. This command will output two files, we are interested mostly in file.imiss. This file will tell us how many SNPs failed genotyping for each sample (N_MISS) as well as the fraction (F_MISS). The .lmiss file contains information sorted by SNP.
# 3a. The data can be visualized in R with the script "qc_figs.R" 
plink2 --bfile salmonella --chr-set 18  --missing --out salmonella.missing


## 4. In preliminary testing, all samples were genotyped at > 96% percent. If this changes, then use the --remove flag to omit samples at this step, and then re-run the missing command to check out SNP level QC. Same concept for SNP level missingness - if there are SNPs that were poorly genotyped across a certain percentage of samples, these can be removed with --exclude (but run the next step first).


## 5. Check genotype call rates between cases/controls.
# 5a. Output is a file with fraction missing in cases and controls and Fisher's exact test used to determine if the difference is significant. Anderson et al 2010 use a cutoff of p < 0.00001, not clear how that was chosen. Can also do multiple testing adjustment and see if any differences are significant (coded in "qc_figs.R")
plink2 --bfile salmonella --chr-set 18  --test-missing --out salmonella.missing.case-control


## 6. Perform LD-pruning - 50 kb window, 5kb step, 0.2 r2
# 6a. The LD pruned dataset will be used to a) generate the GRM and b) perform a PCA
# 6b. First step is performing the LD, and the second step excludes any of the SNPs that were found to be in LD according to the above parameters.
plink2 --bfile salmonella --chr-set 18  --indep-pairwise 50 5 0.2 --make-founders --out salmonella.ld
plink2 --bfile salmonella --exclude salmonella.ld.prune.out --make-bed --out salmonella.pruned 


## 7. PCA. Can specify the number of eigenvalues after the --pca flag if desired (default is 20). The PCA can be visualized with the R script "pca.R"
plink2 --bfile salmonella.pruned --chr-set 18  --pca --out salmonella.pca


## 8. Construct a GRM 
# 8a. using GCTA
gcta64 --bfile <bed-file> --make-grm --out <grm> --autosome-num 18 --maf 0.05
# 8b. or use plink2
plink2 --bfile salm --chr-set 18  --make-rel square --out test-grm


# 8c. Construct an additional GRM using the above GCTA-generated GRM to account for some family-ness.
gcta64 --grm <grm> --make-bK 0.05 --out <out> --autosome-num 18 --maf 0.05


## 9. Perform GREML analysis *note no --bfile option*
# 9a. GREML estimates the variance explained by the SNPs used to construct/estimate the GRM
gcta64 --grm <grm> --reml --grm-adj 0 --pheno <pheno-file> --autosome-num 18 --out <out-file> --maf 0.05


# Perform GREML on family data
# Include the grm-adj??
gcta64 --reml --pheno <pheno-file> --autosome-num 18 --maf 0.05 --mgrm <txt-file with list of grms> --grm-adj 0 


## 10. MLMA analysis - this doesn't seem to accept any of our external variables (fixed/random)
gcta64 --bfile <bed-file> --mlma --grm <grm> --pheno <pheno-file> --autosome-num 18 --out <out-file> --maf 0.05
