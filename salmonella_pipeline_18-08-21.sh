#### Set of commands to process Salmonella GWAS data ####

## Preliminary processes

## set variables:
base=$1
dir=/home/russ/snpChip
out=$2


# cp ~/Dropbox/temp/$1* .
cp ~/snpChip/ref_files/sscrofa11.rsids .

## Create some tracking files
touch excluded_info.txt
echo -e "reason\tnumber" > excluded_info.txt

## Remove SNPs present on unknown chromosomes. In initial testing, there were 105 SNPs mapped to chr0,21,and 22, but the porcine karyotype is 18,XY... can filter those out if desired, but **note** 21 may be pseudoautosomal X and 22 may be chrM. Not sure about 0. Either way these are funky and difficult to reconcile. Probably best to exclude.
echo -e "\nRemoving SNPs that fall on unknown chromosomes (0, 21, and 22)"
awk '$1==0 || $1==21 || $1 ==22' $base.map > exclude.snps
awk '$1==0 || $1==21 || $1 ==22' $base.map > unknown.chr.snps
exclude_no=$(wc -l exclude.snps | cut -d ' ' -f 1)
echo -e "Unknown_chr"'\t'"$exclude_no" >> excluded_info.txt

# Also exclude all SNPs that don't have an rsID (these could be user specific SNPs. BL probably won't like this. Can always remove this step and re-run.)
# First, get AXids from AAS output, found in the .map file:
cut -f 2 $base.map | sed '/#/d' > expt.axid
# List of AXIds and corresponding rsids are found in a file output from AAS that I called "all_snp_info" and found in the ref_files.
cut -f 1,812 "$dir"/testing/ref_files/all_snp_info.txt > ref.ax-rs
# Remove headers from ref.ax-rsid
sed -i '/#/d' ref.ax-rs
sed -i '1d' ref.ax-rs

# Run ax-to-rsid.py
echo -e "\nConverting AXIOM ids to rsIDs, and excluding any SNP that lacks an rsID"
python "$dir"/scripts/ax-to-rsid.py expt.axid ref.ax-rs > "$base".ax-rsids
# Determine which axids do not have an rsID AND are presnt in our data:
no_rsid_no=$(grep -c no_match "$base".ax-rsids)
echo -e "no_rsid"'\t'"$no_rsid_no" >> excluded_info.txt

# Add those to the exclusion list (comment out this line if we want to include SNPs with no rsID):
grep no_match "$base".ax-rsids | cut -f 1 >> exclude.snps
grep no_match "$base".ax-rsids > missing-rsids.txt
# Then, check that all rsIDs map to Sscrofa 11:
# First, create a file with only rsIDs:
grep rs "$base".ax-rsids | cut -f 2 > sscrofa10.rsids
# Run remap_rsIDs.py (since this takes forever, made this optional):
echo -e '\n'
read -p "Re-map rsIDs to Sscrofa11? Warning: this will take awhile. Select no if this has already been done (y or n)" remap
if [[ $remap = "y" ]]
then
	echo -e "\nRe-mapping rsIDs from Sscrofa10 to Sscrofa11. This will take awhile (hours)..."
	python "$dir"/scripts/remap_rsIDs.py sscrofa10.rsids >> sscrofa11.rsids
fi
# How many, and which, SNPs don't match to Sscrofa11?
echo -e "\nExcluding SNPs that do not map to Sscrofa11 (ie do not have an rsID in Sscrofa11)"
not_in_ssc11=$(grep -c FAIL sscrofa11.rsids)
echo -e "not_in_ssc11"'\t'"$not_in_ssc11" >> excluded_info.txt
grep FAIL sscrofa11.rsids > missing-in-sscrofa11.rsids
cat missing-in-sscrofa11.rsids >> exclude.snps

## Now, create the binary files used in PLINK.
## 1. Convert the output of Axiom Analysis Suite to BED file.
# 1a. some samples are missing sex info - thus use the --allow-no-sex flag.
# 1b. plink defaults to 1 for controls, 2 for cases. --1 tells plink that controls are 0 and cases are 1.
# 1c. --chr-set sets the number of autosomes to 18, and the chrX to 18+1, chrY = 18+2, pseudoautosomal region of X to 18+3, and chrM to 18+4. This is producing some odd results - omit for now.
# 1d. Exclude the SNPs from above
# 1e. set the minor allele frequency
plink2 --file "$base" --make-bed --allow-no-sex --1 --out $2 --maf 0.05 --exclude exclude.snps

## 2. If all the data has been exported, you can subset only the specific samples that you want. Subset the data to include only the specific desired samples
# 2b. Keep-file should be two column tab-separated, Family_ID and Individual_ID (include a header)
# 2d. Check-sex (optional, and currently not working). Plink2 has a build in tool for this, but I haven't quite figured out i) how to properly code the pseudoautosomal region and ii) what the coordinates of the pseudoautosomal region actually are. 
# https://doi.org/10.1159/000351310 has identified one border at  X:6743567, but the second (tail) border is still unknown... Could just use the end of X? X:144288218 - tried this, didn't help at all.
# Alternatively, AAS outputs computed sex, which can be compared to phenotype if desired using the script "check_gender.py" - need to do a little file finaggling (eg cut -f 2,5 complete.ped > known.sex, cut -f 1,8 sample_info.txt | sed 's/_.*CEL//' | sed 's/female/2/' sed 's/male/1/' > computed.sex)
# plink2 --bfile complete --make-bed --allow-no-sex --keep salmonella_list.txt --out salmonella

##### INDIVIDUAL LEVEL QC ######
## 3. Calculate missingness. This command will output two files, we are interested mostly in file.imiss. This file will tell us how many SNPs failed genotyping for each sample (N_MISS) as well as the fraction (F_MISS). The .lmiss file contains information sorted by SNP.
# 3a. The data can be visualized in R with the script "qc_figs.R" 
# Note: it's a good idea to re-run these sample and snp level metrics, as removing the SNPs at the earlier stage changes the overall results. For example, a sample that was 90% genotyped before may no longer be after the removal of the weird SNPs above.
plink2 --bfile $2 --missing --out $2.missing

## Assuming a 90% genotyping rate:
awk '$6 >= 0.10' $2.missing.imiss | sed '1d' |  tr -s ' ' '\t' | sed 's/^\t//' | cut -f 1,2 > failed.samples.txt


## 4. Use the --remove flag to omit samples at this step.

plink2 --bfile $2 --remove failed.samples.txt --make-bed --out $2.samples.pruned

##### SNP LEVEL QC ######
## NOTE: All individuals requiring removal should be removed before proceeding to SNP filtering. ##

plink2 --bfile $2.samples.pruned --missing --out $2.samples.pruned.missing

awk '$5 >= 0.05' $2.samples.pruned.missing.lmiss | sed '1d' | tr -s ' ' '\t' | sed 's/^\t//' | cut -f 1,2 > failed.snps.txt

## 5. Check genotype call rates between cases/controls.
# 5a. Output is a file with fraction missing in cases and controls and Fisher's exact test used to determine if the difference is significant. Anderson et al 2010 use a cutoff of p < 0.00001, not clear how that was chosen. Can also do multiple testing adjustment and see if any differences are significant (coded in "qc_figs.R")
# 5b. Not relevant if that data is quantitative. Thus changed into a switch:
read -p "Is the data case-control or quantitative? (cc or q)" type
if [[ $type == "cc" ]]
then
	plink2 --bfile $2.samples.pruned --test-missing --out $2.missing.case-control
	awk '$5 < 0.00001' $2.missing.case-control.missing | tr -s ' ' '\t' | sed 's/^\t//' | cut -f 1,2 >> failed.snps.txt
fi


## Remove failed snps:
plink2 --bfile $2.samples.pruned --exclude failed.snps.txt --maf 0.05 --make-bed --out $2.final

## 6. Perform LD-pruning - 50 kb window, 5kb step, 0.2 r2
# 6a. The LD pruned dataset will be used to perform a PCA
# 6b. First step is performing the LD, and the second step excludes any of the SNPs that were found to be in LD according to the above parameters.
# plink2 --bfile salmonella  --indep-pairwise 50 5 0.2 --make-founders --out salmonella.ld

## 7. PCA. Can specify the number of eigenvalues after the --pca flag if desired (default is 20). The PCA can be visualized with the R script "pca.R"
# plink2 --bfile salmonella.pruned --pca --out salmonella.pca

### COVARIATES
read -p "Pausing to make sure you have a covariates file prepared. Hit enter when you do"


### GEMMA
# make GRM
gemma -bfile $2.final -gk 1 -o $2-grm
#Perform single snp analysis
gemma -bfile $2.final -k output/$2-grm.cXX.txt -lmm 4 -c covariates.txt -o $2
# Covariates can be generated using the covariate_generator.R script. Will have to re-construct a PED file with updated sample numbers.

## Generate a list of sig SNPs only:
awk '$12 <= 0.00005 || $13 <= 0.00005 || $14 <= 0.00005' output/"$2".assoc.txt > output/"$2".sig.snps.txt