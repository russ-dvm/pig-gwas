library(ggplot2)
library(GMMAT)


#### SET VARIABLES ####

dir <- "~/snpChip/round_2/"
sample <- "sal_shed_v2-6"


## Automatically select the covariates based on the data frame in covariate_selector.R
## Check the source file to see what covariates were used for which timepoint, or to alter them.
source("~/snpChip/scripts/covariate_selector.R")
covars <- assign_covars(sample)

## THE MODEL
sal_or_sero <- ifelse(strsplit(sample, split = "_")[[1]][2] == "sero", "sero", "sal")
modelPhenoPos <- grep(sal_or_sero, covars)
modelPheno <- covars[modelPhenoPos]
modelCovarsList <- covars[-modelPhenoPos]
modelCovars <- paste(modelCovarsList, collapse = "+")
modelFinal <- paste(modelPheno, "~", modelCovars, sep = " ")
modelFinal

##DEPRECATED BUT KEPT FOR POSTERITY....the manual way of doing things...
# covars <- c("diet_int", "trial_int", "sal26") ## SHEDDING
# covars <- c("seas2_int", "sal2") ## SEROPOSITIVITY

#### IMPORT DATA ####
### Don't touch
out_name <- paste(sample, "_analysis", sep = "")
read_in_fam <- paste(dir, sample, "/", out_name, ".final.fam", sep = "")
read_in_grm <- paste(dir, sample, "/", "output", "/", out_name, "-grm.cXX.txt", sep = "")
missing_out <- paste(dir, sample,"/", "samples_missing_any_pheno.txt", sep = "")

exp_data <- read.table(read_in_fam, h = F)
exp_data$rank <- 1:nrow(exp_data)

master <- read.table("~/snpChip/testing/ref_files/salmonella_all_info.txt", h = T, sep = "\t", na.strings = c("", NA))

combo <- merge(exp_data, master, by.x = "V2", by.y = "newid")

#### GENERATE PHENO FILE ####
## Reformat categorical variables into integers for subsequent analysis in GEMMA - needed for GMMAT? - yes.

combo$diet_int <- ifelse(combo$diet == "HC", combo$diet <- 1, combo$diet <- 2)
combo$trial_int <- ifelse(combo$trial == "summer", combo$trial_int <- 1, combo$trial_int <- 2)

combo$seas2_int <- ifelse(combo$seas2 == "spring", combo$seas2_int <- 1, ifelse(combo$seas2 == "summer", combo$seas2_int <- 2, ifelse(combo$seas2 == "fall", combo$seas2_int <- 3, combo$seas2_int <- 4)))
combo$seas3_int <- ifelse(combo$seas3 == "spring", combo$seas3_int <- 1, ifelse(combo$seas3 == "summer", combo$seas3_int <- 2, ifelse(combo$seas3 == "fall", combo$seas3_int <- 3, combo$seas3_int <- 4)))
combo$seas4_int <- ifelse(combo$seas4 == "spring", combo$seas4_int <- 1, ifelse(combo$seas4 == "summer", combo$seas4_int <- 2, ifelse(combo$seas4 == "fall", combo$seas4_int <- 3, combo$seas4_int <- 4)))
combo$seas5_int <- ifelse(combo$seas5 == "spring", combo$seas5_int <- 1, ifelse(combo$seas5 == "summer", combo$seas5_int <- 2, ifelse(combo$seas5 == "fall", combo$seas5_int <- 3, combo$seas5_int <- 4)))
combo$seas6_int <- ifelse(combo$seas6 == "spring", combo$seas6_int <- 1, ifelse(combo$seas6 == "summer", combo$seas6_int <- 2, ifelse(combo$seas6 == "fall", combo$seas6_int <- 3, combo$seas6_int <- 4)))

## Pheno for v2-5
combo$sal25 <- ifelse((combo$sal2 == 0 | is.na(combo$sal2)) & (combo$sal3 == 0 | is.na(combo$sal3)) & (combo$sal4 == 0 | is.na(combo$sal4)) & (combo$sal5 == 0 | is.na(combo$sal5)), 0, 1)

# Pheno for v2-6
combo$sal26 <- ifelse((combo$sal2 == 0 | is.na(combo$sal2)) & (combo$sal3 == 0 | is.na(combo$sal3)) & (combo$sal4 == 0 | is.na(combo$sal4)) & (combo$sal5 == 0 | is.na(combo$sal5)) & (combo$sal6 == 0 | is.na(combo$sal6)), 0, 1)

## Remove any samples that are missing any of the covariate information.
## AUTOMATED sample removal
## Determine column values based on the covar names
covar_columns <- lapply(covars, function(x) {which(colnames(combo) == x)})
## Identify missing information in all of the columns
missing_by_covar <- lapply(covar_columns, function(x) {is.na(combo[,x])})
## Collapse the list to one final list in which any missing sample is keps
missingAny <- Reduce("|", missing_by_covar)
## Identify samples with any of the missing info
missing <- combo[missingAny, c(2,1)]
## Old, manual method (depracated): 
# missing <- combo[is.na(combo$sal2) | is.na(combo$seas2_int), c(2,1)]


write.table(missing, file = missing_out, row.names = F, col.names = F, quote = F, sep = "\t")
setwd(paste(dir, sample, sep = ""))
old_plink <- paste(out_name, ".final", sep = "")
new_plink <- paste(out_name, ".final.gmmat", sep = "")
system(paste("plink2 --bfile", old_plink, "--remove samples_missing_any_pheno.txt --make-bed --out", new_plink))

covars_rank <- c(covars, "rank")
pheno1 <- combo[,covars_rank]

pheno2 <- pheno1[order(pheno1$rank),]

pheno3 <- pheno2[,covars]

#### Import the GRM from GEMMA (does not require a covariates file)
grm <- as.matrix(read.table(read_in_grm))

#### SPECIFY THE GENOTYPES ####
geno_wald <- paste(sample, "_analysis.final", sep = "")
gmmat_score <- paste(sample, "_gmmat_scores.txt", sep = "")

#### FIT THE NULL MODEL ####

model <- glmmkin(fixed = modelFinal, data = pheno3, kins = grm, family = binomial(link = "logit"))


glmm.score(model, infile = new_plink, outfile = gmmat_score)


gmmatResults <- read.table(gmmat_score, h = T, sep = "\t")

# glmm.wald(fixed = sero2 ~ age2 + sal2, data = pheno3, kins = grm, infile = geno_wald, family = binomial(link = "logit"), snps = "AX-116662071")

#### PLOTS ####
## To use qqman the data needs to be in a DF with columsn "CHR", "BP", "P", and "SNP"
qqGmmat <- data.frame("CHR" = as.integer(gmmatResults$CHR), "BP" = as.integer(gmmatResults$POS), "P" = as.numeric(gmmatResults$PVAL), "SNP" = as.factor(gmmatResults$SNP))

source("~/scripts/fdr.R")
qqGmmatFdr <- fdr(qqGmmat, "P", 0.1)

qqGmmatFdr[qqGmmatFdr$method.comp == F,]

## Add in start/stop coordinates - these all have NA P values so won't actually be plotted.
for (x in 1:20){
  qqGmmatFdr <- rbind(qqGmmatFdr, c(x, 1, NA, NA))
}
qqGmmatFdr <- rbind(qqGmmatFdr, c(1, 315321322, NA, NA))
qqGmmatFdr <- rbind(qqGmmatFdr, c(2, 162569375, NA, NA))
qqGmmatFdr <- rbind(qqGmmatFdr, c(3, 144787322, NA, NA))
qqGmmatFdr <- rbind(qqGmmatFdr, c(4, 143465943, NA, NA))
qqGmmatFdr <- rbind(qqGmmatFdr, c(5, 111506441, NA, NA))
qqGmmatFdr <- rbind(qqGmmatFdr, c(6, 157765593, NA, NA))
qqGmmatFdr <- rbind(qqGmmatFdr, c(7, 134764511, NA, NA))
qqGmmatFdr <- rbind(qqGmmatFdr, c(8, 148491826, NA, NA))
qqGmmatFdr <- rbind(qqGmmatFdr, c(9, 153670197, NA, NA))
qqGmmatFdr <- rbind(qqGmmatFdr, c(10, 79102373, NA, NA))
qqGmmatFdr <- rbind(qqGmmatFdr, c(11, 87690581, NA, NA))
qqGmmatFdr <- rbind(qqGmmatFdr, c(12, 63558571, NA, NA))
qqGmmatFdr <- rbind(qqGmmatFdr, c(13, 218635234, NA, NA))
qqGmmatFdr <- rbind(qqGmmatFdr, c(14, 153851969, NA, NA))
qqGmmatFdr <- rbind(qqGmmatFdr, c(15, 157681621, NA, NA))
qqGmmatFdr <- rbind(qqGmmatFdr, c(16, 86898991, NA, NA))
qqGmmatFdr <- rbind(qqGmmatFdr, c(17, 69701581, NA, NA))
qqGmmatFdr <- rbind(qqGmmatFdr, c(18, 61220071, NA, NA))
qqGmmatFdr <- rbind(qqGmmatFdr, c(19, 144288218, NA, NA)) #X chrom
qqGmmatFdr <- rbind(qqGmmatFdr, c(20, 1637716, NA, NA)) #Y chrom

## Adjust the data frame to generate points for manhattan plotting
source("~/scripts/manhattan.R")
qqGmmatMan <- manhattan(qqGmmatFdr)


## Make chromosomes into factors (have to be numeric for the manhattan function, but better for coloring as a factor)
qqGmmatMan$CHR <- factor(qqGmmatMan$CHR, levels = unique(qqGmmatMan$CHR))

g_title <- paste(out_name, "gmmat/logistic")

## Manhattan Plot
man <- ggplot(qqGmmatMan, aes(x = pos, y = logp, color = CHR)) + 
  geom_point() +
  theme_classic() +
  scale_x_continuous(minor_breaks = minor_ticks, breaks = ticks, labels = labs) +
  scale_color_identity() +
  scale_y_continuous(expand = c(0,0), limits = c(0,8)) +
  xlab("Chromosome") +
  ylab("-log(p)") +
  geom_hline(yintercept = 6.30103, colour = "red") +
  annotate("text", label = "p < 5x10-7", y = 6.30103, x = 2596603306, vjust = -0.5) +
  geom_hline(yintercept = 4.30103, colour = "red", linetype = "dashed") +
  annotate("text", label = "p < 5x10-5", y = 4.3013, x = 2596603306, vjust = -0.5) + 
  ggtitle(paste(g_title, modelFinal, sep = "  |  "))
# man


write_out_man_gmmat <- paste(dir, sample, "/", out_name, "_gmmat_manhattan1.png", sep = "")
ggsave(man, file = write_out_man_gmmat, units = "in", height = 3, width = 7)

## QQ plot
qqplot <- qqman::qq(qqGmmatMan$P)
lamdaEst <- GenABEL::estlambda(qqGmmatMan$P, method = "median")
main_title <- paste("Lambda estimate:", lamdaEst[1], sep = " ")
title(main = main_title)

write_out_qq_gmmat <- paste(dir, sample, "/", out_name, "_gmmat_qqPlot.png", sep = "")

dev.print(png, write_out_qq_gmmat, height = 500, width = 500)

