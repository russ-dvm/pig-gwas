#### FOR USE WITH GEMMA!
## Set variables
dir <- "~/snpChip/round_3/"
name <- "sal_shed_v2-6"
out_name <- paste(name, "_analysis", sep = "")

######## Don't touch
read_in <- paste(dir, name, "/", out_name, ".final.fam", sep = "")
out_dir <- paste(dir, name,"/", "covariates.txt", sep = "")

exp_data <- read.table(read_in, h = F)
exp_data$rank <- 1:nrow(exp_data)

master <- read.table("~/snpChip/testing/ref_files/salmonella_all_info.txt", h = T, sep = "\t", na.strings = c("", NA))

combo <- merge(exp_data, master, by.x = "V2", by.y = "newid")

## Reformat categorical variables into integers for subsequent analysis in GEMEMA.

combo$diet_int <- ifelse(combo$diet == "HC", combo$diet <- 1, combo$diet <- 2)
combo$trial_int <- ifelse(combo$trial == "summer", combo$trial_int <- 1, combo$trial_int <- 2)

combo$seas2_int <- ifelse(combo$seas2 == "spring", combo$seas2_int <- 1, ifelse(combo$seas2 == "summer", combo$seas2_int <- 2, ifelse(combo$seas2 == "fall", combo$seas2_int <- 3, combo$seas2_int <- 4)))
combo$seas3_int <- ifelse(combo$seas3 == "spring", combo$seas3_int <- 1, ifelse(combo$seas3 == "summer", combo$seas3_int <- 2, ifelse(combo$seas3 == "fall", combo$seas3_int <- 3, combo$seas3_int <- 4)))
combo$seas4_int <- ifelse(combo$seas4 == "spring", combo$seas4_int <- 1, ifelse(combo$seas4 == "summer", combo$seas4_int <- 2, ifelse(combo$seas4 == "fall", combo$seas4_int <- 3, combo$seas4_int <- 4)))
combo$seas5_int <- ifelse(combo$seas5 == "spring", combo$seas5_int <- 1, ifelse(combo$seas5 == "summer", combo$seas5_int <- 2, ifelse(combo$seas5 == "fall", combo$seas5_int <- 3, combo$seas5_int <- 4)))
combo$seas6_int <- ifelse(combo$seas6 == "spring", combo$seas6_int <- 1, ifelse(combo$seas6 == "summer", combo$seas6_int <- 2, ifelse(combo$seas6 == "fall", combo$seas6_int <- 3, combo$seas6_int <- 4)))

# According to C Schutt's analysis of the variables (using linear or logistic regression with V Farzan):
# "The significant variables (univar analysis) were:
### THESE WERE USED FOR ANALYSIS ROUND 1. THE SAME COVARIATES WERE USED FOR ALL TIME POINTS.
## Seropositivity (0/1): sal (categorical) and age (non-normal distribution) were significant.
## S/P (titres): sal, trial, season (all categorical), and age (non-normal distribution) were significant."
## manual covariate definition:
# covars <- c("farm")


## Perhaps more reproducible way is to use the covariate_selector function that has a predefined data frame of the covariates.
## Note that this script was written for GMMAT which requires the phenotype to be defined. GEMMA does not need this. Needs to be removed.
source("~/snpChip/scripts/covariate_selector.R")
covars <- assign_covars(name)
covars <- covars[-1]

## add in the initial ranking to ensure the output is in the same order as the original data
covars_rank <- c(covars, "rank")
gemma1 <- combo[,covars_rank]

## add in the required list of 1s (see gemma manual)
gemma1$int <- 1

## reorder according to the original input order
gemma2 <- gemma1[order(gemma1$rank),]

## select only the vector of 1s and the covariates, and omit the ranks
covars_int <- c("int", covars)
gemma3 <- gemma2[,covars_int]
covars
## write out for incorporation into the gemma model.
write.table(gemma3, file = out_dir, row.names = F, quote = F, col.names = F)

