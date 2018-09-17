library(tidyverse)

#PED
# exp_data <- read.table("~/snpChip/testing/corinne/corinne.ped", sep = "\t")
# ped <- exp_data[,c(1:6)]
#FAM
dir <- "~/snpChip/testing_round_1/"
name <- "sal_titre_v5"
out_name <- "sal_titre_v5_analysis"
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
# # Seropositivity (0/1): sal (categorical) and age (non-normal distribution) were significant.
# # S/P (titres): sal, trial, season (all categorical), and age (non-normal distribution) were significant."


## Enter the name(s) of the covariates
covars <- c("age5", "sal5", "trial_int", "seas5_int")


covars_rank <- c(covars, "rank")
gemma1 <- combo[,covars_rank]
gemma1$int <- 1

gemma2 <- gemma1[order(gemma1$rank),]

covars_int <- c("int", covars)
gemma3 <- gemma2[,covars_int]

write.table(gemma3, file = out_dir, row.names = F, quote = F, col.names = F)

