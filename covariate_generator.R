library(tidyverse)

#PED
# exp_data <- read.table("~/snpChip/testing/corinne/corinne.ped", sep = "\t")
# ped <- exp_data[,c(1:6)]
#FAM
exp_data <- read.table("~/snpChip/testing/russ/salmonella.final.fam", sep = "\t", h = F)

exp_data$rank <- 1:nrow(exp_data)

master <- read.table("~/snpChip/testing/ref_files/salmonella_all_info.txt", h = T, sep = "\t")

combo <- merge(exp_data, master, by.x = "V2", by.y = "newid")

gemma1 <- combo[,c(7,15,79,46)]
gemma1$int <- 1

gemma2 <- gemma1[order(gemma1$rank),]

gemma3 <- gemma2[,c(5,2,3,4)]

write.table(gemma3, file = "~/snpChip/testing/corinne/covariates_age.txt", row.names = F, quote = F, col.names = F)

