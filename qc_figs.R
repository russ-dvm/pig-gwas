library(tidyverse)

## Read in data
imiss <- "~/snpChip/testing/salmonella.missing.imiss"
lmiss <- "~/snpChip/testing/salmonella.missing.lmiss"
missing <- read.table(imiss, h = T)
snp.missing <- read.table(lmiss, h = T)
sex <- read.table("~/snpChip/testing/salmonella.sex.sexcheck", h = T) # plink sex
sex.inhouse <- read.table("~/snpChip/testing/sex.check", h = F)

cc.missing <- read.table("~/snpChip/testing/salmonella.missing.case-control.missing", h = T)

## Sample level figures
ggplot(missing, aes(x = N_MISS)) + geom_histogram(stat = "bin", binwidth = 10) +
  xlab("Number of SNPs not genotyped") + ylab("Number of pigs")

ggplot(missing, aes(x = F_MISS*100)) + geom_histogram(stat = "bin") +
  xlab("Percent SNPs not genotyped") + ylab("Number of pigs")

subset(missing, F_MISS*100 >=3)


## SNP level figures
ggplot(snp.missing, aes(x = N_MISS)) + geom_histogram(stat = "bin")
ggplot(snp.missing, aes(x = F_MISS*100)) + geom_histogram(stat = "bin")

subset(missing, F_MISS*100 >= 2)


## Case-control figures
source("~/scripts/fdr.R")
head(cc.missing)
a <- fdr(cc.missing, "P", 0.05)

## Sex
ggplot(sex, aes(x = F)) + geom_histogram(bins = 100) + scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.05))
ggplot(sex.inhouse[!is.na(sex.inhouse$V2),], aes(x = V2)) + geom_histogram(stat = "count", fill = "light grey") + stat_count(geom="text", colour=" black", size=3.5, aes(label=..count..), position=position_stack(vjust=0.5)) + theme_classic() + scale_y_continuous(expand = c(0,0), limits = c(0, 800))


