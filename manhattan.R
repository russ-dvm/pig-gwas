library(tidyverse)


## Import data
dir <- "~/snpChip/testing_round_1/"
name <- "sal_sero_v5"
out_name <- "sal_sero_v5_analysis"
######## Don't touch

read_in <- paste(dir, name, "/output/", out_name, ".assoc.txt", sep = "")
write_out_man <- paste(dir, name, "/", out_name, "_manhattan.png", sep = "")
write_out_qq <- paste(dir, name, "/", out_name, "_qqPlot.png", sep = "")

salm <- read.table(read_in, h = T, sep = "\t")
g_title <- out_name
# salm <- read.table("~/snpChip/testing/russ/output/test.assoc.txt", h = T, sep = "\t")

## To use qqman the data needs to be in a DF with columsn "CHR", "BP", "P", and "SNP"
# salm_qq <- data.frame("CHR" = as.integer(salm$CHR), "BP" = as.integer(salm$BP), "P" = as.numeric(salm$P), "SNP" = as.factor(salm$SNP))
#For gemma results use this one:
salm_qq <- data.frame("CHR" = as.integer(salm$chr), "BP" = as.integer(salm$ps), "P" = as.numeric(salm$p_wald), "SNP" = as.factor(salm$rs))


###***TMP***###
#remove chr0, 21 and 22....
salm_qq_sub <- subset(salm_qq, salm_qq$CHR != 0 & salm_qq$CHR != 21 & salm_qq$CHR != 22)

## Perform the p adjustments
source("~/scripts/fdr.R")
qq_results <- fdr(salm_qq_sub, "P", 0.1)

#Check to make sure the two methods agreed with each other...
qq_results[qq_results$method.comp == F,]

## Add in start/stop coordinates - these all have NA P values so won't actually be plotted.
for (x in 1:20){
  qq_results <- rbind(qq_results, c(x, 1, NA, NA))
}
qq_results <- rbind(qq_results, c(1, 315321322, NA, NA))
qq_results <- rbind(qq_results, c(2, 162569375, NA, NA))
qq_results <- rbind(qq_results, c(3, 144787322, NA, NA))
qq_results <- rbind(qq_results, c(4, 143465943, NA, NA))
qq_results <- rbind(qq_results, c(5, 111506441, NA, NA))
qq_results <- rbind(qq_results, c(6, 157765593, NA, NA))
qq_results <- rbind(qq_results, c(7, 134764511, NA, NA))
qq_results <- rbind(qq_results, c(8, 148491826, NA, NA))
qq_results <- rbind(qq_results, c(9, 153670197, NA, NA))
qq_results <- rbind(qq_results, c(10, 79102373, NA, NA))
qq_results <- rbind(qq_results, c(11, 87690581, NA, NA))
qq_results <- rbind(qq_results, c(12, 63558571, NA, NA))
qq_results <- rbind(qq_results, c(13, 218635234, NA, NA))
qq_results <- rbind(qq_results, c(14, 153851969, NA, NA))
qq_results <- rbind(qq_results, c(15, 157681621, NA, NA))
qq_results <- rbind(qq_results, c(16, 86898991, NA, NA))
qq_results <- rbind(qq_results, c(17, 69701581, NA, NA))
qq_results <- rbind(qq_results, c(18, 61220071, NA, NA))
qq_results <- rbind(qq_results, c(19, 144288218, NA, NA)) #X chrom
qq_results <- rbind(qq_results, c(20, 1637716, NA, NA)) #Y chrom

## Adjust the data frame to generate points for manhattan plotting
source("~/scripts/manhattan.R")
salm_man <- manhattan(qq_results)


## Make chromosomes into factors (have to be numeric for the manhattan function, but better for coloring as a factor)
salm_man$CHR <- factor(salm_man$CHR, levels = unique(salm_man$CHR))

salm_test <- merge(salm_man, qq_results, by.x = "SNP", by.y = "SNP")

## Plot
man <- ggplot(salm_man, aes(x = pos, y = logp, color = CHR)) + 
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
  ggtitle(g_title)
# man

ggsave(man, file = write_out_man, units = "in", height = 3, width = 7)

 ## QQ plot
qqplot <- qqman::qq(salm_man$P)
lamdaEst <- GenABEL::estlambda(salm_man$P, method = "median")
main_title <- paste("Lambda estimate:", lamdaEst[1], sep = " ")
title(main = main_title)

dev.print(png, write_out_qq, height = 500, width = 500)

