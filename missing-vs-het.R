library(tidyverse)

data <- read.table("~/snpChip/testing/plink_test.het", h = T)
missing <- read.table("~/snpChip/testing/plink_test.imiss", h = T)

data$obsHetRate <- (data$N.NM.-data$O.HOM.)/data$N.NM. 
standdev <- sd(data$obsHetRate)

final <- merge(data, missing)

a <- ggplot(final, aes(x = F_MISS, y = obsHetRate)) + 
  geom_point() + 
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  scale_x_log10(limits = c(0.0001, 1), breaks = c(0.0001,0.001,0.01,0.1,1)) +
  geom_vline(xintercept = 0.03, linetype = "dashed", colour = "blue") + 
  geom_hline(yintercept = (mean(data$obsHetRate) + 3*standdev), linetype = "dashed", colour = "blue") + 
  geom_hline(yintercept = (mean(data$obsHetRate) - 3*standdev), linetype = "dashed", colour = "blue") + 
  scale_y_continuous(limits = c(0, 0.5)) +
  ylab("Observed heterozygosity rate") + 
  xlab("Proportion of missing genotypes") +
  annotation_logticks(sides = "b")
a
ggsave(a, file = "~/snpChip/testing/missing-vs-heterozygosity.pdf")

failed <- subset(final, obsHetRate < mean(data$obsHetRate) - 3*standdev | obsHetRate > mean(data$obsHetRate) + 3*standdev | F_MISS > 0.03)
failed <- failed[c(1,2)]
write.table(failed, file = "~/snpChip/testing/fail-imiss-het-qc.txt", row.names = F, sep = "\t", quote = F)

t <- function(x){
  print(10^x)
}
