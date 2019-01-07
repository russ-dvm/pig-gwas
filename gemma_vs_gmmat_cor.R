library(ggplot2)

dir <- "~/snpChip/round_3/"
sample <- "sal_sero_v5"

gmmat <- paste(dir, sample, "/", sample, "_gmmat_scores.txt", sep = "")
gemma <- paste(dir, sample, "/", "output", "/", sample, "_analysis.assoc.txt", sep = "")

gmmatResults <- read.table(gmmat, sep = "\t", h = T)
gemmaResults <- read.table(gemma, sep = "\t", h = T)

gmmatResultsOrdered <- gmmatResults[order(gmmatResults$PVAL, decreasing = F), ]
gmmatResultsOrdered$gmmat_rank <- c(1:nrow(gmmatResultsOrdered))

gemmaResultsOrdered <- gemmaResults[order(gemmaResults$p_wald, decreasing = F),]
gemmaResultsOrdered$gemma_rank <- c(1:nrow(gemmaResultsOrdered))

both <- merge(gemmaResultsOrdered, gmmatResultsOrdered, by.x = "rs", by.y = "SNP")
both1000 <- subset(both, gemma_rank <= 1000)
both100 <- subset(both, gemma_rank <= 100)
both10 <- subset(both, gemma_rank <= 10)


bothCor <- cor.test(both$gmmat_rank, both$gemma_rank)
cor.test.value <- signif(bothCor$estimate, 3)
cor.pvalue <- signif(bothCor$p.value, 3)

bothCor1000 <- cor.test(both1000$gmmat_rank, both1000$gemma_rank)
cor.test.value1000 <- signif(bothCor1000$estimate, 3)
cor.pvalue1000 <- signif(bothCor1000$p.value, 3)

bothCor100 <- cor.test(both100$gmmat_rank, both100$gemma_rank)
cor.test.value100 <- signif(bothCor100$estimate, 3)
cor.pvalue100 <- signif(bothCor100$p.value, 3)

bothCor10 <- cor.test(both10$gmmat_rank, both10$gemma_rank)
cor.test.value10 <- signif(bothCor10$estimate, 3)
cor.pvalue10 <- signif(bothCor10$p.value, 3)


comp <- ggplot(both, aes(x = gemma_rank, y = gmmat_rank)) + 
  geom_point() + 
  geom_smooth(method = lm) +
  ggtitle(paste(sample, ", ", "R = ", cor.test.value, ", ", "p = ", cor.pvalue)) +
  theme_bw() +
  theme(text = element_text(size = 5))

comp_1000 <- ggplot(both1000, aes(x = gemma_rank, y = gmmat_rank)) + 
  geom_point() + 
  geom_smooth(method = lm) +
  ggtitle(paste(sample, ", ", "R = ", cor.test.value1000, ", ", "p = ", cor.pvalue1000)) +
  theme_bw() +
  theme(text = element_text(size = 5))

comp_100 <- ggplot(both100, aes(x = gemma_rank, y = gmmat_rank)) + 
  geom_point() + 
  geom_smooth(method = lm) +
  ggtitle(paste(sample, ", ", "R = ", cor.test.value100, ", ", "p = ", cor.pvalue100)) +
  theme_bw() +
  theme(text = element_text(size = 5))

comp_10 <- ggplot(both10, aes(x = gemma_rank, y = gmmat_rank)) + 
  geom_point() + 
  geom_smooth(method = lm) +
  ggtitle(paste(sample, ", ", "R = ", cor.test.value10, ", ", "p = ", cor.pvalue10)) +
  theme_bw() +
  theme(text = element_text(size = 5))

# ggsave(paste(dir, sample, "/", sample, "_cor.png", sep = ""), comp, units = "in", width = 6, height = 6)
library(cowplot)
final <- plot_grid(comp, comp_1000, comp_100, comp_10)
save_plot(paste(dir, sample, "/", sample, "_cor.png", sep = ""), final, base_aspect_ratio = 1)  

