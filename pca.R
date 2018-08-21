library(tidyverse)

#### Read in data
pca <- read.table("~/snpChip/testing/russ/plink-pca.eigenvec", h = F)
ei <- read.table("~/snpChip/testing/russ/plink-pca.eigenval", h = F)
deets <- read.table("~/snpChip/testing/ref_files/salmonella_all_info.txt", h = T, sep = "\t")

pca <- read.table("~/snpChip/testing/russ/gcta-pca.eigenvec", h = F)
ei <- read.table("~/snpChip/testing/russ/gcta-pca.eigenval", h = F)

#### Calculate the fraction of the total for each eigenvalue
ei$prop <- ei$V1/sum(ei$V1)*100

#### Figures
## Plot the PCA components
ggplot(ei, aes(x = rev(as.factor(V1)), y = prop)) + geom_bar(stat = "identity")  + theme_classic() +  scale_y_continuous(limits = c(0,100), expand = c(0,0)) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab("PCA component")

## Take a look at the first two eigenvectors. First create a dataframe that includes the eigenvectors along with phenotypic information (aka deets) about the samples.
a <- merge(pca, deets, by.x = "V2", by.y = "newid")
#Basic
ggplot(a, aes(x = V3, y = V4)) + geom_point() + xlab("PCA 1") + ylab("PCA 2") +theme_bw()
ggplot(a, aes(x = V3, y = V4, color = as.factor(farm), shape = as.factor(genetics))) + geom_point() + xlab("PCA 1") + ylab("PCA 2") +theme_bw()
ggplot(a, aes(x = V3, y = V4, color = sex)) + geom_point()+ xlab("PCA 1") + ylab("PCA 2") +theme_bw()
ggplot(a, aes(x = V3, y = V4, color = diet)) + geom_point()+ xlab("PCA 1") + ylab("PCA 2") +theme_bw()
ggplot(a, aes(x = V3, y = V4, color = trial)) + geom_point()+ xlab("PCA 1") + ylab("PCA 2") +theme_bw()
ggplot(a, aes(x = V3, y = V4, color = as.factor(pigpos))) + geom_point()+ xlab("PCA 1") + ylab("PCA 2") +theme_bw()
ggplot(a, aes(x = V3, y = V4, color = genetics)) + geom_point()+ xlab("PCA 1") + ylab("PCA 2") +theme_bw()
ggplot(a, aes(x = V3, y = V4, color = as.factor(farm))) + geom_point()+ xlab("PCA 1") + ylab("PCA 2") +theme_bw()

