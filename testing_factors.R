library(tidyverse)
library(lme4)

## Test the significance/associations of the various factors of the expt against salmonella status to decide which to include.

## Read in data

info <- read.table("~/snpChip/testing/salmonella_all_info.txt", h = T, sep = "\t")

## Interpret salmonella status as a factor
info$pigpos_fact <- as.factor(info$pigpos)

## Remove pigs with no salmonella info
paste("Num of pigs without Salmonella info:", sum(is.na(info$pigpos)))
info <- info[!is.na(info$pigpos),]

## Let's have a look at the percentage of each variable that is positive for disease.
sex1 <- ggplot(info, aes(x = pigpos_fact, fill = sex)) + geom_histogram(stat = "count", position = "dodge")
sex1
sex_percentage <- info %>% group_by(sex) %>% summarise(dis_perc = mean(pigpos) * 100) %>%
  ungroup()
sex2 <- ggplot(sex_percentage, aes(x = sex, y = dis_perc)) + geom_bar(stat = "identity") + ylim(0,100) + ylab("Percent positive for Salmonella")
sex2
## These figures show us that a slightly higher percentage of female animals were positive for salmonella, but looks like there were a higher number of males in total in the study.

## Perform the regression on sex:
sex.glm <- glm(info$pigpos_fact ~ info$sex, family = "binomial")

## Info on interpretation from pagepiccini.com

## Remember this is in logit space, so we are looking at the probability that y will happen given x. Value of 0 = 50%, negative = less than 50%, positive value = greater than 50% chance.

## Looking at the estimate of the intercept will tell you the probability of the first factor (in the case of sex, female) being Salmonella positive.

## So in this case, (intercept) = -0.4031 means that females have a slightly less than 50% chance of being Salmonella positive.
## Looking at our variable, sex, the value of -0.325 suggests that males have a lower probability of being positive than females.

summary(sex.glm)


## Now let's check some other variables (farm, parity, trial, diet, etc)
diet1 <- ggplot(info, aes(x=pigpos_fact, fill = diet)) + geom_histogram(stat = "count", position = "dodge")
diet_percentage <- info %>% group_by(diet) %>% summarise(dis_perc = mean(pigpos) * 100) %>%
  ungroup()
diet2 <- ggplot(diet_percentage, aes(x = diet, y = dis_perc)) + geom_bar(stat = "identity") + ylim(0,100) + ylab("Percent positive for Salmonella")

## Make the model
diet.glm <- glm(info$pigpos_fact ~ info$diet, family = "binomial")
summary(diet.glm)

## Interpretation: our variable, diet, does not have a signiifcant effect on Salmonella status. I'm not sure what to do with the intercept in this case.

## Variable: trial
trial1 <- ggplot(info, aes(x=pigpos_fact, fill = trial)) + geom_histogram(stat = "count", position = "dodge")
trial_percentage <- info %>% group_by(trial) %>% summarise(dis_perc = mean(pigpos) * 100) %>%
  ungroup()
trial2 <- ggplot(trial_percentage, aes(x = trial, y = dis_perc)) + geom_bar(stat = "identity") + ylim(0,100) + ylab("Percent positive for Salmonella")

## Make the model
trial.glm <- glm(info$pigpos_fact ~ info$trial, family = "binomial")
summary(trial.glm)


## Variable: farm
farm1 <- ggplot(info, aes(x=pigpos_fact, fill = as.factor(farm))) + geom_histogram(stat = "count", position = "dodge")
farm_percentage <- info %>% group_by(farm) %>% summarise(dis_perc = mean(pigpos) * 100) %>%
  ungroup()
farm2 <- ggplot(farm_percentage, aes(x = farm, y = dis_perc)) + geom_bar(stat = "identity") + ylim(0,100) + ylab("Percent positive for Salmonella")

## Make the model
farm.glm <- glm(info$pigpos_fact ~ as.factor(info$farm), family = "binomial")
summary(farm.glm)

## Interpretation: major effects of farm.

## Variable: genetics company 
## This is likely related to the farm data - ie I think each farm has one genetics company
gen1 <- ggplot(info, aes(x=pigpos_fact, fill = genetics)) + geom_histogram(stat = "count", position = "dodge")
genetics_percentage <- info %>% group_by(genetics) %>% summarise(dis_perc = mean(pigpos) * 100) %>%
  ungroup()
gen2 <- ggplot(genetics_percentage, aes(x = genetics, y = dis_perc)) + geom_bar(stat = "identity") + ylim(0,100) + ylab("Percent positive for Salmonella")

## Make the model
genetics.glm <- glm(info$pigpos_fact ~ info$genetics, family = "binomial")
summary(genetics.glm)

ggplot(info, aes(x = as.factor(info$farm), y = genetics, color = pigpos_fact, alpha = 0.6)) + geom_jitter(width = 0.2, height = 0.2)

#### What about interactions? These should be biologically relevant.
## Sex and diet?

sex_diet.glm <- glm(info$pigpos_fact ~ info$sex + info$diet + info$diet:info$sex, family = "binomial")
summary(sex_diet.glm)

trial_diet.glm <- glm(info$pigpos_fact ~ info$trial + info$diet + info$trial:info$diet, family = "binomial")
summary(trial_diet.glm)

genetics_farm.glm <- glm(info$pigpos_fact ~ info$genetics + as.factor(info$farm) + info$genetics:as.factor(info$farm), family = "binomial")
summary(genetics_farm.glm)



## Export
a <- arrangeGrob(sex1, sex2, trial1, trial2, gen1, gen2, farm1, farm2, ncol = 2)
b <- arrangeGrob(diet1, diet2, ncol = 2, nrow = 3)
ggsave(a, file = "~/Dropbox/temp1.pdf", units = "in", height = 11, width = 8.5)
ggsave(b, file = "~/Dropbox/temp2.pdf", units = "in", height = 11, width = 8.5)
