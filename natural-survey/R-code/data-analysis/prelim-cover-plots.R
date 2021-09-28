# created 2021-09-27
# preliminary analyses for plot cover/richness
library(vegan)
library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(MuMIn)
library(Rmisc)
library(ggplot2)
library(dplyr)
library(sjPlot)

# set working directory to cleaned data folder
setwd("C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/Natural Site Survey/Reed_USDA-NIFA_NaturalSurvey/Data-cleaned")

# set path to R-output folder within Data-analysis
Rout <- "C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/Natural Site Survey/Reed_USDA-NIFA_NaturalSurvey/R-code/Data-analysis/R-output"

# read in cover data
cov.long <- read.csv('cover_plots-rows.csv') # plots as rows
cov.wide <- read.csv('cover_plots-columns.csv', header = FALSE) # plots as columns

# subset data to subsample plots
cov.long <- subset(cov.long, !is.na(transect))
cov.wide <- cov.wide[,-c(2,6,10,14,18,22,26,30,34,38,
                        42,46,50,54,58,62,66,70,74,78)]

# subsample plots were 1m2. there were 3 plots per stand (meaning 6 plots per block; 60 total across sites).
# within plots, I estimated cover to species level, estimated litter and bareground, and measured litter depth.

# create a matrix of only the species as columns
mtrx <- cov.long[,11:116]

# read in functional group metadata to group species
fxgrps <- read.csv('functional_groups.csv')

# prep a cover dataframe with metadata; factor all the categorical variables
cover <- lapply(cov.long[,1:6], factor) %>%
  as.data.frame(.)

# prep functional group dataframe
fxgrps <- cbind(fxgrps, cov.wide[11:116,2:61])
colnames(fxgrps)[12:71] <- cov.wide[3,2:61]

# make cover measurements numeric
fxgrps[,12:71] <- sapply(fxgrps[,12:71],as.numeric)

#### calculate total and functional group covers ####
cover$covTotal <- colSums(fxgrps[,12:71])
cover$covNative <- colSums(fxgrps[fxgrps$nativity=="Native",12:71])
cover$covIntro <- colSums(fxgrps[fxgrps$nativity=="Introduced",12:71])
cover$covAnnual <- colSums(fxgrps[fxgrps$duration=="Annual",12:71])
cover$covPeren <- colSums(fxgrps[fxgrps$duration=="Perennial",12:71])
cover$covGram <- colSums(fxgrps[fxgrps$form2=="Graminoid",12:71])
cover$covForb <- colSums(fxgrps[fxgrps$form2=="Forb",12:71])
cover$covNatAnn <- colSums(fxgrps[fxgrps$nativity=="Native"&fxgrps$duration=="Annual",12:71])
cover$covNatPer <- colSums(fxgrps[fxgrps$nativity=="Native"&fxgrps$duration=="Perennial",12:71])
cover$covIntAnn <- colSums(fxgrps[fxgrps$nativity=="Introduced"&fxgrps$duration=="Annual",12:71])
cover$covIntPer <- colSums(fxgrps[fxgrps$nativity=="Introduced"&fxgrps$duration=="Perennial",12:71])
cover$covAG <- colSums(fxgrps[fxgrps$form2=="Graminoid"&fxgrps$duration=="Annual",12:71])
cover$covPG <- colSums(fxgrps[fxgrps$form2=="Graminoid"&fxgrps$duration=="Perennial",12:71])
cover$covAF <- colSums(fxgrps[fxgrps$form2=="Forb"&fxgrps$duration=="Annual",12:71])
cover$covPF <- colSums(fxgrps[fxgrps$form2=="Forb"&fxgrps$duration=="Perennial",12:71])
cover$covNAG <- colSums(fxgrps[fxgrps$funcgroup=="NAG",12:71])
cover$covNAF <- colSums(fxgrps[fxgrps$funcgroup=="NAF",12:71])
cover$covNPG <- colSums(fxgrps[fxgrps$funcgroup=="NPG",12:71])
cover$covNPF <- colSums(fxgrps[fxgrps$funcgroup=="NPF",12:71])
cover$covIAG <- colSums(fxgrps[fxgrps$funcgroup=="IAG",12:71])
cover$covIAF <- colSums(fxgrps[fxgrps$funcgroup=="IAF",12:71])
cover$covIPG <- colSums(fxgrps[fxgrps$funcgroup=="IPG",12:71])
cover$covIPF <- colSums(fxgrps[fxgrps$funcgroup=="IPF",12:71])
cover$covWoody <- colSums(fxgrps[fxgrps$funcgroup=="Woody",12:71])
cover$covUnk <- colSums(fxgrps[fxgrps$funcgroup=="Unk",12:71])

#### calculate relative covers ####
relcover<-cover[,1:6]
relcover$relNative<-cover$covNative/cover$covTotal
relcover$relIntro<-cover$covIntro/cover$covTotal
relcover$relAnnual<-cover$covAnnual/cover$covTotal
relcover$relPeren<-cover$covPeren/cover$covTotal
relcover$relGram<-cover$covGram/cover$covTotal
relcover$relForb<-cover$covForb/cover$covTotal
relcover$relNatAnn<-cover$covNatAnn/cover$covTotal
relcover$relNatPer<-cover$covNatPer/cover$covTotal
relcover$relIntAnn<-cover$covIntAnn/cover$covTotal
relcover$relIntPer<-cover$covIntPer/cover$covTotal
relcover$relNAG<-cover$covNAG/cover$covTotal
relcover$relNAF<-cover$covNAF/cover$covTotal
relcover$relNPG<-cover$covNPG/cover$covTotal
relcover$relNPF<-cover$covNPF/cover$covTotal
relcover$relIAG<-cover$covIAG/cover$covTotal
relcover$relIAF<-cover$covIAF/cover$covTotal
relcover$relIPG<-cover$covIPG/cover$covTotal
relcover$relIPF<-cover$covIPF/cover$covTotal
relcover$relWoody<-cover$covWoody/cover$covTotal
relcover$relUnk<-cover$covUnk/cover$covTotal

# set up logit back-transform function for summary table purposes
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

#### Richness ####
cover$rich <- specnumber(mtrx)

# histogram of richness
hist((cover$rich))
hist(sqrt(cover$rich))
hist(log(cover$rich))

# calculate sample means and CIs
sumSE.richness <- summarySE(cover, measurevar = "rich", groupvars = 'stand') # just aggregating to stand
sumSE.richness2 <- summarySE(cover, measurevar = "rich", groupvars = c('stand','block','site')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-plots/richness-plot_sampledata.pdf"), height = 2.5, width = 4)
ggplot(data = cover) +
  geom_jitter(aes(x = block, y = rich, col = stand, shape = site), height = 0, width = 0.1, alpha = 0.5, size = 1.5) +
  geom_errorbar(data = sumSE.richness2, aes(x = block, ymin = rich - ci, ymax = rich + ci, col = stand), width = 0) +
  geom_point(data = sumSE.richness2, aes(x = block, y = rich, col = stand), size = 3) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank()) +
  labs(x = 'Block', y = 'Richness (plots)') +
  guides(shape = FALSE) +
  facet_grid(.~ site,scales = 'free')
dev.off()

# model richness
m.rich1 <- lmer(log(rich) ~ stand + (1|site/block/transect), data = cover)
summary(m.rich1)
AICc(m.rich1)
plot(m.rich1)

# getting a singular fit... perhaps we don't need to include transect random effect? 
# Example #1 from this page (https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/models-with-random-effects-blocking-and-pseudoreplication.html)
# is essentially our exact case (except we also have random sites). they do not include a random
# effect to account for the 3 subsamples within a treatment within a block.

# at the same time, a singular fit isn't necessarily despair; it actually still fits a more
# complex random effects model structure than you would be dropping transect completely.

# drop transect completely.
m.rich1 <- lmer(log(rich) ~ stand + (1|site/block), data = cover)
summary(m.rich1)
AICc(m.rich1)
# still get singular fit. site is the problem.

# include the interaction of block and transect nested within site
m.rich1 <- lmer(log(rich) ~ stand + (1|site/block:transect), data = cover)
summary(m.rich1)
AICc(m.rich1)
Anova(m.rich1,type=3,test.statistic = 'F')
ranova(m.rich1)

# # check site as fixed interaction
# m.rich1 <- lmer(log(rich) ~ stand * site + (1|block:transect), data = cover)
# summary(m.rich1)
# AICc(m.rich1)
# # no significant interaction. just drop site completely?
m.rich1 <- lmer(log(rich) ~ stand + (1|block:transect), data = cover)
summary(m.rich1)
AICc(m.rich1)

# calculate modeled means and confidence intervals
emm.rich1 <- data.frame(emmeans(m.rich1, ~ stand, type = 'response'))
emm.rich1

# make a table of the results
tab_model(m.rich1, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

# plot the modeled estimates! 
pdf(paste0(Rout,"/cover-plots/richness-plot_modeled.pdf"), height = 2.5, width = 2.5)
ggplot(data = cover) +
  geom_jitter(aes(x = stand, y = rich, col = stand, shape = site), height = 0, width = 0.1, alpha = 0.5, size = 1.5) +
  geom_errorbar(data = emm.rich1, aes(x = stand, ymin = lower.CL, ymax = upper.CL, col = stand), width = 0) +
  geom_point(data = emm.rich1, aes(x = stand, y = response, col = stand), size = 3.5) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = '', y = 'Richness (plots)') +
  guides(col = FALSE)
dev.off()

#### Simpson's Diversity ####
cover$diversity <- diversity(mtrx, index = "simpson")

# histogram of diversity
hist((cover$diversity))
hist(logit(cover$diversity))

# calculate sample means and CIs
sumSE.diversity <- summarySE(cover, measurevar = "diversity", groupvars = 'stand') # just aggregating to stand
sumSE.diversity2 <- summarySE(cover, measurevar = "diversity", groupvars = c('stand','block','site')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-plots/diversity-plot_sampledata.pdf"), height = 2.5, width = 4)
ggplot(data = cover) +
  geom_jitter(aes(x = block, y = diversity, col = stand, shape = site), height = 0, width = 0.1, alpha = 0.5, size = 1.5) +
  geom_errorbar(data = sumSE.diversity2, aes(x = block, ymin = diversity - ci, ymax = diversity + ci, col = stand), width = 0) +
  geom_point(data = sumSE.diversity2, aes(x = block, y = diversity, col = stand), size = 3) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank()) +
  labs(x = 'Block', y = "Simpson's Diversity (plots)") +
  guides(shape = FALSE) +
  facet_grid(.~ site,scales = 'free')
dev.off()

# model diversity
m.diversity1 <- lmer(logit(diversity) ~ stand + (1|site/block/transect), data = cover)
summary(m.diversity1)
AICc(m.diversity1)
plot(m.diversity1)

# getting a singular fit... perhaps we don't need to include transect random effect? 
# Example #1 from this page (https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/models-with-random-effects-blocking-and-pseudoreplication.html)
# is essentially our exact case (except we also have random sites). they do not include a random
# effect to account for the 3 subsamples within a treatment within a block.

# at the same time, a singular fit isn't necessarily despair; it actually still fits a more
# complex random effects model structure than you would be dropping transect completely.

# drop transect completely.
m.diversity1 <- lmer(logit(diversity) ~ stand + (1|site/block), data = cover)
summary(m.diversity1)
AICc(m.diversity1)
# still get singular fit. site is the problem.

# include the interaction of block and transect nested within site
m.diversity1 <- lmer(logit(diversity) ~ stand + (1|site/block:transect), data = cover)
summary(m.diversity1)
AICc(m.diversity1)
Anova(m.diversity1,type=2,test.statistic = 'F')
ranova(m.diversity1)

# # check site as fixed interaction
# m.diversity1 <- lmer(logit(diversity) ~ stand * site + (1|block:transect), data = cover)
# summary(m.diversity1)
# AICc(m.diversity1)
# no significant interaction. just drop site completely?
m.diversity1 <- lmer(logit(diversity) ~ stand + (1|block:transect), data = cover)
summary(m.diversity1)
AICc(m.diversity1)

# calculate modeled means and confidence intervals
emm.diversity1 <- data.frame(emmeans(m.diversity1, ~ stand, type = 'response'))
emm.diversity1

# make a table of the results
tab_model(m.diversity1, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

# plot the modeled estimates! 
pdf(paste0(Rout,"/cover-plots/diversity-plot_modeled.pdf"), height = 2.5, width = 2.5)
ggplot(data = cover) +
  geom_jitter(aes(x = stand, y = diversity, col = stand, shape = site), height = 0, width = 0.1, alpha = 0.5, size = 1.5) +
  geom_errorbar(data = emm.diversity1, aes(x = stand, ymin = lower.CL, ymax = upper.CL, col = stand), width = 0) +
  geom_point(data = emm.diversity1, aes(x = stand, y = response, col = stand), size = 3.5) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = '', y = "Simpson's Diversity (plots)") +
  guides(col = FALSE)
dev.off()

#### Evenness ####
cover$dominance <- 1 - cover$diversity # Simpson's dominance
cover$invsimp <- 1/cover$dominance
cover$evenness <- cover$invsimp/cover$rich

# histogram of evenness
hist((cover$evenness))
hist(logit(cover$evenness))

# calculate sample means and CIs
sumSE.evenness <- summarySE(cover, measurevar = "evenness", groupvars = 'stand') # just aggregating to stand
sumSE.evenness2 <- summarySE(cover, measurevar = "evenness", groupvars = c('stand','block','site')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-plots/evenness-plot_sampledata.pdf"), height = 2.5, width = 4)
ggplot(data = cover) +
  geom_jitter(aes(x = block, y = evenness, col = stand, shape = site), height = 0, width = 0.1, alpha = 0.5, size = 1.5) +
  geom_errorbar(data = sumSE.evenness2, aes(x = block, ymin = evenness - ci, ymax = evenness + ci, col = stand), width = 0) +
  geom_point(data = sumSE.evenness2, aes(x = block, y = evenness, col = stand), size = 3) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank()) +
  labs(x = 'Block', y = "Simpson's Evenness (plots)") +
  guides(shape = FALSE) +
  facet_grid(.~ site,scales = 'free')
dev.off()

# model evenness
m.evenness1 <- lmer(logit(evenness) ~ stand + (1|site/block/transect), data = cover)
summary(m.evenness1)
AICc(m.evenness1)
plot(m.evenness1)

# getting a singular fit... perhaps we don't need to include transect random effect? 
# Example #1 from this page (https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/models-with-random-effects-blocking-and-pseudoreplication.html)
# is essentially our exact case (except we also have random sites). they do not include a random
# effect to account for the 3 subsamples within a treatment within a block.

# at the same time, a singular fit isn't necessarily despair; it actually still fits a more
# complex random effects model structure than you would be dropping transect completely.

# drop transect completely.
m.evenness1 <- lmer(logit(evenness) ~ stand + (1|site/block), data = cover)
summary(m.evenness1)
AICc(m.evenness1)
# that fixes it.

# include the interaction of block and transect nested within site
m.evenness1 <- lmer(logit(evenness) ~ stand + (1|site/block:transect), data = cover)
summary(m.evenness1)
AICc(m.evenness1)
Anova(m.evenness1,type=2,test.statistic = 'F')
ranova(m.evenness1)

# # check site as fixed interaction
# m.evenness1 <- lmer(logit(evenness) ~ stand * site + (1|block:transect), data = cover)
# summary(m.evenness1)
# AICc(m.evenness1)

# calculate modeled means and confidence intervals
emm.evenness1 <- data.frame(emmeans(m.evenness1, ~ stand, type = 'response'))
emm.evenness1

# make a table of the results
tab_model(m.evenness1, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

# plot the modeled estimates! 
pdf(paste0(Rout,"/cover-plots/evenness-plot_modeled.pdf"), height = 2.5, width = 2.5)
ggplot(data = cover) +
  geom_jitter(aes(x = stand, y = evenness, col = stand, shape = site), height = 0, width = 0.1, alpha = 0.5, size = 1.5) +
  geom_errorbar(data = emm.evenness1, aes(x = stand, ymin = lower.CL, ymax = upper.CL, col = stand), width = 0) +
  geom_point(data = emm.evenness1, aes(x = stand, y = response, col = stand), size = 3.5) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = '', y = "Simpson's Evenness (plots)") +
  guides(col = FALSE)
dev.off()
