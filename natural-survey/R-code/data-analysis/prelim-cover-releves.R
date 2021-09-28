# created 2021-09-27
# preliminary analyses for releve cover/richness
library(vegan)
library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(MuMIn)
library(Rmisc)
library(ggplot2)
library(dplyr)

# set working directory to cleaned data folder
setwd("C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/Natural Site Survey/Reed_USDA-NIFA_NaturalSurvey/Data-cleaned")

# set path to R-output folder within Data-analysis
Rout <- "C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/Natural Site Survey/Reed_USDA-NIFA_NaturalSurvey/R-code/Data-analysis/R-output"

# read in cover data
cov.long <- read.csv('cover_plots-rows.csv') # plots as rows
cov.wide <- read.csv('cover_plots-columns.csv', header = FALSE) # plots as columns

# subset data to large releve plots
cov.long <- subset(cov.long, is.na(transect))
cov.wide <- cov.wide[,c(1,2,6,10,14,18,22,26,30,34,38,
                        42,46,50,54,58,62,66,70,74,78)]

# large releve plots were 100m2. there was 1 releve per stand (meaning 2 releves per block; 20 total across sites).
# within releves, I estimated cover to species level, estimated litter and bareground. did not measure litter depth.
# despite estimating cover to species in releves, I do not think we should analyze cover at the releve level.
# these estimates were very challenging over such a large area, and are going to be more reliable at the subsampled
# plot level. Instead, releve cover should be primarily used for richness.

# create a matrix of only the species as columns
mtrx <- cov.long[,11:116]

# read in functional group metadata to group species
fxgrps <- read.csv('functional_groups.csv')

# prep a cover dataframe with metadata; factor all the categorical variables
cover <- lapply(cov.long[,1:6], factor) %>%
  as.data.frame(.)

# prep functional group dataframe
fxgrps <- cbind(fxgrps, cov.wide[11:116,2:21])
colnames(fxgrps)[12:31] <- cov.wide[3,2:21]

# make cover measurements numeric
fxgrps[,12:31] <- sapply(fxgrps[,12:31],as.numeric)

#### calculate total and functional group covers ####
cover$covTotal <- colSums(fxgrps[,12:31])
cover$covNative <- colSums(fxgrps[fxgrps$nativity=="Native",12:31])
cover$covIntro <- colSums(fxgrps[fxgrps$nativity=="Introduced",12:31])
cover$covAnnual <- colSums(fxgrps[fxgrps$duration=="Annual",12:31])
cover$covPeren <- colSums(fxgrps[fxgrps$duration=="Perennial",12:31])
cover$covGram <- colSums(fxgrps[fxgrps$form2=="Graminoid",12:31])
cover$covForb <- colSums(fxgrps[fxgrps$form2=="Forb",12:31])
cover$covNatAnn <- colSums(fxgrps[fxgrps$nativity=="Native"&fxgrps$duration=="Annual",12:31])
cover$covNatPer <- colSums(fxgrps[fxgrps$nativity=="Native"&fxgrps$duration=="Perennial",12:31])
cover$covIntAnn <- colSums(fxgrps[fxgrps$nativity=="Introduced"&fxgrps$duration=="Annual",12:31])
cover$covIntPer <- colSums(fxgrps[fxgrps$nativity=="Introduced"&fxgrps$duration=="Perennial",12:31])
cover$covAG <- colSums(fxgrps[fxgrps$form2=="Graminoid"&fxgrps$duration=="Annual",12:31])
cover$covPG <- colSums(fxgrps[fxgrps$form2=="Graminoid"&fxgrps$duration=="Perennial",12:31])
cover$covAF <- colSums(fxgrps[fxgrps$form2=="Forb"&fxgrps$duration=="Annual",12:31])
cover$covPF <- colSums(fxgrps[fxgrps$form2=="Forb"&fxgrps$duration=="Perennial",12:31])
cover$covNAG <- colSums(fxgrps[fxgrps$funcgroup=="NAG",12:31])
cover$covNAF <- colSums(fxgrps[fxgrps$funcgroup=="NAF",12:31])
cover$covNPG <- colSums(fxgrps[fxgrps$funcgroup=="NPG",12:31])
cover$covNPF <- colSums(fxgrps[fxgrps$funcgroup=="NPF",12:31])
cover$covIAG <- colSums(fxgrps[fxgrps$funcgroup=="IAG",12:31])
cover$covIAF <- colSums(fxgrps[fxgrps$funcgroup=="IAF",12:31])
cover$covIPG <- colSums(fxgrps[fxgrps$funcgroup=="IPG",12:31])
cover$covIPF <- colSums(fxgrps[fxgrps$funcgroup=="IPF",12:31])
cover$covWoody <- colSums(fxgrps[fxgrps$funcgroup=="Woody",12:31])
cover$covUnk <- colSums(fxgrps[fxgrps$funcgroup=="Unk",12:31])

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



#### Richness ####
cover$rich <- specnumber(mtrx)

# histogram of richness
hist((cover$rich))

# calculate stand means and CIs
sumSE.rich <- summarySE(cover, measurevar = "rich", groupvars = c('site','stand')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-releves/richness-releve_sampledata.pdf"), height = 2.5, width = 4)
ggplot(data = cover) +
  geom_line(aes(x = stand, y = rich, group = block), alpha = 0.2)+
  geom_point(aes(x = stand, y = rich, shape = site), alpha = 0.2, size = 1.5) +
  geom_errorbar(data = sumSE.rich, aes(x = stand, ymin = rich - ci, ymax = rich + ci, col = stand), width = 0) +
  geom_point(data = sumSE.rich, aes(x = stand, y = rich, col = stand), size = 3) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank()) +
  labs(x = 'Stand', y = 'Richness (relevés)') +
  guides(shape = FALSE) +
  facet_grid(.~ site,scales = 'free')
dev.off()

# model richness
m.rich1 <- lmer(rich ~ stand + (1|site/block), data = cover)
summary(m.rich1)
AICc(m.rich1)
plot(m.rich1)
# singular fit issue coming from block... remove from model

m.rich1 <- lmer(rich ~ stand + (1|site), data = cover)
summary(m.rich1)

# calculate modeled means and confidence intervals
emm.rich1 <- data.frame(emmeans(m.rich1, ~ stand, type = 'response'))
emm.rich1

# make a table of the results
tab_model(m.rich1, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

# plot the modeled estimates! 
pdf(paste0(Rout,"/cover-releves/richness-releve_modeled.pdf"), height = 2.5, width = 2.5)
ggplot(data = cover) +
  geom_line(aes(x = stand, y = rich, group = block), alpha = 0.2)+
  geom_point(aes(x = stand, y = rich, shape = site), alpha = 0.2, size = 1.5) +
  geom_errorbar(data = emm.rich1, aes(x = stand, ymin = lower.CL, ymax = upper.CL, col = stand), width = 0) +
  geom_point(data = emm.rich1, aes(x = stand, y = emmean, col = stand), size = 3.5) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = '', y = 'Richness (relevés)') +
  guides(col = FALSE)
dev.off()

#### Simpson's Diversity ####
# (probably not reliable for releves since it depends on cover estimates) 
cover$diversity <- diversity(mtrx, index = "simpson")

# histogram of diversity
hist((cover$diversity))
hist(logit(cover$diversity))

# calculate stand means and CIs
sumSE.diversity <- summarySE(cover, measurevar = "diversity", groupvars = c('site','stand')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-releves/diversity-releve_sampledata.pdf"), height = 2.5, width = 4)
ggplot(data = cover) +
  geom_line(aes(x = stand, y = diversity, group = block), alpha = 0.2)+
  geom_point(aes(x = stand, y = diversity, shape = site), alpha = 0.2, size = 1.5) +
  geom_errorbar(data = sumSE.diversity, aes(x = stand, ymin = diversity - ci, ymax = diversity + ci, col = stand), width = 0) +
  geom_point(data = sumSE.diversity, aes(x = stand, y = diversity, col = stand), size = 3) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank()) +
  labs(x = 'Stand', y = "Simpson's Diversity (relevés)") +
  guides(shape = FALSE) +
  facet_grid(.~ site,scales = 'free')
dev.off()

# model diversity
m.diversity1 <- lmer(logit(diversity) ~ stand + (1|site/block), data = cover)
summary(m.diversity1)
AICc(m.diversity1)
plot(m.diversity1)
# singular fit issue coming from site... remove from model

m.diversity1 <- lmer(logit(diversity) ~ stand + (1|block), data = cover)
summary(m.diversity1)

# calculate modeled means and confidence intervals
emm.diversity1 <- data.frame(emmeans(m.diversity1, ~ stand, type = 'response'))
emm.diversity1

# make a table of the results
tab_model(m.diversity1, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

# plot the modeled estimates! 
pdf(paste0(Rout,"/cover-releves/diversity-releve_modeled.pdf"), height = 2.5, width = 2.5)
ggplot(data = cover) +
  geom_line(aes(x = stand, y = diversity, group = block), alpha = 0.2)+
  geom_point(aes(x = stand, y = diversity, shape = site), alpha = 0.2, size = 1.5) +
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
  labs(x = '', y = 'Diversity (relevés)') +
  guides(col = FALSE)
dev.off()

#### Evenness ####
cover$dominance <- 1 - cover$diversity # Simpson's dominance
cover$invsimp <- 1/cover$dominance
cover$evenness <- cover$invsimp/cover$rich

# histogram of evenness
hist((cover$evenness))
hist(logit(cover$evenness))

# calculate stand means and CIs
sumSE.evenness <- summarySE(cover, measurevar = "evenness", groupvars = c('site','stand')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-releves/evenness-releve_sampledata.pdf"), height = 2.5, width = 4)
ggplot(data = cover) +
  geom_line(aes(x = stand, y = evenness, group = block), alpha = 0.2)+
  geom_point(aes(x = stand, y = evenness, shape = site), alpha = 0.2, size = 1.5) +
  geom_errorbar(data = sumSE.evenness, aes(x = stand, ymin = evenness - ci, ymax = evenness + ci, col = stand), width = 0) +
  geom_point(data = sumSE.evenness, aes(x = stand, y = evenness, col = stand), size = 3) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank()) +
  labs(x = 'Stand', y = "Simpson's evenness (relevés)") +
  guides(shape = FALSE) +
  facet_grid(.~ site,scales = 'free')
dev.off()

# model evenness
m.evenness1 <- lmer(logit(evenness) ~ stand + (1|site/block), data = cover)
summary(m.evenness1)
AICc(m.evenness1)
plot(m.evenness1)

# calculate modeled means and confidence intervals
emm.evenness1 <- data.frame(emmeans(m.evenness1, ~ stand, type = 'response'))
emm.evenness1

# make a table of the results
tab_model(m.evenness1, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

# plot the modeled estimates! 
pdf(paste0(Rout,"/cover-releves/evenness-releve_modeled.pdf"), height = 2.5, width = 2.5)
ggplot(data = cover) +
  geom_line(aes(x = stand, y = evenness, group = block), alpha = 0.2)+
  geom_point(aes(x = stand, y = evenness, shape = site), alpha = 0.2, size = 1.5) +
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
  labs(x = '', y = 'Evenness (relevés)') +
  guides(col = FALSE)
dev.off()
