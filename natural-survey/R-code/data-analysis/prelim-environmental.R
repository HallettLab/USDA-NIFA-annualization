# updated 2021-09-24
# preliminary analyses for environmental data
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

# read in environmental data
env <- read.csv('environmental.csv')

# factor all the categorical variables
env[,1:7] <- lapply(env[,1:7], factor)

# experiment was set up as a randomized complete block design with subsampling. 
# 10 blocks across 2 sites (5 blocks per site); each block had a paired annual and perennial stand.
# within each stand, I collected data from 3 subsamples (plots) along a transect.
# Thus, the best model structure should have transect nested within block nested within site as random effects.

#### soil moisture ####
# histogram of soil moisture
hist((env$moist_mean))

# considering soil moisture is proportion data bound by 0 and 1, we should either
# use a logit transformation, or beta regression from glmmTMB
hist(logit(env$moist_mean)) # looks good! logit transform is a lot easier.

# calculate sample means and CIs
sumSE.moist <- summarySE(env, measurevar = "moist_mean", groupvars = 'stand')
sumSE.moist2 <- summarySE(env, measurevar = "moist_mean", groupvars = c('stand','block','site'))

# plot the sample data with blocks as x-axis
pdf(paste0(Rout,"/moist_sampledata.pdf"), height = 2.5, width = 4)
ggplot(data = env) +
  geom_jitter(aes(x = block, y = moist_mean, col = stand, shape = site), height = 0, width = 0.1, alpha = 0.5, size = 1.5) +
  geom_errorbar(data = sumSE.moist2, aes(x = block, ymin = moist_mean - ci, ymax = moist_mean + ci, col = stand), width = 0) +
  geom_point(data = sumSE.moist2, aes(x = block, y = moist_mean, col = stand), size = 3) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank()) +
  labs(x = 'Block', y = 'Soil moisture (VWC)') +
  guides(shape = FALSE) +
  facet_grid(.~ site,scales = 'free')
dev.off()


#### model soil moisture
m.moist1 <- lmer(logit(moist_mean) ~ stand + (1|site/block/transect), data = env)
summary(m.moist1)
AICc(m.moist1)

# getting a singular fit... perhaps we don't need to include transect random effect? 
# Example #1 from this page (https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/models-with-random-effects-blocking-and-pseudoreplication.html)
# is essentially our exact case (except we also have random sites). they do not include a random
# effect to account for the 3 subsamples within a treatment within a block.

# at the same time, a singular fit isn't necessarily despair; it actually still fits a more
# complex random effects model structure than you would be dropping transect completely.

# drop transect completely.
m.moist1 <- lmer(logit(moist_mean) ~ stand + (1|site/block), data = env)
summary(m.moist1)
AICc(m.moist1)

# include the interaction of block and transect nested within site
m.moist1 <- lmer(logit(moist_mean) ~ stand + (1|site/block:transect), data = env)
summary(m.moist1)
# aha! this is the SAME model result as the original (with full nesting), EXCEPT it does not include the
# block:site random term and does not have a singular fit!
AICc(m.moist1)
Anova(m.moist1,type=2,test.statistic = 'F')
ranova(m.moist1)

# calculate modeled means and confidence intervals
emm.moist1 <- data.frame(emmeans(m.moist1, ~ stand, type = 'response'))
emm.moist1

# plot the modeled estimates! 
pdf(paste0(Rout,"/moist_modeled.pdf"), height = 2.5, width = 2.5)
ggplot(data = env) +
  geom_jitter(aes(x = stand, y = moist_mean, col = stand, shape = site), height = 0, width = 0.1, alpha = 0.5, size = 1.5) +
  geom_errorbar(data = emm.moist1, aes(x = stand, ymin = lower.CL, ymax = upper.CL, col = stand), width = 0) +
  geom_point(data = emm.moist1, aes(x = stand, y = response, col = stand), size = 3.5) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = '', y = 'Soil moisture (VWC)') +
  guides(col = FALSE)
dev.off()

#### soil depth ####
# histogram of soil depth
hist((env$depth_mean))

# calculate sample means and CIs
sumSE.depth <- summarySE(env, measurevar = "depth_mean", groupvars = 'stand')
sumSE.depth2 <- summarySE(env, measurevar = "depth_mean", groupvars = c('stand','block','site'))

# plot the sample data with blocks as x-axis
pdf(paste0(Rout,"/depth_sampledata.pdf"), height = 2.5, width = 4)
ggplot(data = env) +
  geom_jitter(aes(x = block, y = depth_mean, col = stand, shape = site), height = 0, width = 0.1, alpha = 0.5, size = 1.5) +
  geom_errorbar(data = sumSE.depth2, aes(x = block, ymin = depth_mean - ci, ymax = depth_mean + ci, col = stand), width = 0) +
  geom_point(data = sumSE.depth2, aes(x = block, y = depth_mean, col = stand), size = 3) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank()) +
  labs(x = 'Block', y = 'Soil depth (cm)') +
  guides(shape = FALSE) +
  facet_grid(.~ site,scales = 'free')
dev.off()


#### model soil depth
m.depth1 <- lmer(depth_mean ~ stand + (1|site/block/transect), data = env)
summary(m.depth1)
AICc(m.depth1)
# same singular fit issue as above.

# drop transect completely.
m.depth1 <- lmer(depth_mean ~ stand + (1|site/block), data = env)
summary(m.depth1)
AICc(m.depth1)
# still have singular fit?

# include the interaction of block and transect nested within site
m.depth1 <- lmer(depth_mean ~ stand + (1|site/block:transect), data = env)
summary(m.depth1)
AICc(m.depth1)
# still have singular fit... but AICc is lowest.
Anova(m.depth1,type=2,test.statistic = 'F')
ranova(m.depth1)

# calculate modeled means and confidence intervals
emm.depth1 <- data.frame(emmeans(m.depth1, ~ stand, type = 'response'))
emm.depth1

# plot the modeled estimates! 
pdf(paste0(Rout,"/depth_modeled.pdf"), height = 2.5, width = 2.5)
ggplot(data = env) +
  geom_jitter(aes(x = stand, y = depth_mean, col = stand, shape = site), height = 0, width = 0.1, alpha = 0.5, size = 1.5) +
  geom_errorbar(data = emm.depth1, aes(x = stand, ymin = lower.CL, ymax = upper.CL, col = stand), width = 0) +
  geom_point(data = emm.depth1, aes(x = stand, y = emmean, col = stand), size = 3.5) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = '', y = 'Soil depth (cm)') +
  guides(col = FALSE)
dev.off()

#### regress soil moisture on depth ####
m.moist2 <- lmer(logit(moist_mean) ~ depth_mean * stand + (1|site/block:transect), data = env)
AICc(m.moist2)
m.moist2 <- lmer(logit(moist_mean) ~ depth_mean + stand + (1|site/block:transect), data = env)
AICc(m.moist2)
summary(m.moist2)

# calculate predicted moisture
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

emm.moist2 <- predict(m.moist2, type = 'response') %>%
  logit2prob(.)
df <- cbind(env,emm.moist2)

# plot moisture as a function of depth and stand
pdf(paste0(Rout,"/moist_depth.pdf"), height = 2.5, width = 4)
ggplot(data = df) +
  geom_point(aes(x = depth_mean, y = moist_mean, col = stand, shape = site), alpha = 0.5, size = 2) +
  # geom_line(aes(x = depth_mean, y = emm.moist2, col = stand)) +
  # geom_point(data = emm.moist1, aes(x = stand, y = response, col = stand), size = 4) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank()) +
  labs(x = 'Soil depth (cm)', y = 'Soil moisture (VWC)') +
  guides(col = FALSE) +
  facet_grid(.~ site)
dev.off()
