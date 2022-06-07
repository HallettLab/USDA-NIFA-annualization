# Analyzing Pisgah and SEM sites separately
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

# read in biomass data
bms <- read.csv('biomass.csv')

# factor all the categorical variables
bms[,1:7] <- lapply(bms[,1:7], factor)

# subsample plots were 1m2. there were 3 plots per stand (meaning 6 plots per block; 60 total across sites).

# separate data by sites
bms.Pisgah <- subset(bms, site == 'Pisgah')
bms.SEM <- subset(bms, site == 'South Eugene Meadows')

#### biomass ####
# histogram of biomass
hist((bms$biomass_gm2))
hist(sqrt(bms$biomass_gm2))
hist(log(bms$biomass_gm2))

# calculate sample means and CIs
sumSE.biomass <- summarySE(bms, measurevar = "biomass_gm2", groupvars = 'stand') # just aggregating to stand
sumSE.biomass2 <- summarySE(bms, measurevar = "biomass_gm2", groupvars = c('stand','block','site')) # aggregating to stand, block, site

# plot the sample data with blocks as x-axis
# pdf(paste0(Rout,"/biomass/biomass_sampledata.pdf"), height = 2.5, width = 4)
ggplot(data = bms) +
  geom_jitter(aes(x = block, y = biomass_gm2, col = stand, shape = site), height = 0, width = 0.1, alpha = 0.5, size = 1.5) +
  geom_errorbar(data = sumSE.biomass2, aes(x = block, ymin = biomass_gm2 - ci, ymax = biomass_gm2 + ci, col = stand), width = 0) +
  geom_point(data = sumSE.biomass2, aes(x = block, y = biomass_gm2, col = stand), size = 3) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank()) +
  labs(x = 'Block', y = 'Biomass (g/m2)') +
  guides(shape = FALSE) +
  facet_grid(.~ site,scales = 'free')
# dev.off()


#### model biomass
m.biomass1 <- lmer(log(biomass_gm2) ~ stand + (1|site/block/transect), data = bms)
summary(m.biomass1)
AICc(m.biomass1)

# getting a singular fit... perhaps we don't need to include transect random effect? 
# Example #1 from this page (https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/models-with-random-effects-blocking-and-pseudoreplication.html)
# is essentially our exact case (except we also have random sites). they do not include a random
# effect to account for the 3 subsamples within a treatment within a block.

# at the same time, a singular fit isn't necessarily despair; it actually still fits a more
# complex random effects model structure than you would be dropping transect completely.

# drop transect completely.
m.biomass1 <- lmer(log(biomass_gm2) ~ stand + (1|site/block), data = bms)
summary(m.biomass1)
AICc(m.biomass1)
# still get singular fit. site is the problem.

# include the interaction of block and transect nested within site
m.biomass1 <- lmer(log(biomass_gm2) ~ stand + (1|site/block:transect), data = bms)
summary(m.biomass1)
AICc(m.biomass1)
Anova(m.biomass1,type=2,test.statistic = 'F')
ranova(m.biomass1)

# # check site as fixed interaction
# m.biomass1 <- lmer(log(biomass_gm2) ~ stand * site + (1|block:transect), data = bms)
# summary(m.biomass1)
# AICc(m.biomass1)
# no significant interaction. just drop site completely?
m.biomass1 <- lmer(log(biomass_gm2) ~ stand + (1|block:transect), data = bms)
summary(m.biomass1)
AICc(m.biomass1)

# calculate modeled means and confidence intervals
emm.biomass1 <- data.frame(emmeans(m.biomass1, ~ stand, type = 'response'))
emm.biomass1

# make a table of the results
tab_model(m.biomass1, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

# plot the modeled estimates! 
# pdf(paste0(Rout,"/biomass/biomass_modeled.pdf"), height = 2.5, width = 2.5)
ggplot(data = bms) +
  geom_jitter(aes(x = stand, y = biomass_gm2, col = stand, shape = site), height = 0, width = 0.1, alpha = 0.5, size = 1.5) +
  geom_errorbar(data = emm.biomass1, aes(x = stand, ymin = lower.CL, ymax = upper.CL, col = stand), width = 0) +
  geom_point(data = emm.biomass1, aes(x = stand, y = response, col = stand), size = 3.5) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = '', y = 'Biomass (g/m2)') +
  guides(col = FALSE)
# dev.off()

#### regress biomass on soil moisture and depth ####
# bring in environmental data
env <- read.csv('environmental.csv')
bms$moist <- env$moist_mean; bms$depth <- env$depth_mean

m.biomass2 <- lmer(log(biomass_gm2) ~ stand * (depth + moist) + (1|site/block:transect), data = bms)
AICc(m.biomass2)
summary(m.biomass2)

m.biomass2 <- lmer(log(biomass_gm2) ~ stand * (depth + moist) + (1|block:transect), data = bms)
AICc(m.biomass2)
summary(m.biomass2)

# make a table of the results
tab_model(m.biomass2, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

# calculate predicted moisture
nd1 <- expand.grid(depth = seq(min(bms$depth), max(bms$depth), length.out = 20),
                  moist = mean(bms$moist),
                  stand = levels(bms$stand))
nd1$prd <- exp(predict(m.biomass2, newdata = nd1, re.form = NA, type = "response"))

nd2 <- expand.grid(moist = seq(min(bms$moist), max(bms$moist), length.out = 20),
                   depth = mean(bms$depth),
                   stand = levels(bms$stand))
nd2$prd <- exp(predict(m.biomass2, newdata = nd2, re.form = NA, type = "response"))

# plot biomass as a function of depth and stand
# pdf(paste0(Rout,"/biomass/biomass_depth.pdf"), height = 3.5, width = 3.5)
ggplot(data = df) +
  geom_point(aes(x = depth, y = biomass_gm2, col = stand, shape = site), alpha = 0.5, size = 2) +
  geom_line(data = nd1, aes(x = depth, y = prd, col = stand)) +
  # geom_point(data = emm.biomass1, aes(x = stand, y = response, col = stand), size = 4) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank()) +
  labs(x = 'Soil depth (cm)', y = 'Biomass (g/m2)') +
  guides()
# dev.off()

# plot biomass as a function of moisture and stand
# pdf(paste0(Rout,"/biomass/biomass_moist.pdf"), height = 3.5, width = 3.5)
ggplot(data = df) +
  geom_point(aes(x = moist, y = biomass_gm2, col = stand, shape = site), alpha = 0.5, size = 2) +
  geom_line(data = nd2, aes(x = moist, y = prd, col = stand)) +
  # geom_point(data = emm.biomass1, aes(x = stand, y = response, col = stand), size = 4) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank()) +
  labs(x = 'Soil moisture (VWC)', y = 'Biomass (g/m2)') +
  guides()
# dev.off()
