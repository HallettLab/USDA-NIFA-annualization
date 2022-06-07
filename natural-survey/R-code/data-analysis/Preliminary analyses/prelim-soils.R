# created 2021-09-24
# preliminary analyses for soils data
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

# read in soils data
soils <- read.csv('soils.csv')

# factor all the categorical variables
soils[,1:8] <- lapply(soils[,1:8], factor)

# experiment was set up as a randomized complete block design with subsampling. 
# 10 blocks across 2 sites (5 blocks per site); each block had a paired annual and perennial stand.
# within each stand, I collected data from 3 subsamples (plots) along a transect.
# within each plot, I composited 3 soil samples of 0-10cm depth, and 3 samples of 10-30cm depth IF i could harvest from that depth (uncommon).
# Thus, the best model structure should have transect nested within block nested within site as random effects.
# separate samples of different depths
# I also had texture run on a subset of the 0-10cm depth samples (one per stand per block).

# subset texture data
texture <- subset(soils, !is.na(sand_pct))
depth10 <- subset(soils, depth=='0-10cm')
depth30 <- subset(soils, depth=='10-30cm')

#### percent sand ####
hist(texture$sand_pct)

# calculate sample means and CIs
sumSE.sand <- summarySE(texture, measurevar = "sand_pct", groupvars = 'stand') # just aggregating to stand
sumSE.sand2 <- summarySE(texture, measurevar = "sand_pct", groupvars = c('stand','site')) # aggregating to stand, site

# plot the sample data with blocks as x-axis
# pdf(paste0(Rout,"/soils/sand_sampledata.pdf"), height = 2.5, width = 2.5)
ggplot(data = texture) +
  geom_jitter(aes(x = stand, y = sand_pct, col = stand, shape = site), height = 0, width = 0.1, alpha = 0.5, size = 1.5) +
  geom_errorbar(data = sumSE.sand, aes(x = stand, ymin = sand_pct - se, ymax = sand_pct + se, col = stand), width = 0) +
  geom_point(data = sumSE.sand, aes(x = stand, y = sand_pct, col = stand), size = 3.5) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = '', y = 'Percent sand') +
  guides(col = FALSE)
# dev.off()

#### model sand
m.sand1 <- lmer(sand_pct ~ stand + (1|site/block), data = texture)
summary(m.sand1)
AICc(m.sand1)

# calculate modeled means and confidence intervals
emm.sand1 <- data.frame(emmeans(m.sand1, ~ stand, type = 'response'))
emm.sand1

# make a table of the results
tab_model(m.sand1, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

# plot the modeled estimates! 
pdf(paste0(Rout,"/soils/sand_modeled.pdf"), height = 2.5, width = 2.5)
ggplot(data = texture) +
  geom_jitter(aes(x = stand, y = sand_pct, col = stand, shape = site), height = 0, width = 0.1, alpha = 0.5, size = 1.5) +
  geom_errorbar(data = emm.sand1, aes(x = stand, ymin = lower.CL, ymax = upper.CL, col = stand), width = 0) +
  geom_point(data = emm.sand1, aes(x = stand, y = emmean, col = stand), size = 3.5) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = '', y = 'Percent sand') +
  guides(col = FALSE)
dev.off()

#### percent silt ####
hist(texture$silt_pct)

# calculate sample means and CIs
sumSE.silt <- summarySE(texture, measurevar = "silt_pct", groupvars = 'stand') # just aggregating to stand
sumSE.silt2 <- summarySE(texture, measurevar = "silt_pct", groupvars = c('stand','site')) # aggregating to stand, site

# plot the sample data with blocks as x-axis
# pdf(paste0(Rout,"/soils/silt_sampledata.pdf"), height = 2.5, width = 2.5)
ggplot(data = texture) +
  geom_jitter(aes(x = stand, y = silt_pct, col = stand, shape = site), height = 0, width = 0.1, alpha = 0.5, size = 1.5) +
  geom_errorbar(data = sumSE.silt, aes(x = stand, ymin = silt_pct - se, ymax = silt_pct + se, col = stand), width = 0) +
  geom_point(data = sumSE.silt, aes(x = stand, y = silt_pct, col = stand), size = 3.5) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = '', y = 'Percent silt') +
  guides(col = FALSE)
# dev.off()

#### model silt
m.silt1 <- lmer(logit(silt_pct/100) ~ stand + (1|site/block), data = texture)
summary(m.silt1)
AICc(m.silt1)

# calculate modeled means and confidence intervals
emm.silt1 <- data.frame(emmeans(m.silt1, ~ stand, type = 'response'))
emm.silt1

# make a table of the results
tab_model(m.silt1, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

# plot the modeled estimates! 
# pdf(paste0(Rout,"/soils/silt_modeled.pdf"), height = 2.5, width = 2.5)
ggplot(data = texture) +
  geom_jitter(aes(x = stand, y = silt_pct, col = stand, shape = site), height = 0, width = 0.1, alpha = 0.5, size = 1.5) +
  geom_errorbar(data = emm.silt1, aes(x = stand, ymin = lower.CL*100, ymax = upper.CL*100, col = stand), width = 0) +
  geom_point(data = emm.silt1, aes(x = stand, y = response*100, col = stand), size = 3.5) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = '', y = 'Percent silt') +
  guides(col = FALSE)
# dev.off()

#### percent clay ####
hist(texture$clay_pct)

# calculate sample means and CIs
sumSE.clay <- summarySE(texture, measurevar = "clay_pct", groupvars = 'stand') # just aggregating to stand
sumSE.clay2 <- summarySE(texture, measurevar = "clay_pct", groupvars = c('stand','site')) # aggregating to stand, site

# plot the sample data with blocks as x-axis
# pdf(paste0(Rout,"/soils/clay_sampledata.pdf"), height = 2.5, width = 2.5)
ggplot(data = texture) +
  geom_jitter(aes(x = stand, y = clay_pct, col = stand, shape = site), height = 0, width = 0.1, alpha = 0.5, size = 1.5) +
  geom_errorbar(data = sumSE.clay, aes(x = stand, ymin = clay_pct - se, ymax = clay_pct + se, col = stand), width = 0) +
  geom_point(data = sumSE.clay, aes(x = stand, y = clay_pct, col = stand), size = 3.5) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = '', y = 'Percent clay') +
  guides(col = FALSE)
# dev.off()

#### model clay
m.clay1 <- lmer(clay_pct ~ stand + (1|site/block), data = texture)
summary(m.clay1)
AICc(m.clay1)

# calculate modeled means and confidence intervals
emm.clay1 <- data.frame(emmeans(m.clay1, ~ stand, type = 'response'))
emm.clay1

# make a table of the results
tab_model(m.clay1, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

# plot the modeled estimates! 
# pdf(paste0(Rout,"/soils/clay_modeled.pdf"), height = 2.5, width = 2.5)
ggplot(data = texture) +
  geom_jitter(aes(x = stand, y = clay_pct, col = stand, shape = site), height = 0, width = 0.1, alpha = 0.5, size = 1.5) +
  geom_errorbar(data = emm.clay1, aes(x = stand, ymin = lower.CL, ymax = upper.CL, col = stand), width = 0) +
  geom_point(data = emm.clay1, aes(x = stand, y = emmean, col = stand), size = 3.5) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = '', y = 'Percent clay') +
  guides(col = FALSE)
# dev.off()

#### percent carbon, 0-10cm ####
hist((depth10$C_pct))
hist(logit(depth10$C_pct/100))

# calculate sample means and CIs
sumSE.C_pct <- summarySE(depth10, measurevar = "C_pct", groupvars = 'stand', na.rm=T) # just aggregating to stand
sumSE.C_pct2 <- summarySE(depth10, measurevar = "C_pct", groupvars = c('stand','block','site'), na.rm=T) # aggregating to stand, block, site

# plot the sample data with blocks as x-axis
# pdf(paste0(Rout,"/soils/C_pct_sampledata.pdf"), height = 2.5, width = 4)
ggplot(data = depth10) +
  geom_jitter(aes(x = block, y = C_pct, col = stand, shape = site), height = 0, width = 0.1, alpha = 0.5, size = 1.5) +
  geom_errorbar(data = sumSE.C_pct2, aes(x = block, ymin = C_pct - ci, ymax = C_pct + ci, col = stand), width = 0) +
  geom_point(data = sumSE.C_pct2, aes(x = block, y = C_pct, col = stand), size = 3) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank()) +
  labs(x = 'Block', y = 'Percent Carbon') +
  guides(shape = FALSE) +
  facet_grid(.~ site,scales = 'free')
# dev.off()

#### model C_pct
m.C_pct1 <- lmer(logit(C_pct/100) ~ stand + (1|site/block/transect), data = depth10)
summary(m.C_pct1)
AICc(m.C_pct1)

# m.C_pct1 <- lmer(logit(C_pct/100) ~ stand*site + (1|block/transect), data = depth10)
# summary(m.C_pct1)


# calculate modeled means and confidence intervals
emm.C_pct1 <- data.frame(emmeans(m.C_pct1, ~ stand, type = 'response'))
emm.C_pct1

# make a table of the results
tab_model(m.C_pct1, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

# plot the modeled estimates! 
# pdf(paste0(Rout,"/C_pct/C_pct_modeled.pdf"), height = 2.5, width = 2.5)
ggplot(data = depth10) +
  geom_jitter(aes(x = stand, y = C_pct, col = stand, shape = site), height = 0, width = 0.1, alpha = 0.5, size = 1.5) +
  geom_errorbar(data = emm.C_pct1, aes(x = stand, ymin = lower.CL*100, ymax = upper.CL*100, col = stand), width = 0) +
  geom_point(data = emm.C_pct1, aes(x = stand, y = response*100, col = stand), size = 3.5) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = '', y = 'Percent Carbon') +
  guides(col = FALSE)
# dev.off()


#### percent Nitrogen, 0-10cm ####
hist((depth10$N_pct))
hist(logit(depth10$N_pct/100))

# calculate sample means and CIs
sumSE.N_pct <- summarySE(depth10, measurevar = "N_pct", groupvars = 'stand', na.rm=T) # just aggregating to stand
sumSE.N_pct2 <- summarySE(depth10, measurevar = "N_pct", groupvars = c('stand','block','site'), na.rm=T) # aggregating to stand, block, site

# plot the sample data with blocks as x-axis
# pdf(paste0(Rout,"/soils/N_pct_sampledata.pdf"), height = 2.5, width = 4)
ggplot(data = depth10) +
  geom_jitter(aes(x = block, y = N_pct, col = stand, shape = site), height = 0, width = 0.1, alpha = 0.5, size = 1.5) +
  geom_errorbar(data = sumSE.N_pct2, aes(x = block, ymin = N_pct - ci, ymax = N_pct + ci, col = stand), width = 0) +
  geom_point(data = sumSE.N_pct2, aes(x = block, y = N_pct, col = stand), size = 3) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank()) +
  labs(x = 'Block', y = 'Percent Nitrogen') +
  guides(shape = FALSE) +
  facet_grid(.~ site,scales = 'free')
# dev.off()

#### model N_pct
m.N_pct1 <- lmer(logit(N_pct/100) ~ stand + (1|site/block/transect), data = depth10)
summary(m.N_pct1)
AICc(m.N_pct1)

# m.N_pct1 <- lmer(logit(N_pct/100) ~ stand + (1|site), data = sumSE.N_pct2)
# summary(m.N_pct1)
# AICc(m.N_pct1)

# m.N_pct1 <- lm(logit(N_pct/100) ~ stand * site, data = sumSE.N_pct2)
# summary(m.N_pct1)
# AICc(m.N_pct1)

# m.N_pct1 <- lmer(logit(N_pct/100) ~ stand*site + (1|block/transect), data = depth10)
# summary(m.N_pct1)

# calculate modeled means and confidence intervals
emm.N_pct1 <- data.frame(emmeans(m.N_pct1, ~ stand, type = 'response'))
emm.N_pct1

# make a table of the results
tab_model(m.N_pct1, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

# plot the modeled estimates! 
# pdf(paste0(Rout,"/N_pct/N_pct_modeled.pdf"), height = 2.5, width = 2.5)
ggplot(data = depth10) +
  geom_jitter(aes(x = stand, y = N_pct, col = stand, shape = site), height = 0, width = 0.1, alpha = 0.5, size = 1.5) +
  geom_errorbar(data = emm.N_pct1, aes(x = stand, ymin = lower.CL*100, ymax = upper.CL*100, col = stand), width = 0) +
  geom_point(data = emm.N_pct1, aes(x = stand, y = response*100, col = stand), size = 3.5) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = '', y = 'Percent Nitrogen') +
  guides(col = FALSE)
# dev.off()

#### percent OM, 0-10cm ####
hist((depth10$OM_pct))
hist(logit(depth10$OM_pct/100))

# calculate sample means and CIs
sumSE.OM_pct <- summarySE(depth10, measurevar = "OM_pct", groupvars = 'stand', na.rm=T) # just aggregating to stand
sumSE.OM_pct2 <- summarySE(depth10, measurevar = "OM_pct", groupvars = c('stand','block','site'), na.rm=T) # aggregating to stand, block, site

# plot the sample data with blocks as x-axis
# pdf(paste0(Rout,"/soils/OM_pct_sampledata.pdf"), height = 2.5, width = 4)
ggplot(data = depth10) +
  geom_jitter(aes(x = block, y = OM_pct, col = stand, shape = site), height = 0, width = 0.1, alpha = 0.5, size = 1.5) +
  geom_errorbar(data = sumSE.OM_pct2, aes(x = block, ymin = OM_pct - ci, ymax = OM_pct + ci, col = stand), width = 0) +
  geom_point(data = sumSE.OM_pct2, aes(x = block, y = OM_pct, col = stand), size = 3) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank()) +
  labs(x = 'Block', y = 'Percent OM') +
  guides(shape = FALSE) +
  facet_grid(.~ site,scales = 'free')
# dev.off()

#### model OM_pct
m.OM_pct1 <- lmer(logit(OM_pct/100) ~ stand + (1|site/block/transect), data = depth10)
summary(m.OM_pct1)
AICc(m.OM_pct1)

# m.OM_pct1 <- lmer(logit(OM_pct/100) ~ stand + (1|site), data = sumSE.OM_pct2)
# summary(m.OM_pct1)
# AICc(m.OM_pct1)

# m.OM_pct1 <- lm(logit(OM_pct/100) ~ stand * site, data = sumSE.OM_pct2)
# summary(m.OM_pct1)
# AICc(m.OM_pct1)

# m.OM_pct1 <- lmer(logit(OM_pct/100) ~ stand*site + (1|block/transect), data = depth10)
# summary(m.OM_pct1)

# calculate modeled means and confidence intervals
emm.OM_pct1 <- data.frame(emmeans(m.OM_pct1, ~ stand, type = 'response'))
emm.OM_pct1

# make a table of the results
tab_model(m.OM_pct1, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

# plot the modeled estimates! 
# pdf(paste0(Rout,"/OM_pct/OM_pct_modeled.pdf"), height = 2.5, width = 2.5)
ggplot(data = depth10) +
  geom_jitter(aes(x = stand, y = OM_pct, col = stand, shape = site), height = 0, width = 0.1, alpha = 0.5, size = 1.5) +
  geom_errorbar(data = emm.OM_pct1, aes(x = stand, ymin = lower.CL*100, ymax = upper.CL*100, col = stand), width = 0) +
  geom_point(data = emm.OM_pct1, aes(x = stand, y = response*100, col = stand), size = 3.5) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = '', y = 'Percent OM') +
  guides(col = FALSE)
# dev.off()
