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
bms$block <- factor(bms$block, levels = c('1','2','3','4','5',
                                              '6','7','8','9','10'))
bms$block2 <- factor(rep(c(rep("a",6),rep("b",6),rep("c",6),rep("d",6),rep("e",6)),2))
bms$transect <- factor(bms$transect, levels = c('1a','1b','2a','2b','3a','3b','4a','4b','5a','5b',
                                                    '6a','6b','7a','7b','8a','8b','9a','9b','10a','10b'))

# separate data by sites
bms.Pisgah <- subset(bms, site == 'Pisgah')
bms.SEM <- subset(bms, site == 'South Eugene Meadows')

#### biomass ####
hist(log(bms.Pisgah$biomass_gm2))
hist(log(bms.SEM$biomass_gm2))

# calculate transect means and CIs
sumSE.biomass_gm2 <- summarySE(bms, measurevar = "biomass_gm2", groupvars = c('block','stand','site')) # aggregating to stand, block, site

pdf(paste0(Rout,"/biomass/biomass.pdf"), height = 2.5, width = 4)
ggplot() +
  geom_boxplot(data = bms, aes(x = stand, y = biomass_gm2, fill = stand)) +
  geom_point(data = bms, aes(x = stand, y = biomass_gm2, group = block2, shape = block2), alpha = 0.1, position = position_dodge(0.8), size = 2) +
  geom_line(data = bms, aes(x = stand, y = biomass_gm2, group = transect), alpha = 0.1, position = position_dodge(0.8))+
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_shape_manual(values = rep(c(16,17,15,3,7),2)) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  labs(x = 'Stand', y = "Biomass [g m2] (1 m2 subplots)") +
  facet_grid(.~ site, scales = 'free')
dev.off()

# model biomass
m.biomass.Pisgah <- lmer((biomass_gm2) ~ stand + (1|block), data = bms.Pisgah)
m.biomass.Pisgah <- lmer(log(biomass_gm2) ~ stand + (1|transect), data = bms.Pisgah)
tab_model(m.biomass.Pisgah, linebreak = FALSE, p.val = 'satterthwaite', transform = 'exp')
AICc(m.biomass.Pisgah)

m.biomass.SEM <- lmer(log(biomass_gm2) ~ stand + (1|block), data = bms.SEM)
m.biomass.SEM <- lmer(log(biomass_gm2) ~ stand + (1|transect), data = bms.SEM)
tab_model(m.biomass.SEM, linebreak = FALSE, p.val = 'satterthwaite', transform = 'exp')
AICc(m.biomass.SEM)

# make a table of the results
tab_model(m.biomass.Pisgah, m.biomass.SEM, linebreak = FALSE, p.val = 'satterthwaite', transform = 'exp')

#### regress biomass on soil moisture and depth ####
# bring in environmental data
env <- read.csv('environmental.csv')
bms$moist <- env$moist_mean; bms$depth <- env$depth_mean
bms.Pisgah <- subset(bms, site == 'Pisgah')
bms.SEM <- subset(bms, site == 'South Eugene Meadows')

hist(bms.Pisgah$moist)
hist(bms.SEM$moist)
hist(bms.Pisgah$depth)
hist(bms.SEM$depth)

# plot biomass as a function of depth and stand
# pdf(paste0(Rout,"/biomass/biomass_depth.pdf"), height = 3.5, width = 3.5)
ggplot() +
  geom_point(data = bms, aes(x = depth, y = biomass_gm2, col = stand, shape = block2), size = 2) +
  # geom_smooth(data = bms, aes(x = depth, y = biomass_gm2), method = 'lm') +
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_shape_manual(values = rep(c(16,17,15,3,7),2)) +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  labs(x = 'Soil depth (cm)', y = "Biomass [g m2] (1 m2 subplots)") +
  guides(shape = 'none') +
  facet_grid(.~ site)
# dev.off()

# plot biomass as a function of moisture and stand
# pdf(paste0(Rout,"/biomass/biomass_moist.pdf"), height = 3.5, width = 3.5)
ggplot() +
  geom_point(data = bms, aes(x = moist, y = biomass_gm2, col = stand, shape = block2), size = 2) +
  # geom_smooth(data = bms, aes(x = moist, y = biomass_gm2), method = 'lm') +
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_shape_manual(values = rep(c(16,17,15,3,7),2)) +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  labs(x = 'Soil moisture (vwc)', y = "Biomass [g m2] (1 m2 subplots)") +
  guides(shape = 'none') +
  facet_grid(.~ site)
# dev.off()


m.biomass.Pisgah.depth <- lmer(log(biomass_gm2) ~ stand + depth + (1|block), data = bms.Pisgah)
m.biomass.Pisgah.depth <- lmer(log(biomass_gm2) ~ stand + depth + (1|transect), data = bms.Pisgah)
tab_model(m.biomass.Pisgah.depth, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.biomass.Pisgah.depth)

m.biomass.SEM.depth <- lmer(log(biomass_gm2) ~ depth + (1|block), data = bms.SEM)
m.biomass.SEM.depth <- lmer(log(biomass_gm2) ~ stand + depth + (1|transect), data = bms.SEM)
tab_model(m.biomass.SEM.depth, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.biomass.SEM.depth)

# adding depth does NOT improve model fit over a stand-only model.


m.biomass.Pisgah.moist <- lmer(log(biomass_gm2) ~ stand + moist + (1|block), data = bms.Pisgah)
m.biomass.Pisgah.moist <- lmer(log(biomass_gm2) ~ stand * moist + (1|transect), data = bms.Pisgah)
tab_model(m.biomass.Pisgah.moist, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.biomass.Pisgah.moist)

m.biomass.SEM.moist <- lmer(log(biomass_gm2) ~ stand + moist + (1|block), data = bms.SEM)
m.biomass.SEM.moist <- lmer(log(biomass_gm2) ~ stand * + moist + (1|transect), data = bms.SEM)
tab_model(m.biomass.SEM.moist, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.biomass.SEM.moist)


# try dredging to see which model is best fit
global.Pisgah <- lmer(log(biomass_gm2) ~ stand * (depth + moist) + (1|transect), data = bms.Pisgah, na.action = 'na.fail', REML = FALSE)
global.SEM <- lmer(log(biomass_gm2) ~ stand * (depth + moist) + (1|transect), data = bms.SEM, na.action = 'na.fail', REML = FALSE)
tab_model(global.Pisgah, linebreak = FALSE, p.val = 'satterthwaite', show.aicc = TRUE)

m.biomass.Pisgah.dredge <- dredge(global.Pisgah)
m.biomass.SEM.dredge <- dredge(global.SEM)
