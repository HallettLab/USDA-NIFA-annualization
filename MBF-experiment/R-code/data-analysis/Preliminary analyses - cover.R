#### 6/10/2022: Preliminary analyses MBF experiment
library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(MuMIn)
library(Rmisc)
library(ggplot2)
library(dplyr)
library(sjPlot)
library(tidyr)
options(contrasts = c("contr.sum", "contr.poly"))

# set working directory to cleaned data folder
setwd("C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/MBF Experiment/Reed_USDA-NIFA_MBF Experiment/Data-cleaned")

# set path to R-output folder within Data-analysis
Rout <- "C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/MBF Experiment/Reed_USDA-NIFA_MBF Experiment/R-code/Data-analysis/R-output"

# read in processed cover data
cover <- read.csv('cover-processed.csv')
cover.avg <- read.csv('cover_avg-processed.csv')

# factor all the categorical variables
cover[,1:9] <- lapply(cover[,1:9], factor)
cover.avg[,1:7] <- lapply(cover.avg[,1:7], factor)

levels(cover$pasture)
levels(cover$pasture) <- c('Ryegrass', 'Perennial forage mix')
levels(cover$month)
levels(cover$month) <- c('Nov','Jan','Feb','Mar','Apr','May')

levels(cover.avg$pasture)
levels(cover.avg$pasture) <- c('Ryegrass', 'Perennial forage mix')

cover.seeded <- subset(cover, seeded == 'Seeded')
# throw away november sampling since we measured cover ~1 week after seeding
cover.seeded <- subset(cover.seeded, !month == 'Nov')

cover.avg.seeded <- subset(cover.avg, seeded == 'Seeded')

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}


#### Total native cover ####
# calculate means and SEs
sumSE.nativecov <- summarySE(cover.seeded, measurevar = "Seeded.native", groupvars = c('month','pasture','grazing','seeded'))

png(paste0(Rout,"/natcov.png"),height = 3.5, width = 6, units = "in", res = 600)
ggplot(cover.seeded, aes(x = month, y = Seeded.native, fill = grazing)) +
  geom_boxplot() +
  # geom_vline(xintercept = )
  labs(x = 'Month', y = 'Native (%) cover', fill = '') +
  # scale_y_log10() +
  theme_bw() +
  theme(legend.margin = margin(0,0,0,0),
        legend.position = 'bottom') +
  # legend.title = element_blank())
  facet_grid(. ~ pasture)
dev.off()

# ggplot() +
#   geom_point(data = cover.seeded, aes(x = month, y = Seeded.native, group = block, shape = grazing, col = grazing), alpha = 0.2, position = position_dodge(0.8), size = 2) +
#   geom_errorbar(data = sumSE.nativecov, aes(x = month, ymin = Seeded.native - se, ymax = Seeded.native + se, col = grazing), width = 0) +
#   geom_point(data = sumSE.nativecov, aes(x = month, y = Seeded.native, col = grazing, shape = grazing), size = 4) +
#   # geom_line(data = cover.seeded, aes(x = month, y = Seeded.native, group = block), alpha = 0.1, position = position_dodge(0.8))+
#   labs(x = 'Month', y = 'Total native cover', col = '', shape = '') +
#   theme_bw() +
#   theme(legend.margin = margin(0,0,0,0)) +
#   # legend.title = element_blank())
#   facet_grid(seeded ~ pasture)

# run a model
cover.seeded$Seeded.native.logit <- logit(cover.seeded$Seeded.native, adjust = 0.01)
hist(cover.seeded$Seeded.native.logit)
hist(cover.seeded$Seeded.native)

m.natcov <- lmer(Seeded.native.logit ~ grazing * pasture * month + (1|block) + (1|plot), data = cover.seeded, REML = FALSE)
m.natcov <- lmer(Seeded.native.logit ~ grazing * pasture * month + (1|plot), data = cover.seeded, REML = FALSE)
AICc(m.natcov)
summary(m.natcov)
Anova(m.natcov, type = 3)

# january
m.natcov.jan <- lm(Seeded.native.logit ~ grazing * pasture, data = cover.seeded[cover.seeded$month == 'Jan',])
Anova(m.natcov.jan, type = 3)
summary(m.natcov.jan)
emm.natcov.jan <- emmeans(m.natcov.jan, pairwise ~ grazing|pasture)

# february
m.natcov.feb <- lm(Seeded.native.logit ~ grazing * pasture, data = cover.seeded[cover.seeded$month == 'Feb',])
Anova(m.natcov.feb, type = 3)
summary(m.natcov.feb)
emm.natcov.feb <- emmeans(m.natcov.feb, pairwise ~ grazing|pasture)

# march
m.natcov.mar <- lm(Seeded.native.logit ~ grazing * pasture, data = cover.seeded[cover.seeded$month == 'Mar',])
Anova(m.natcov.mar, type = 3)
summary(m.natcov.mar)
emm.natcov.mar <- emmeans(m.natcov.mar, pairwise ~ grazing)

# april
m.natcov.apr <- lm(Seeded.native.logit ~ grazing * pasture, data = cover.seeded[cover.seeded$month == 'Apr',])
Anova(m.natcov.apr, type = 3)
summary(m.natcov.apr)
emm.natcov.apr <- emmeans(m.natcov.apr, pairwise ~ grazing)

# may
m.natcov.may <- lm(Seeded.native.logit ~ grazing * pasture, data = cover.seeded[cover.seeded$month == 'May',])
Anova(m.natcov.may, type = 3)
summary(m.natcov.may)
emm.natcov.may <- emmeans(m.natcov.may, pairwise ~ grazing)

# #### Total native cover: average two plots together for each block ####
# sumSE.nativecov2 <- summarySE(cover.seeded, measurevar = "Seeded.native", groupvars = c('month','block','pasture','grazing','seeded'))
# hist(logit(sumSE.nativecov2$Seeded.native, adjust = 0.01))
# sumSE.nativecov2$Seeded.native.logit <- logit(sumSE.nativecov2$Seeded.native, adjust = 0.01)
# 
# ggplot(sumSE.nativecov2, aes(x = month, y = Seeded.native, fill = grazing)) +
#   geom_boxplot() +
#   labs(x = 'Month', y = 'Total native cover', fill = '') +
#   theme_bw() +
#   theme(legend.margin = margin(0,0,0,0)) +
#   # legend.title = element_blank())
#   facet_grid(seeded ~ pasture)
# 
# m.natcov <- lmer(Seeded.native.logit ~ grazing * pasture * month + (1|block), data = sumSE.nativecov2, REML = FALSE)
# AICc(m.natcov)
# summary(m.natcov)
# Anova(m.natcov)
# 
# # january
# m.natcov.jan <- lm(Seeded.native.logit ~ grazing * pasture, data = sumSE.nativecov2[sumSE.nativecov2$month == 'Jan',])
# summary(m.natcov.jan)
# emm.natcov.jan <- emmeans(m.natcov.jan, pairwise ~ grazing|pasture)
# 
# # february
# m.natcov.feb <- lm(Seeded.native.logit ~ grazing * pasture, data = sumSE.nativecov2[sumSE.nativecov2$month == 'Feb',])
# summary(m.natcov.feb)
# emm.natcov.feb <- emmeans(m.natcov.feb, pairwise ~ grazing|pasture)
# 
# # march
# m.natcov.mar <- lm(Seeded.native.logit ~ grazing * pasture, data = sumSE.nativecov2[sumSE.nativecov2$month == 'Mar',])
# summary(m.natcov.mar)
# emm.natcov.mar <- emmeans(m.natcov.mar, pairwise ~ grazing|pasture)
# 
# # april
# m.natcov.apr <- lm(Seeded.native.logit ~ grazing * pasture, data = sumSE.nativecov2[sumSE.nativecov2$month == 'Apr',])
# summary(m.natcov.apr)
# emm.natcov.apr <- emmeans(m.natcov.apr, pairwise ~ grazing|pasture)
# 
# # may
# m.natcov.may <- lm(Seeded.native.logit ~ grazing * pasture, data = sumSE.nativecov2[sumSE.nativecov2$month == 'May',])
# summary(m.natcov.may)
# emm.natcov.may <- emmeans(m.natcov.may, pairwise ~ grazing|pasture)


#### Native richness ####
# calculate means and SEs
sumSE.natrich <- summarySE(cover.seeded, measurevar = "natrich", groupvars = c('month','pasture','grazing','seeded'))

# png(paste0(Rout,"/natrich.png"),height = 3.5, width = 6, units = "in", res = 600)
ggplot(cover.seeded, aes(x = month, y = natrich, fill = grazing)) +
  geom_boxplot() +
  labs(x = 'Month', y = 'Native richness', fill = '') +
  theme_bw() +
  theme(legend.margin = margin(0,0,0,0),
        legend.position = 'bottom') +
  # legend.title = element_blank())
  facet_grid(seeded ~ pasture)
dev.off()
# ggplot() +
#   geom_point(data = cover.seeded, aes(x = month, y = natrich, group = block, shape = grazing, col = grazing), alpha = 0.2, position = position_dodge(0.8), size = 2) +
#   geom_errorbar(data = sumSE.natrich, aes(x = month, ymin = natrich - se, ymax = natrich + se, col = grazing), width = 0) +
#   geom_point(data = sumSE.natrich, aes(x = month, y = natrich, col = grazing, shape = grazing), size = 4) +
#   # geom_line(data = cover.seeded, aes(x = month, y = natrich, group = block), alpha = 0.1, position = position_dodge(0.8))+
#   labs(x = 'Month', y = 'Native richness', col = '', shape = '') +
#   theme_bw() +
#   theme(legend.margin = margin(0,0,0,0)) +
#   # legend.title = element_blank())
#   facet_grid(seeded ~ pasture)

# run a model
hist(cover.seeded$natrich)

m.natrich <- glmer(natrich ~ grazing * pasture * month + (1|block) + (1|plot),family = 'poisson', data = cover.seeded)
m.natrich <- glmer(natrich ~ grazing * pasture * month + (1|plot), data = cover.seeded,family = 'poisson')
AICc(m.natrich)
summary(m.natrich)
Anova(m.natrich, type = 3)
tab_model(m.natrich)

# january
m.natrich.jan <- glm(natrich ~ grazing * pasture, data = cover.seeded[cover.seeded$month == 'Jan',],family = 'poisson')
summary(m.natrich.jan)
emm.natrich.jan <- emmeans(m.natrich.jan, pairwise ~ grazing|pasture)

# february
m.natrich.feb <- glm(natrich ~ grazing * pasture, data = cover.seeded[cover.seeded$month == 'Feb',],family = 'poisson')
summary(m.natrich.feb)
emm.natrich.feb <- emmeans(m.natrich.feb, pairwise ~ grazing|pasture)

# march
m.natrich.mar <- glm(natrich ~ grazing * pasture, data = cover.seeded[cover.seeded$month == 'Mar',],family = 'poisson')
summary(m.natrich.mar)
emm.natrich.mar <- emmeans(m.natrich.mar, pairwise ~ grazing|pasture)

# april
m.natrich.apr <- glm(natrich ~ grazing * pasture, data = cover.seeded[cover.seeded$month == 'Apr',],family = 'poisson')
summary(m.natrich.apr)
emm.natrich.apr <- emmeans(m.natrich.apr, pairwise ~ grazing|pasture)

# may
m.natrich.may <- glm(natrich ~ grazing * pasture, data = cover.seeded[cover.seeded$month == 'May',],family = 'poisson')
summary(m.natrich.may)
emm.natrich.may <- emmeans(m.natrich.may, pairwise ~ grazing|pasture)

#### Native cover (averaged across time) ####
hist(logit(cover.avg.seeded$Seeded.native, adjust = 0.01))
cover.avg.seeded$Seeded.native.logit <- logit(cover.avg.seeded$Seeded.native, adjust = 0.01)

# calculate means and SEs
sumSE.nativecov.avg <- summarySE(cover.avg.seeded, measurevar = "Seeded.native", groupvars = c('pasture','grazing','seeded'))

ggplot(cover.avg.seeded, aes(x = grazing, y = Seeded.native, fill = grazing)) +
  geom_boxplot() +
  labs(x = '', y = 'Avg. total native cover', fill = '') +
  theme_bw() +
  theme(legend.margin = margin(0,0,0,0)) +
  # legend.title = element_blank())
  facet_grid(seeded ~ pasture)

ggplot() +
  geom_point(data = cover.avg.seeded, aes(x = grazing, y = Seeded.native, group = block, shape = grazing, col = grazing), alpha = 0.2, position = position_dodge(0.8), size = 2) +
  geom_errorbar(data = sumSE.nativecov.avg, aes(x = grazing, ymin = Seeded.native - se, ymax = Seeded.native + se, col = grazing), width = 0) +
  geom_point(data = sumSE.nativecov.avg, aes(x = grazing, y = Seeded.native, col = grazing, shape = grazing), size = 4) +
  # geom_line(data = cover.seeded, aes(x = month, y = Seeded.native, group = block), alpha = 0.1, position = position_dodge(0.8))+
  labs(x = '', y = 'Avg. total native cover', col = '', shape = '') +
  theme_bw() +
  theme(legend.margin = margin(0,0,0,0)) +
  # legend.title = element_blank())
  facet_grid(seeded ~ pasture)

m.natcov.avg <- lmer(Seeded.native.logit ~ grazing * pasture + (1|block), data = cover.avg.seeded, REML = FALSE)
m.natcov.avg <- lm(Seeded.native.logit ~ grazing * pasture, data = cover.avg.seeded)
AICc(m.natcov.avg)
summary(m.natcov.avg)
Anova(m.natcov.avg, type = 2)

#### Native richness (cumulative across time) ####
# calculate means and SEs
sumSE.natrich.avg <- summarySE(cover.avg.seeded, measurevar = "natrich", groupvars = c('pasture','grazing','seeded'))

png(paste0(Rout,"/natrich.png"),height = 3.5, width = 6, units = "in", res = 600)
ggplot(cover.avg.seeded, aes(x = grazing, y = natrich, fill = grazing)) +
  geom_boxplot() +
  labs(x = '', y = 'Cumulative native richness', fill = '') +
  theme_bw() +
  theme(legend.margin = margin(0,0,0,0),
        legend.position = 'bottom') +
  # legend.title = element_blank())
  facet_grid(. ~ pasture)
dev.off()

ggplot() +
  geom_point(data = cover.avg.seeded, aes(x = grazing, y = natrich, group = block, shape = grazing, col = grazing), alpha = 0.2, position = position_dodge(0.8), size = 2) +
  geom_errorbar(data = sumSE.natrich.avg, aes(x = grazing, ymin = natrich - se, ymax = natrich + se, col = grazing), width = 0) +
  geom_point(data = sumSE.natrich.avg, aes(x = grazing, y = natrich, col = grazing, shape = grazing), size = 4) +
  # geom_line(data = cover.seeded, aes(x = month, y = natrich, group = block), alpha = 0.1, position = position_dodge(0.8))+
  labs(x = '', y = 'Cumulative native richness', col = '', shape = '') +
  theme_bw() +
  theme(legend.margin = margin(0,0,0,0)) +
  # legend.title = element_blank())
  facet_grid(seeded ~ pasture)

m.natrich.avg <- glmer(natrich ~ grazing * pasture + (1|block), data = cover.avg.seeded, family = 'poisson')
m.natrich.avg <- glm(natrich ~ grazing * pasture, data = cover.avg.seeded, family = 'poisson')
AICc(m.natrich.avg)
summary(m.natrich.avg)
Anova(m.natrich.avg, type = 3)

#### native cover ~ bareground ####
ggplot() +
  geom_point(data = cover.avg.seeded, aes(x = Bareground_November, y = Seeded.native, shape = pasture, col = grazing), position = position_dodge(0.8), size = 2) +
  # geom_smooth(data = cover.avg.seeded, aes(x = Bareground_November, y = Seeded.native, shape = grazing, col = grazing), method = 'lm') +
  labs(x = 'November (%) bareground', y = 'Avg. total native cover', col = '', shape = 'Pasture type') +
  theme_bw() +
  theme(legend.margin = margin(0,0,0,0))
# legend.title = element_blank())

m.natcov.avg_bare <- lmer(Seeded.native.logit ~ grazing * pasture * Bareground_November + (1|block), data = cover.avg.seeded, REML = FALSE)
m.natcov.avg_bare <- lm(Seeded.native.logit ~ grazing * pasture * Bareground_November, data = cover.avg.seeded)
AICc(m.natcov.avg_bare)
tab_model(m.natcov.avg_bare)
Anova(m.natcov.avg_bare, type = 2)
m.natcov.avg_bare <- lm(Seeded.native.logit ~ Bareground_November, data = cover.avg.seeded)

# plot predictions
# nd = expand.grid(grazing = levels(cover.avg.seeded$grazing),
#                  pasture = levels(cover.avg.seeded$pasture),
#                  Bareground_November = seq(min(cover.avg.seeded$Bareground_November),
#                                            max(cover.avg.seeded$Bareground_November), length.out = 20))
nd = expand.grid(Bareground_November = seq(min(cover.avg.seeded$Bareground_November),
                                           max(cover.avg.seeded$Bareground_November), length.out = 20))
nd$pred <- predict(m.natcov.avg_bare, newdata = nd)
nd$pred.back <- logit2prob(nd$pred) *100

png(paste0(Rout,"/natcov_bareground.png"),height = 3, width = 4.5, units = "in", res = 600)
ggplot() +
  geom_point(data = cover.avg.seeded, aes(x = Bareground_November, y = Seeded.native, shape = pasture, col = grazing), size = 2, alpha = 0.5) +
  geom_line(data = nd, aes(x = Bareground_November, y = pred.back)) +
  labs(x = 'November (%) bareground', y = 'Avg. native (%) cover', col = '', shape = 'Pasture type') +
  theme_bw() +
  theme(legend.margin = margin(0,0,0,0))
dev.off()


#### native cover ~ grass ####
ggplot() +
  geom_point(data = cover.avg.seeded, aes(x = Forage.grass, y = Seeded.native, shape = pasture, col = grazing), position = position_dodge(0.8), size = 2) +
  # geom_smooth(data = cover.avg.seeded, aes(x = Bareground_November, y = Seeded.native, shape = grazing, col = grazing), method = 'lm') +
  labs(x = 'Grass (%) cover', y = 'Avg. total native cover', col = '', shape = '') +
  theme_bw() +
  theme(legend.margin = margin(0,0,0,0)) +
  # legend.title = element_blank())
  facet_grid(seeded ~ .)

m.natcov.avg_grass <- lmer(Seeded.native.logit ~ grazing * pasture * Forage.grass + (1|block), data = cover.avg.seeded, REML = FALSE)
m.natcov.avg_grass <- lm(Seeded.native.logit ~ grazing * pasture * Forage.grass, data = cover.avg.seeded)
AICc(m.natcov.avg_grass)
summary(m.natcov.avg_grass)
tab_model(m.natcov.avg_grass)
Anova(m.natcov.avg_grass, type = 2)
m.natcov.avg_grass <- lm(Seeded.native.logit ~ Forage.grass, data = cover.avg.seeded)

nd = expand.grid(Forage.grass = seq(min(cover.avg.seeded$Forage.grass),
                                    max(cover.avg.seeded$Forage.grass), length.out = 20))
nd$pred <- predict(m.natcov.avg_grass, newdata = nd)
nd$pred.back <- logit2prob(nd$pred) *100

png(paste0(Rout,"/natcov_grass.png"),height = 3, width = 4.5, units = "in", res = 600)
ggplot() +
  geom_point(data = cover.avg.seeded, aes(x = Forage.grass, y = Seeded.native, shape = pasture, col = grazing), size = 2, alpha = 0.5) +
  geom_line(data = nd, aes(x = Forage.grass, y = pred.back)) +
  labs(x = 'Avg. grass (%) cover', y = 'Avg. native (%) cover', col = '', shape = 'Pasture type') +
  theme_bw() +
  theme(legend.margin = margin(0,0,0,0))
dev.off()










#### native rich ~ bareground ####
ggplot() +
  geom_point(data = cover.avg.seeded, aes(x = Bareground_November, y = natrich, shape = pasture, col = grazing), position = position_dodge(0.8), size = 2) +
  # geom_smooth(data = cover.avg.seeded, aes(x = Bareground_November, y = Seeded.native, shape = grazing, col = grazing), method = 'lm') +
  labs(x = 'November (%) bareground', y = 'Cumulative native richness', col = '', shape = '') +
  theme_bw() +
  theme(legend.margin = margin(0,0,0,0)) +
  # legend.title = element_blank())
  facet_grid(seeded ~ .)

m.natrich.avg_bare <- glmer(natrich ~ grazing * pasture * Bareground_November + (1|block), data = cover.avg.seeded, family = 'poisson')
m.natrich.avg_bare <- glm(natrich ~ grazing * pasture * Bareground_November, data = cover.avg.seeded, family = 'poisson')
AICc(m.natrich.avg_bare)
summary(m.natrich.avg_bare)
tab_model(m.natrich.avg_bare)
Anova(m.natrich.avg_bare, type = 2)
m.natrich.avg_bare <- glm(natrich ~ Bareground_November, data = cover.avg.seeded, family = 'poisson')

# plot predictions
# nd = expand.grid(grazing = levels(richer.avg.seeded$grazing),
#                  pasture = levels(richer.avg.seeded$pasture),
#                  Bareground_November = seq(min(richer.avg.seeded$Bareground_November),
#                                            max(richer.avg.seeded$Bareground_November), length.out = 20))
nd = expand.grid(Bareground_November = seq(min(cover.avg.seeded$Bareground_November),
                                           max(cover.avg.seeded$Bareground_November), length.out = 20))
nd$pred <- exp(predict(m.natrich.avg_bare, newdata = nd))
# nd$pred.back <- logit2prob(nd$pred) *100

png(paste0(Rout,"/natrich_bareground2.png"),height = 3, width = 4.5, units = "in", res = 600)
ggplot() +
  geom_point(data = cover.avg.seeded, aes(x = Bareground_November, y = natrich, shape = pasture, col = grazing), alpha = 0.5, size = 2) +
  geom_line(data = nd, aes(x = Bareground_November, y = pred)) +
  labs(x = 'November (%) bareground', y = 'Cumulative native richness', col = '', shape = 'Pasture type') +
  theme_bw() +
  theme(legend.margin = margin(0,0,0,0))
dev.off()

#### native rich ~ grass ####
ggplot() +
  geom_point(data = cover.avg.seeded, aes(x = Forage.grass, y = natrich, shape = pasture, col = grazing), size = 2) +
  # geom_smooth(data = cover.avg.seeded, aes(x = Forage.grass, y = natrich, shape = grazing, col = grazing), method = 'lm') +
  labs(x = 'Grass (%) cover', y = 'Cumulative native richness', col = '', shape = '') +
  theme_bw() +
  theme(legend.margin = margin(0,0,0,0)) +
  # legend.title = element_blank())
  facet_grid(seeded ~ .)

m.natrich.avg_grass <- glmer(natrich ~ grazing * pasture * Forage.grass + (1|block), data = cover.avg.seeded, family = 'poisson')
m.natrich.avg_grass <- glm(natrich ~ grazing * pasture * Forage.grass, data = cover.avg.seeded, family = 'poisson')
AICc(m.natrich.avg_grass)
summary(m.natrich.avg_grass)
tab_model(m.natrich.avg_grass)
Anova(m.natrich.avg_grass, type = 2)
m.natrich.avg_grass <- glm(natrich ~ Forage.grass, data = cover.avg.seeded, family = 'poisson')

# plot predictions
# nd = expand.grid(grazing = levels(richer.avg.seeded$grazing),
#                  pasture = levels(richer.avg.seeded$pasture),
#                  Forage.grass = seq(min(richer.avg.seeded$Forage.grass),
#                                            max(richer.avg.seeded$Forage.grass), length.out = 20))
nd = expand.grid(Forage.grass = seq(min(cover.avg.seeded$Forage.grass),
                                           max(cover.avg.seeded$Forage.grass), length.out = 20))
nd$pred <- exp(predict(m.natrich.avg_grass, newdata = nd))
# nd$pred.back <- logit2prob(nd$pred) *100

png(paste0(Rout,"/natrich_grass2.png"),height = 3, width = 4.5, units = "in", res = 600)
ggplot() +
  geom_point(data = cover.avg.seeded, aes(x = Forage.grass, y = natrich, shape = pasture, col = grazing), alpha = 0.5, size = 2) +
  geom_line(data = nd, aes(x = Forage.grass, y = pred)) +
  labs(x = 'Avg. grass (%) cover', y = 'Cumulative native richness', col = '', shape = 'Pasture type') +
  theme_bw() +
  theme(legend.margin = margin(0,0,0,0))
dev.off()

#### chronology of individual native species ####
# natives <- cover.seeded[,-26] # get rid of native seedling unknown
natives <- gather(cover.seeded, species, cover, Achillea.millefolium:Prunella.vulgaris, factor_key=TRUE)
natives$species <- gsub("\\."," ",as.character(natives$species))
natives <- subset(natives, !month == 'Nov')

natives$species <- factor(natives$species, levels = c('Acmispon americanus','Clarkia purpurea','Collinsia grandiflora','Collomia grandiflora','Gilia capitata','Plectritis congesta',
                                                      'Achillea millefolium','Eriophyllum lanatum','Lomatium utriculatum','Prunella vulgaris','Native seedling  unknown'))
levels(natives$species) <- c('Acmispon','Clarkia','Collinsia','Collomia','Gilia','Plectritis',
                                                      'Achillea','Eriophyllum','Lomatium','Prunella','UNKNOWN')
natives[is.na(natives$species),]$species <- 'UNKNOWN'

# calculate means and SEs
sumSE.species <- summarySE(natives, measurevar = "cover", groupvars = c('month','pasture','grazing','seeded','species'))

ggplot(natives, aes(x = month, y = cover, fill = species)) +
  geom_boxplot() +
  # geom_vline(xintercept = )
  labs(x = 'Month', y = 'Native (%) cover', fill = '') +
  # scale_y_log10() +
  theme_bw() +
  theme(legend.margin = margin(0,0,0,0),
        legend.position = 'bottom') +
  # legend.title = element_blank())
  facet_grid(pasture ~ grazing)
dev.off()

ggplot() +
  # geom_point(data = natives, aes(x = month, y = cover, group = block, shape = grazing, col = species), alpha = 0.2, position = position_dodge(0.8), size = 2) +
  # geom_errorbar(data = sumSE.species, aes(x = month, ymin = cover - se, ymax = cover + se, col = species), width = 0) +
  geom_point(data = sumSE.species, aes(x = month, y = cover, col = species)) +
  geom_line(data = sumSE.species, aes(x = month, y = cover, col = species, group = species), alpha = 0.5)+
  scale_color_manual(values = c('forestgreen','deeppink','blue3','tan1','lightskyblue','lightpink','chartreuse3','gold','darkred','darkorchid','darkgrey')) +
  labs(x = 'Month', y = 'Mean (%) cover', col = '', shape = '') +
  theme_bw() +
  theme(legend.margin = margin(0,0,0,0)) +
  # legend.title = element_blank())
  facet_grid(pasture ~ grazing)

#### try presence/absence of natives through time ####
natives$PA <- 0
natives[natives$cover > 0,]$PA <- 1

natives.PA <- aggregate(PA ~ month + pasture + grazing + species, data = natives, FUN = sum)

png(paste0(Rout,"/species_presence.png"),height = 4.5, width = 6, units = "in", res = 600)
ggplot() +
  # geom_point(data = natives, aes(x = month, y = cover, group = block, shape = grazing, col = species), alpha = 0.2, position = position_dodge(0.8), size = 2) +
  # geom_errorbar(data = sumSE.species, aes(x = month, ymin = cover - se, ymax = cover + se, col = species), width = 0) +
  geom_point(data = natives.PA, aes(x = month, y = PA, col = species)) +
  geom_line(data = natives.PA, aes(x = month, y = PA, col = species, group = species), alpha = 0.5)+
  scale_color_manual(values = c('forestgreen','deeppink','blue3','tan1','lightskyblue','lightpink','chartreuse3','gold','darkred','darkorchid','darkgrey')) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10,12)) +
  labs(x = 'Month', y = 'Number of plots present', col = '', shape = '') +
  theme_bw() +
  theme(legend.margin = margin(0,0,0,0)) +
  facet_grid(pasture ~ grazing)
dev.off()


#### chronology of forage mix ####
forage <- gather(cover, species, cover, Brassica.spp.:x.Festulolium, factor_key=TRUE)
# forage$species <- gsub("\\."," ",as.character(forage$species))
# forage$species <- gsub("spp","spp.",as.character(forage$species))
forage <- subset(forage, !species == 'x.Festulolium' & !species == 'Festuca.arundinacea')

forage$species <- factor(forage$species, levels = c('Forage.grass.1','Brassica.spp.','Chicorium.endivia','Plantago.lanceolata','Trifolium.spp.'))
levels(forage$species) <- c('Grass','Brassica spp.','Chicorium','Plantago','Trifolium spp.')

# subset to just perennial forage mix zone
forage.forage <- subset(forage, pasture == 'Perennial forage mix')

ggplot(forage.forage, aes(x = month, y = cover, fill = species)) +
  geom_boxplot() +
  # geom_vline(xintercept = )
  labs(x = 'Month', y = '(%) Cover', fill = '') +
  scale_y_log10() +
  theme_bw() +
  theme(legend.margin = margin(0,0,0,0),
        legend.position = 'bottom') +
  # legend.title = element_blank())
  facet_grid(pasture ~ grazing)
dev.off()

# calculate means and SEs
sumSE.forage <- summarySE(forage.forage, measurevar = "cover", groupvars = c('month','pasture','grazing','species'))
# median.forage <- aggregate(cover ~ month + pasture + grazing + species, data = forage.forage, FUN = median)
sumSE.forage2 <- subset(sumSE.forage, !species == 'Grass')


png(paste0(Rout,"/forage_forbs.png"),height = 2.5, width = 6, units = "in", res = 600)
ggplot() +
  # geom_point(data = natives, aes(x = month, y = cover, group = block, shape = grazing, col = species), alpha = 0.2, position = position_dodge(0.8), size = 2) +
  # geom_errorbar(data = sumSE.species, aes(x = month, ymin = cover - se, ymax = cover + se, col = species), width = 0) +
  geom_point(data = sumSE.forage2, aes(x = month, y = cover, col = species)) +
  geom_line(data = sumSE.forage2, aes(x = month, y = cover, col = species, group = species), alpha = 0.5)+
  scale_color_manual(values = c('forestgreen','deeppink','blue3','tan1','lightskyblue','lightpink','chartreuse3','gold','darkred','darkorchid','darkgrey')) +
  labs(x = 'Month', y = 'Mean (%) cover', col = '', shape = '') +
  # scale_y_log10()+
  theme_bw() +
  theme(legend.margin = margin(0,0,0,0)) +
  # legend.title = element_blank())
  facet_grid(pasture ~ grazing)
dev.off()

#### grass, bareground, litter ####
background <- cover[,-c(13:14)]
background <- gather(background, category, cover, Bareground:Forage.grass, factor_key=TRUE)

background$category <- factor(background$category, levels = c('Bareground','Litter','Moss','Forage.grass'))
levels(background$category) <- c('Bare ground','Litter','Moss','Grass')


ggplot(background, aes(x = month, y = cover, fill = category)) +
  geom_boxplot() +
  # geom_vline(xintercept = )
  labs(x = 'Month', y = '(%) Cover', fill = '') +
  scale_y_log10() +
  theme_bw() +
  theme(legend.margin = margin(0,0,0,0),
        legend.position = 'bottom') +
  # legend.title = element_blank())
  facet_grid(pasture ~ grazing)
dev.off()

# calculate means and SEs
sumSE.background <- summarySE(background, measurevar = "cover", groupvars = c('month','pasture','grazing','category'))

png(paste0(Rout,"/background.png"),height = 4.5, width = 6, units = "in", res = 600)
ggplot() +
  # geom_point(data = natives, aes(x = month, y = cover, group = block, shape = grazing, col = species), alpha = 0.2, position = position_dodge(0.8), size = 2) +
  geom_errorbar(data = sumSE.background, aes(x = month, ymin = cover - sd, ymax = cover + sd, col = category), width = 0.2) +
  geom_point(data = sumSE.background, aes(x = month, y = cover, col = category)) +
  geom_line(data = sumSE.background, aes(x = month, y = cover, col = category, group = category), alpha = 0.5)+
  scale_color_manual(values = c('tan1','deeppink','blue3','forestgreen')) +
  labs(x = 'Month', y = 'Mean (%) cover', col = '', shape = '') +
  # scale_y_log10()+
  theme_bw() +
  theme(legend.margin = margin(0,0,0,0)) +
  facet_grid(pasture ~ grazing)
dev.off()
