# created 2021-09-28
# preliminary analyses for environmental data
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

# read in pollinator data
pol <- read.csv('pollinator.csv')

# factor all the categorical variables
pol[,1:6] <- lapply(pol[,1:6], factor)

# there were so few pollinators, let's aggregate to transect (get rid of plot-level data)
pol2 <- aggregate(visits ~ site + block + stand + transect, data = pol, sum)

# histogram of pollinator visits
hist(pol2$visits)

# calculate stand means and CIs
sumSE.visits <- summarySE(pol2, measurevar = "visits", groupvars = c('site','stand')) # aggregating to stand, block, site

pdf(paste0(Rout,"/pollinators/visits_sampledata.pdf"), height = 2.5, width = 4)
ggplot(data = pol2) +
  geom_line(aes(x = stand, y = visits, group = block), alpha = 0.2)+
  geom_point(aes(x = stand, y = visits, shape = site), alpha = 0.2, size = 1.5) +
  geom_errorbar(data = sumSE.visits, aes(x = stand, ymin = visits - se, ymax = visits + se, col = stand), width = 0) +
  geom_point(data = sumSE.visits, aes(x = stand, y = visits, col = stand), size = 3) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank()) +
  labs(x = 'Stand', y = 'Pollinator visits') +
  guides(shape = FALSE) +
  facet_grid(.~ site,scales = 'free')
dev.off()

# model visits
m.visits1 <- glmer(visits ~ stand + (1|site/block), family = poisson, data = pol2)
summary(m.visits1)
AICc(m.visits1)
plot(m.visits1)
# singular fit issue coming from both site and block... remove from model

m.visits1 <- glm(visits ~ stand, family = poisson, data = pol2)
summary(m.visits1)

# calculate modeled means and confidence intervals
emm.visits1 <- data.frame(emmeans(m.visits1, ~ stand, type = 'response'))
emm.visits1

# make a table of the results
tab_model(m.visits1, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

# plot the modeled estimates! 
pdf(paste0(Rout,"/pollinators/visits_modeled.pdf"), height = 2.5, width = 2.5)
ggplot(data = pol2) +
  geom_line(aes(x = stand, y = visits, group = block), alpha = 0.2)+
  geom_point(aes(x = stand, y = visits, shape = site), alpha = 0.2, size = 1.5) +
  geom_errorbar(data = emm.visits1, aes(x = stand, ymin = asymp.LCL, ymax = asymp.UCL, col = stand), width = 0) +
  geom_point(data = emm.visits1, aes(x = stand, y = rate, col = stand), size = 3.5) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = '', y = 'Pollinator visits') +
  guides(col = FALSE)
dev.off()
