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

# set working directory to cleaned data folder
setwd("C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/MBF Experiment/Reed_USDA-NIFA_MBF Experiment/Data-cleaned")

# set path to R-output folder within Data-analysis
Rout <- "C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/MBF Experiment/Reed_USDA-NIFA_MBF Experiment/R-code/Data-analysis/R-output"

# read in veg height data
height <- read.csv('VegHeight_MBF_compiled.csv')

# factor all the categorical variables
height[,1:9] <- lapply(height[,1:9], factor)

ggplot(height, aes(x = month, y = ht_mean, fill = grazing)) +
  geom_boxplot() +
  labs(x = 'Month', y = 'Vegetation height (cm)', fill = 'Pasture') +
  theme_bw() +
  theme(legend.margin = margin(0,0,0,0)) +
        # legend.title = element_blank())
  facet_grid(. ~ pasture)

##### model vegetation height #####
m.veght <- lmer(ht_mean ~ grazing * pasture * seeded * month + (1|subplot) + (1|block), data = height, na.action = na.fail, REML = FALSE)
summary(m.veght) 
Anova(m.veght, type = 3)

drg.veght <- dredge(m.veght, fixed = c('grazing', 'pasture'), rank = 'AICc')
m.veght1 <- get.models(drg.veght, subset=delta==0)[[1]] # grab the best-supported model
summary(m.veght1)
Anova(m.veght1, type = 3)

emm.veght <- emmeans(m.veght1, pairwise ~ grazing|month*pasture)

##### cover data #####
cov <- read.csv('Cover_MBF_compiled.csv')
cov.long <- read.csv('Cover_MBF_compiled_long.csv')
