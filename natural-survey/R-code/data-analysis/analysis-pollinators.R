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

# read in pollinator data
pol <- read.csv('pollinator.csv')

# factor all the categorical variables
pol[,1:6] <- lapply(pol[,1:6], factor)

# subsample plots were 1m2. there were 3 plots per stand (meaning 6 plots per block; 60 total across sites).
pol$block <- factor(pol$block, levels = c('1','2','3','4','5',
                                          '6','7','8','9','10'))
pol$transect <- factor(pol$transect, levels = c('1a','1b','2a','2b','3a','3b','4a','4b','5a','5b',
                                                '6a','6b','7a','7b','8a','8b','9a','9b','10a','10b'))


# there were so few pollinators, let's aggregate to transect (get rid of plot-level data)
pol1 <- aggregate(visits ~ site + block + stand + transect + plot, data = pol, sum)
pol2 <- aggregate(visits ~ site + block + stand + transect, data = pol, sum)

# separate data by sites
pol2.Pisgah <- subset(pol2, site == 'Pisgah')
pol2.SEM <- subset(pol2, site == 'South Eugene Meadows')

# histogram of pollinator visits
hist((pol2.Pisgah$visits))
hist((pol2.SEM$visits))

# calculate stand means and CIs
sumSE.visits <- summarySE(pol2, measurevar = "visits", groupvars = c('site','stand')) # aggregating to stand, block, site
sumSE.visits2 <- summarySE(pol, measurevar = "visits", groupvars = c('site','stand','transect')) # aggregating to stand, block, site

pdf(paste0(Rout,"/pollinators/visits.pdf"), height = 2.5, width = 4)
ggplot() +
  geom_point(data = pol2, aes(x = stand, y = visits, group = block, shape = block), alpha = 0.1, position = position_dodge(0.6), size = 2) +
  geom_errorbar(data = sumSE.visits, aes(x = stand, ymin = visits - se, ymax = visits + se, col = stand), width = 0) +
  geom_point(data = sumSE.visits, aes(x = stand, y = visits, col = stand), shape = 16, size = 4) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_shape_manual(values = rep(c(16,17,15,3,7),2)) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank()) +
  labs(x = 'Stand', y = 'Pollinator visits') +
  guides(shape = 'none') +
  facet_grid(.~ site,scales = 'free')
dev.off()

# model visits
m.visits.Pisgah <- lm(visits ~ stand, data = pol2.Pisgah)
tab_model(m.visits.Pisgah, linebreak = FALSE)

m.visits.SEM <- lm(visits ~ stand, data = pol2.SEM)
tab_model(m.visits.SEM, linebreak = FALSE)

# make a table of the results
tab_model(m.visits.Pisgah, m.visits.SEM, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)
