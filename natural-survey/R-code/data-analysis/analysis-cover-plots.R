# Analyzing Pisgah and SEM sites separately
library(vegan)
library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(MuMIn)
library(Rmisc)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(cowplot)
library(sjPlot)

# set working directory to cleaned data folder
setwd("C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/Natural Site Survey/Reed_USDA-NIFA_NaturalSurvey/Data-cleaned")

# set path to R-output folder within Data-analysis
Rout <- "C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/Natural Site Survey/Reed_USDA-NIFA_NaturalSurvey/R-code/Data-analysis/R-output"

# load in .csv files from 'processing-cover-plots' R script
cover <- read.csv('community_cover_plots-processed.csv')
relcover <- read.csv('community_relcover_plots-processed.csv')

cover$block <- factor(cover$block, levels = c('1','2','3','4','5',
                                              '6','7','8','9','10'))
cover$block2 <- factor(rep(c(rep("a",6),rep("b",6),rep("c",6),rep("d",6),rep("e",6)),2))
relcover$block <- factor(relcover$block, levels = c('1','2','3','4','5',
                                              '6','7','8','9','10'))
relcover$block2 <- factor(rep(c(rep("a",6),rep("b",6),rep("c",6),rep("d",6),rep("e",6)),2))

cover$transect <- factor(cover$transect, levels = c('1a','1b','2a','2b','3a','3b','4a','4b','5a','5b',
                                                    '6a','6b','7a','7b','8a','8b','9a','9b','10a','10b'))
relcover$transect <- factor(relcover$transect, levels = c('1a','1b','2a','2b','3a','3b','4a','4b','5a','5b',
                                                          '6a','6b','7a','7b','8a','8b','9a','9b','10a','10b'))

# separate data by sites
cover.Pisgah <- subset(cover, site == 'Pisgah')
cover.SEM <- subset(cover, site == 'South Eugene Meadows')
relcover.Pisgah <- subset(relcover, site == 'Pisgah')
relcover.SEM <- subset(relcover, site == 'South Eugene Meadows')

# subsample plots were 1m2. there were 3 plots per stand (meaning 6 plots per block; 60 total across sites).
# within plots, I estimated cover to species level, estimated litter and bareground, and measured litter depth.

# set up logit back-transform function for summary table purposes
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

#### Richness ####
# histogram of richness
hist(log(cover.Pisgah$rich))
hist((cover.SEM$rich))

# calculate transect means and CIs
sumSE.rich <- summarySE(cover, measurevar = "rich", groupvars = c('block','stand','site')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-plots/richness-plots.pdf"), height = 2.5, width = 4)
ggplot() +
  geom_boxplot(data = cover, aes(x = stand, y = rich, fill = stand)) +
  geom_point(data = cover, aes(x = stand, y = rich, group = block2, shape = block2), alpha = 0.1, position = position_dodge(0.8), size = 2) +
  geom_line(data = cover, aes(x = stand, y = rich, group = transect), alpha = 0.1, position = position_dodge(0.8))+
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_shape_manual(values = rep(c(16,17,15,3,7),2)) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  labs(x = 'Stand', y = "Richness (1 m2 subplots)") +
  facet_grid(.~ site, scales = 'free')
dev.off()

# model richness
m.rich.Pisgah <- lmer((rich) ~ stand + (1|transect), data = cover.Pisgah)
tab_model(m.rich.Pisgah, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.rich.Pisgah)

m.rich.SEM <- lmer((rich) ~ stand + (1|transect), data = cover.SEM)
tab_model(m.rich.SEM, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.rich.SEM)

# make a table of the results
tab_model(m.rich.Pisgah, m.rich.SEM, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

#### Native richness ####
# histogram of native richness
hist((cover.Pisgah$rich.nat))
hist(log(cover.SEM$rich.nat))

# calculate stand means and CIs
sumSE.rich.nat <- summarySE(cover, measurevar = "rich.nat", groupvars = c('block','stand','site')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-plots/richness-native-plots.pdf"), height = 2.5, width = 4)
ggplot() +
  geom_boxplot(data = cover, aes(x = stand, y = rich.nat, fill = stand)) +
  geom_point(data = cover, aes(x = stand, y = rich.nat, group = block2, shape = block2), alpha = 0.1, position = position_dodge(0.8), size = 2) +
  geom_line(data = cover, aes(x = stand, y = rich.nat, group = transect), alpha = 0.1, position = position_dodge(0.8))+
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_shape_manual(values = rep(c(16,17,15,3,7),2)) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  labs(x = 'Stand', y = "Native richness (1 m2 subplots)") +
  facet_grid(.~ site, scales = 'free')
dev.off()

# model native richness
m.rich.nat.Pisgah <- lmer((rich.nat) ~ stand + (1|block), data = cover.Pisgah)
m.rich.nat.Pisgah <- lmer((rich.nat) ~ stand + (1|transect), data = cover.Pisgah)
tab_model(m.rich.nat.Pisgah, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.rich.nat.Pisgah)

m.rich.nat.SEM <- lmer((rich.nat) ~ stand + (1|block), data = cover.SEM)
m.rich.nat.SEM <- lmer((rich.nat) ~ stand + (1|transect), data = cover.SEM)
tab_model(m.rich.nat.SEM, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.rich.nat.SEM)

# make a table of the results
tab_model(m.rich.nat.Pisgah, m.rich.nat.SEM, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

#### Forb richness ####
# histogram of Forb richness
hist((cover.Pisgah$rich.forb))
hist((cover.SEM$rich.forb))

# calculate stand means and CIs
sumSE.rich.forb <- summarySE(cover, measurevar = "rich.forb", groupvars = c('block','stand','site')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-plots/richness-forb-plots.pdf"), height = 2.5, width = 4)
ggplot() +
  geom_boxplot(data = cover, aes(x = stand, y = rich.forb, fill = stand)) +
  geom_point(data = cover, aes(x = stand, y = rich.forb, group = block2, shape = block2), alpha = 0.1, position = position_dodge(0.8), size = 2) +
  geom_line(data = cover, aes(x = stand, y = rich.forb, group = transect), alpha = 0.1, position = position_dodge(0.8))+
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_shape_manual(values = rep(c(16,17,15,3,7),2)) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  labs(x = 'Stand', y = 'Forb richness (1 m2 subplots)') +
  guides(shape = FALSE) +
  facet_grid(.~ site,scales = 'free')
dev.off()

# model forb richness
m.rich.forb.Pisgah <- lmer((rich.forb) ~ stand + (1|block), data = cover.Pisgah)
m.rich.forb.Pisgah <- lmer((rich.forb) ~ stand + (1|transect), data = cover.Pisgah)
tab_model(m.rich.forb.Pisgah, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.rich.forb.Pisgah)

m.rich.forb.SEM <- lmer((rich.forb) ~ stand + (1|block), data = cover.SEM)
m.rich.forb.SEM <- lmer((rich.forb) ~ stand + (1|transect), data = cover.SEM)
tab_model(m.rich.forb.SEM, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.rich.forb.SEM)

# make a table of the results
tab_model(m.rich.forb.Pisgah, m.rich.forb.SEM, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

#### Native forb richness ####
# histogram of natforb richness
hist((cover.Pisgah$rich.natforb))
hist((cover.SEM$rich.natforb))

# calculate stand means and CIs
sumSE.rich.natforb <- summarySE(cover, measurevar = "rich.natforb", groupvars = c('block','stand','site')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-plots/richness-natforb-plots.pdf"), height = 2.5, width = 4)
ggplot() +
  geom_boxplot(data = cover, aes(x = stand, y = rich.natforb, fill = stand)) +
  geom_point(data = cover, aes(x = stand, y = rich.natforb, group = block2, shape = block2), alpha = 0.1, position = position_dodge(0.8), size = 2) +
  geom_line(data = cover, aes(x = stand, y = rich.natforb, group = transect), alpha = 0.1, position = position_dodge(0.8))+
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_shape_manual(values = rep(c(16,17,15,3,7),2)) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  labs(x = 'Stand', y = 'Native forb richness (1 m2 subplots)') +
  guides(shape = FALSE) +
  facet_grid(.~ site,scales = 'free')
dev.off()

# model natforb richness
m.rich.natforb.Pisgah <- lmer((rich.natforb) ~ stand + (1|block), data = cover.Pisgah)
m.rich.natforb.Pisgah <- lmer((rich.natforb) ~ stand + (1|transect), data = cover.Pisgah)
tab_model(m.rich.natforb.Pisgah, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.rich.natforb.Pisgah)

m.rich.natforb.SEM <- lmer((rich.natforb) ~ stand + (1|block), data = cover.SEM)
m.rich.natforb.SEM <- lmer((rich.natforb) ~ stand + (1|transect), data = cover.SEM)
tab_model(m.rich.natforb.SEM, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.rich.natforb.SEM)

# make a table of the results
tab_model(m.rich.natforb.Pisgah, m.rich.natforb.SEM, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

#### Simpson's Diversity ####
# histogram of diversity
hist(logit(cover.Pisgah$diversity))
hist(logit(cover.SEM$diversity))

# calculate stand means and CIs
sumSE.diversity <- summarySE(cover, measurevar = "diversity", groupvars = c('block','stand','site')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-plots/diversity-plots.pdf"), height = 2.5, width = 4)
ggplot() +
  geom_boxplot(data = cover, aes(x = stand, y = diversity, fill = stand)) +
  geom_point(data = cover, aes(x = stand, y = diversity, group = block2, shape = block2), alpha = 0.1, position = position_dodge(0.8), size = 2) +
  geom_line(data = cover, aes(x = stand, y = diversity, group = transect), alpha = 0.1, position = position_dodge(0.8))+
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_shape_manual(values = rep(c(16,17,15,3,7),2)) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  labs(x = 'Stand', y = "Simpson's diversity (1 m2 subplots)") +
  # guides(shape = FALSE) +
  facet_grid(.~ site,scales = 'free')
dev.off()

# model diversity
m.diversity.Pisgah <- lmer((diversity) ~ stand + (1|block), data = cover.Pisgah)
m.diversity.Pisgah <- lmer((diversity) ~ stand + (1|transect), data = cover.Pisgah)
tab_model(m.diversity.Pisgah, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.diversity.Pisgah)

m.diversity.SEM <- lmer((diversity) ~ stand + (1|block), data = cover.SEM)
m.diversity.SEM <- lmer((diversity) ~ stand + (1|transect), data = cover.SEM)
tab_model(m.diversity.SEM, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.diversity.SEM)

# make a table of the results
tab_model(m.diversity.Pisgah, m.diversity.SEM, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

#### Evenness ####
# histogram of evenness
hist(logit(cover.Pisgah$evenness))
hist(logit(cover.SEM$evenness))

# calculate stand means and CIs
sumSE.evenness <- summarySE(cover, measurevar = "evenness", groupvars = c('block','stand','site')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-plots/evenness-plots.pdf"), height = 2.5, width = 4)
ggplot() +
  geom_boxplot(data = cover, aes(x = stand, y = evenness, fill = stand)) +
  geom_point(data = cover, aes(x = stand, y = evenness, group = block2, shape = block2), alpha = 0.1, position = position_dodge(0.8), size = 2) +
  geom_line(data = cover, aes(x = stand, y = evenness, group = transect), alpha = 0.1, position = position_dodge(0.8))+
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_shape_manual(values = rep(c(16,17,15,3,7),2)) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  labs(x = 'Stand', y = "Simpson's evenness (1 m2 subplots)") +
  guides(shape = FALSE) +
  facet_grid(.~ site,scales = 'free')
dev.off()

# model evenness
m.evenness.Pisgah <- lmer((evenness) ~ stand + (1|block), data = cover.Pisgah)
m.evenness.Pisgah <- lmer((evenness) ~ stand + (1|transect), data = cover.Pisgah)
tab_model(m.evenness.Pisgah, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.evenness.Pisgah)

m.evenness.SEM <- lmer((evenness) ~ stand + (1|block), data = cover.SEM)
m.evenness.SEM <- lmer((evenness) ~ stand + (1|transect), data = cover.SEM)
tab_model(m.evenness.SEM, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.evenness.SEM)

# make a table of the results
tab_model(m.evenness.Pisgah, m.evenness.SEM, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

#### Percent bare ground ####
hist(log(cover.Pisgah$bareground))
hist(log(cover.SEM$bareground))

# calculate stand means and CIs
sumSE.bareground <- summarySE(cover, measurevar = "bareground", groupvars = c('block','stand','site')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-plots/bareground-plots.pdf"), height = 2.5, width = 4)
ggplot() +
  geom_boxplot(data = cover, aes(x = stand, y = bareground, fill = stand)) +
  geom_point(data = cover, aes(x = stand, y = bareground, group = block2, shape = block2), alpha = 0.1, position = position_dodge(0.8), size = 2) +
  geom_line(data = cover, aes(x = stand, y = bareground, group = transect), alpha = 0.1, position = position_dodge(0.8))+
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_shape_manual(values = rep(c(16,17,15,3,7),2)) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  labs(x = 'Stand', y = "% Bareground (1 m2 subplots)") +
  facet_grid(.~ site, scales = 'free')
dev.off()

# model bareground
m.bareground.Pisgah <- lmer((bareground) ~ stand + (1|block), data = cover.Pisgah)
m.bareground.Pisgah <- lmer((bareground) ~ stand + (1|transect), data = cover.Pisgah)
tab_model(m.bareground.Pisgah, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.bareground.Pisgah)

m.bareground.SEM <- lmer((bareground) ~ stand + (1|block), data = cover.SEM)
m.bareground.SEM <- lmer((bareground) ~ stand + (1|transect), data = cover.SEM)
tab_model(m.bareground.SEM, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.bareground.SEM)

# make a table of the results
tab_model(m.bareground.Pisgah, m.bareground.SEM, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

#### Percent litter ####
hist(cover.Pisgah$litter_cover)
hist(cover.SEM$litter_cover)

# calculate stand means and CIs
sumSE.litter <- summarySE(cover, measurevar = "litter_cover", groupvars = c('block','stand','site')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-plots-litter-pct.pdf"), height = 2.5, width = 4)
ggplot() +
  geom_boxplot(data = cover, aes(x = stand, y = litter_cover, fill = stand)) +
  geom_point(data = cover, aes(x = stand, y = litter_cover, group = block2, shape = block2), alpha = 0.1, position = position_dodge(0.8), size = 2) +
  geom_line(data = cover, aes(x = stand, y = litter_cover, group = transect), alpha = 0.1, position = position_dodge(0.8))+
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_shape_manual(values = rep(c(16,17,15,3,7),2)) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  labs(x = 'Stand', y = "% Litter (1 m2 subplots)") +
  facet_grid(.~ site, scales = 'free')
dev.off()

# model percent litter
m.litter.Pisgah <- lmer((litter_cover) ~ stand + (1|block), data = cover.Pisgah)
m.litter.Pisgah <- lmer((litter_cover) ~ stand + (1|transect), data = cover.Pisgah)
tab_model(m.litter.Pisgah, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.litter.Pisgah)

m.litter.SEM <- lmer((litter_cover) ~ stand + (1|block), data = cover.SEM)
m.litter.SEM <- lmer((litter_cover) ~ stand + (1|transect), data = cover.SEM)
tab_model(m.litter.SEM, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.litter.SEM)

# make a table of the results
tab_model(m.litter.Pisgah, m.litter.SEM, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)


#### Litter depth ####
hist(cover.Pisgah$litter_depth)
hist(cover.SEM$litter_depth)

# calculate stand means and CIs
sumSE.litter <- summarySE(cover, measurevar = "litter_depth", groupvars = c('block','stand','site')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-plots/litter-plots.pdf"), height = 2.5, width = 4)
ggplot() +
  geom_boxplot(data = cover, aes(x = stand, y = litter_depth, fill = stand)) +
  geom_point(data = cover, aes(x = stand, y = litter_depth, group = block2, shape = block2), alpha = 0.1, position = position_dodge(0.8), size = 2) +
  geom_line(data = cover, aes(x = stand, y = litter_depth, group = transect), alpha = 0.1, position = position_dodge(0.8))+
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_shape_manual(values = rep(c(16,17,15,3,7),2)) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  labs(x = 'Stand', y = "Litter depth [cm] (1 m2 subplots)") +
  facet_grid(.~ site, scales = 'free')
dev.off()

# model litter depth
m.litter.Pisgah <- lmer((litter_depth) ~ stand + (1|block), data = cover.Pisgah)
m.litter.Pisgah <- lmer((litter_depth) ~ stand + (1|transect), data = cover.Pisgah)
tab_model(m.litter.Pisgah, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.litter.Pisgah)

m.litter.SEM <- lmer((litter_depth) ~ stand + (1|block), data = cover.SEM)
m.litter.SEM <- lmer((litter_depth) ~ stand + (1|transect), data = cover.SEM)
tab_model(m.litter.SEM, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.litter.SEM)

# make a table of the results
tab_model(m.litter.Pisgah, m.litter.SEM, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

#### Percent moss ####
hist(log(cover.Pisgah$moss))
hist(log(cover.SEM$moss))

# calculate stand means and CIs
sumSE.moss <- summarySE(cover, measurevar = "moss", groupvars = c('block','stand','site')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-plots/moss-plots.pdf"), height = 2.5, width = 4)
ggplot(data = cover) +
  geom_boxplot(data = cover, aes(x = stand, y = moss, fill = stand)) +
  geom_point(data = cover, aes(x = stand, y = moss, group = block2, shape = block2), alpha = 0.1, position = position_dodge(0.8), size = 2) +
  geom_line(data = cover, aes(x = stand, y = moss, group = transect), alpha = 0.1, position = position_dodge(0.8))+
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_shape_manual(values = rep(c(16,17,15,3,7),2)) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  labs(x = 'Stand', y = "% Moss (1 m2 subplots)") +
  facet_grid(.~ site, scales = 'free')
dev.off()

# model moss
m.moss.Pisgah <- lmer((moss) ~ stand + (1|block), data = cover.Pisgah)
m.moss.Pisgah <- lmer((moss) ~ stand + (1|transect), data = cover.Pisgah)
tab_model(m.moss.Pisgah, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.moss.Pisgah)

m.moss.SEM <- lmer((moss) ~ stand + (1|block), data = cover.SEM)
m.moss.SEM <- lmer((moss) ~ stand + (1|transect), data = cover.SEM)
tab_model(m.moss.SEM, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.moss.SEM)

# make a table of the results
tab_model(m.moss.Pisgah, m.moss.SEM, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

#### Percent native (relative cover) ####
hist(relcover.Pisgah$relNative)
hist(relcover.SEM$relNative)

# calculate stand means and CIs
sumSE.relcovNative <- summarySE(relcover, measurevar = "relNative", groupvars = c('block','stand','site')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-plots/relcov-Native-plots.pdf"), height = 2.5, width = 4)
ggplot(data = relcover) +
  geom_boxplot(data = relcover, aes(x = stand, y = relNative, fill = stand)) +
  geom_point(data = relcover, aes(x = stand, y = relNative, group = block2, shape = block2), alpha = 0.1, position = position_dodge(0.8), size = 2) +
  geom_line(data = relcover, aes(x = stand, y = relNative, group = transect), alpha = 0.1, position = position_dodge(0.8))+
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_shape_manual(values = rep(c(16,17,15,3,7),2)) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  labs(x = 'Stand', y = "Relative % Native (1 m2 subplots)") +
  facet_grid(.~ site, scales = 'free')
dev.off()

# model relative cover native
m.relcovNative.Pisgah <- lmer((relNative) ~ stand + (1|block), data = relcover.Pisgah)
m.relcovNative.Pisgah <- lmer((relNative) ~ stand + (1|transect), data = relcover.Pisgah)
tab_model(m.relcovNative.Pisgah, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.relcovNative.Pisgah)

m.relcovNative.SEM <- lmer((relNative) ~ stand + (1|block), data = relcover.SEM)
m.relcovNative.SEM <- lmer((relNative) ~ stand + (1|transect), data = relcover.SEM)
tab_model(m.relcovNative.SEM, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.relcovNative.SEM)

# make a table of the results
tab_model(m.relcovNative.Pisgah, m.relcovNative.SEM, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)
