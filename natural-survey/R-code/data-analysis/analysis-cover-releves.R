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

# load in .csv files from 'processing-cover-releves' R script
cover <- read.csv('community_cover_releves-processed.csv')
relcover <- read.csv('community_relcover_releves-processed.csv')

# separate data by sites
cover.Pisgah <- subset(cover, site == 'Pisgah')
cover.SEM <- subset(cover, site == 'South Eugene Meadows')
relcover.Pisgah <- subset(relcover, site == 'Pisgah')
relcover.SEM <- subset(relcover, site == 'South Eugene Meadows')

# large releve plots were 100m2. there was 1 releve per stand (meaning 2 releves per block; 20 total across sites).
# within releves, I estimated cover to species level, estimated litter and bareground. did not measure litter depth.
# releves are mostly going to be useful for richness; maybe diversity.

#### Richness ####
# histogram of richness
hist(cover.Pisgah$rich)
hist(cover.SEM$rich)

# calculate stand means and CIs
sumSE.rich <- summarySE(cover, measurevar = "rich", groupvars = c('site','stand')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-releves/richness-releves.pdf"), height = 2.5, width = 4)
ggplot(data = cover) +
  # geom_line(aes(x = stand, y = rich, group = block), alpha = 0.2)+
  # geom_point(aes(x = stand, y = rich, shape = site, group = block), alpha = 0.2, size = 1.5) +
  geom_boxplot(aes(x = stand, y = rich, fill = stand)) +
  # geom_errorbar(data = sumSE.rich, aes(x = stand, ymin = rich - ci, ymax = rich + ci, col = stand), width = 0) +
  # geom_point(data = sumSE.rich, aes(x = stand, y = rich, col = stand), size = 3) +
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank()) +
  labs(x = 'Stand', y = "Richness (100 m2 relevés)") +
  facet_grid(.~ site, scales = 'free')
dev.off()

# model richness
m.rich.Pisgah <- lm(rich ~ stand, data = cover.Pisgah)
summary(m.rich.Pisgah)

m.rich.SEM <- lm(rich ~ stand, data = cover.SEM)
summary(m.rich.SEM)

# make a table of the results
tab_model(m.rich.Pisgah, m.rich.SEM, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

#### Native richness ####
# histogram of native richness
hist((cover.Pisgah$rich.nat))
hist((cover.SEM$rich.nat))

# calculate stand means and CIs
sumSE.rich.nat <- summarySE(cover, measurevar = "rich.nat", groupvars = c('site','stand')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-releves/richness-native-releves.pdf"), height = 2.5, width = 4)
ggplot(data = cover) +
  # geom_line(aes(x = stand, y = rich.nat, group = block), alpha = 0.2)+
  # geom_point(aes(x = stand, y = rich.nat, shape = site), alpha = 0.2, size = 1.5) +
  geom_boxplot(aes(x = stand, y = rich.nat, fill = stand)) +
  # geom_errorbar(data = sumSE.rich.nat, aes(x = stand, ymin = rich.nat - ci, ymax = rich.nat + ci, col = stand), width = 0) +
  # geom_point(data = sumSE.rich.nat, aes(x = stand, y = rich.nat, col = stand), size = 3) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank()) +
  labs(x = 'Stand', y = 'Native Richness (100 m2 relevés)') +
  guides(shape = FALSE) +
  facet_grid(.~ site,scales = 'free')
dev.off()

# model native richness
m.rich.nat.Pisgah <- lm(rich.nat ~ stand, data = cover.Pisgah)
summary(m.rich.nat.Pisgah)

m.rich.nat.SEM <- lm(rich.nat ~ stand, data = cover.SEM)
summary(m.rich.nat.SEM)

# make a table of the results
tab_model(m.rich.nat.Pisgah, m.rich.nat.SEM, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

#### Forb richness ####
# histogram of Forb richness
hist((cover.Pisgah$rich.forb))
hist((cover.SEM$rich.forb))

# calculate stand means and CIs
sumSE.rich.forb <- summarySE(cover, measurevar = "rich.forb", groupvars = c('site','stand')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-releves/richness-forb-releves.pdf"), height = 2.5, width = 4)
ggplot(data = cover) +
  # geom_line(aes(x = stand, y = rich.forb, group = block), alpha = 0.2)+
  # geom_point(aes(x = stand, y = rich.forb, shape = site), alpha = 0.2, size = 1.5) +
  geom_boxplot(aes(x = stand, y = rich.forb, fill = stand)) +
  # geom_errorbar(data = sumSE.rich.forb, aes(x = stand, ymin = rich.forb - ci, ymax = rich.forb + ci, col = stand), width = 0) +
  # geom_point(data = sumSE.rich.forb, aes(x = stand, y = rich.forb, col = stand), size = 3) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank()) +
  labs(x = 'Stand', y = 'Forb Richness (100 m2 relevés)') +
  guides(shape = FALSE) +
  facet_grid(.~ site,scales = 'free')
dev.off()

# model Forb richness
m.rich.forb.Pisgah <- lm((rich.forb) ~ stand, data = cover.Pisgah)
summary(m.rich.forb.Pisgah)

m.rich.forb.SEM <- lm((rich.forb) ~ stand, data = cover.SEM)
summary(m.rich.forb.SEM)

# make a table of the results
tab_model(m.rich.forb.Pisgah, m.rich.forb.SEM, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

#### Native forb richness ####
# histogram of natforb richness
hist((cover.Pisgah$rich.natforb))
hist((cover.SEM$rich.natforb))

# calculate stand means and CIs
sumSE.rich.natforb <- summarySE(cover, measurevar = "rich.natforb", groupvars = c('site','stand')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-releves/richness-natforb-releves.pdf"), height = 2.5, width = 4)
ggplot(data = cover) +
  # geom_line(aes(x = stand, y = rich.natforb, group = block), alpha = 0.2)+
  # geom_point(aes(x = stand, y = rich.natforb, shape = site), alpha = 0.2, size = 1.5) +
  geom_boxplot(aes(x = stand, y = rich.natforb, fill = stand)) +
  # geom_errorbar(data = sumSE.rich.natforb, aes(x = stand, ymin = rich.natforb - ci, ymax = rich.natforb + ci, col = stand), width = 0) +
  # geom_point(data = sumSE.rich.natforb, aes(x = stand, y = rich.natforb, col = stand), size = 3) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank()) +
  labs(x = 'Stand', y = 'Native Forb Richness (100 m2 relevés)') +
  guides(shape = FALSE) +
  facet_grid(.~ site,scales = 'free')
dev.off()

# model natforb richness
m.rich.natforb.Pisgah <- lm((rich.natforb) ~ stand, data = cover.Pisgah)
summary(m.rich.natforb.Pisgah)

m.rich.natforb.SEM <- lm((rich.natforb) ~ stand, data = cover.SEM)
summary(m.rich.natforb.SEM)

# make a table of the results
tab_model(m.rich.natforb.Pisgah, m.rich.natforb.SEM, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

#### Simpson's Diversity ####
# (probably not reliable for releves since it depends on cover estimates) 

# histogram of diversity
hist(logit(cover.Pisgah$diversity))
hist(logit(cover.SEM$diversity))

# calculate stand means and CIs
sumSE.diversity <- summarySE(cover, measurevar = "diversity", groupvars = c('site','stand')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-releves/diversity-releves.pdf"), height = 2.5, width = 4)
ggplot(data = cover) +
  # geom_line(aes(x = stand, y = diversity, group = block), alpha = 0.2)+
  # geom_point(aes(x = stand, y = diversity, shape = site), alpha = 0.2, size = 1.5) +
  geom_boxplot(aes(x = stand, y = diversity, fill = stand)) +
  # geom_errorbar(data = sumSE.diversity, aes(x = stand, ymin = diversity - ci, ymax = diversity + ci, col = stand), width = 0) +
  # geom_point(data = sumSE.diversity, aes(x = stand, y = diversity, col = stand), size = 3) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank()) +
  labs(x = 'Stand', y = "Simpson's Diversity (100 m2 relevés)") +
  # guides(shape = FALSE) +
  facet_grid(.~ site,scales = 'free')
dev.off()

# model diversity
m.diversity.Pisgah <- lm((diversity) ~ stand, data = cover.Pisgah)
summary(m.diversity.Pisgah)

m.diversity.SEM <- lm((diversity) ~ stand, data = cover.SEM)
summary(m.diversity.SEM)

# make a table of the results
tab_model(m.diversity.Pisgah, m.diversity.SEM, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)


#### Evenness ####
# histogram of evenness
hist(logit(cover.Pisgah$evenness))
hist(logit(cover.SEM$evenness))

# calculate stand means and CIs
sumSE.evenness <- summarySE(cover, measurevar = "evenness", groupvars = c('site','stand')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-releves/evenness-releves.pdf"), height = 2.5, width = 4)
ggplot(data = cover) +
  # geom_line(aes(x = stand, y = evenness, group = block), alpha = 0.2)+
  # geom_point(aes(x = stand, y = evenness, shape = site), alpha = 0.2, size = 1.5) +
  geom_boxplot(aes(x = stand, y = evenness, fill = stand)) +
  # geom_errorbar(data = sumSE.evenness, aes(x = stand, ymin = evenness - ci, ymax = evenness + ci, col = stand), width = 0) +
  # geom_point(data = sumSE.evenness, aes(x = stand, y = evenness, col = stand), size = 3) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank()) +
  labs(x = 'Stand', y = "Simpson's evenness (100 m2 relevés)") +
  guides(shape = FALSE) +
  facet_grid(.~ site,scales = 'free')
dev.off()

# model evenness
m.evenness.Pisgah <- lm((evenness) ~ stand, data = cover.Pisgah)
summary(m.evenness.Pisgah)

m.evenness.SEM <- lm((evenness) ~ stand, data = cover.SEM)
summary(m.evenness.SEM)

# make a table of the results
tab_model(m.evenness.Pisgah, m.evenness.SEM, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

#### Percent bare ground ####
hist(cover.Pisgah$bareground)
hist(cover.SEM$bareground)

# calculate stand means and CIs
sumSE.bareground <- summarySE(cover, measurevar = "bareground", groupvars = c('site','stand')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-releves/bareground-releves.pdf"), height = 2.5, width = 4)
ggplot(data = cover) +
  # geom_line(aes(x = stand, y = bareground, group = block), alpha = 0.2)+
  # geom_point(aes(x = stand, y = bareground, shape = site, group = block), alpha = 0.2, size = 1.5) +
  geom_boxplot(aes(x = stand, y = bareground, fill = stand)) +
  # geom_errorbar(data = sumSE.bareground, aes(x = stand, ymin = bareground - ci, ymax = bareground + ci, col = stand), width = 0) +
  # geom_point(data = sumSE.bareground, aes(x = stand, y = bareground, col = stand), size = 3) +
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank()) +
  labs(x = 'Stand', y = "% Bareground (100 m2 relevés)") +
  facet_grid(.~ site, scales = 'free')
dev.off()

# model bareground
m.bareground.Pisgah <- lm((bareground) ~ stand, data = cover.Pisgah)
summary(m.bareground.Pisgah)

m.bareground.SEM <- lm(bareground ~ stand, data = cover.SEM)
summary(m.bareground.SEM)

# make a table of the results
tab_model(m.bareground.Pisgah, m.bareground.SEM, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

#### Percent litter ####
hist(cover.Pisgah$litter)
hist(cover.SEM$litter)

# calculate stand means and CIs
sumSE.litter <- summarySE(cover, measurevar = "litter_cover", groupvars = c('site','stand')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-releves/litter-releves.pdf"), height = 2.5, width = 4)
ggplot(data = cover) +
  # geom_line(aes(x = stand, y = litter_cover, group = block), alpha = 0.2)+
  # geom_point(aes(x = stand, y = litter_cover, shape = site, group = block), alpha = 0.2, size = 1.5) +
  geom_boxplot(aes(x = stand, y = litter_cover, fill = stand)) +
  # geom_errorbar(data = sumSE.litter, aes(x = stand, ymin = litter_cover - ci, ymax = litter_cover + ci, col = stand), width = 0) +
  # geom_point(data = sumSE.litter, aes(x = stand, y = litter_cover, col = stand), size = 3) +
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank()) +
  labs(x = 'Stand', y = "% Litter (100 m2 relevés)") +
  facet_grid(.~ site, scales = 'free')
dev.off()

# model litter
m.litter.Pisgah <- lm((litter_cover) ~ stand, data = cover.Pisgah)
summary(m.litter.Pisgah)

m.litter.SEM <- lm(litter_cover ~ stand, data = cover.SEM)
summary(m.litter.SEM)

# make a table of the results
tab_model(m.litter.Pisgah, m.litter.SEM, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

#### Percent moss ####
hist(cover.Pisgah$moss)
hist(cover.SEM$moss)

# calculate stand means and CIs
sumSE.moss <- summarySE(cover, measurevar = "moss", groupvars = c('site','stand')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-releves/moss-releves.pdf"), height = 2.5, width = 4)
ggplot(data = cover) +
  # geom_line(aes(x = stand, y = moss, group = block), alpha = 0.2)+
  # geom_point(aes(x = stand, y = moss, shape = site, group = block), alpha = 0.2, size = 1.5) +
  geom_boxplot(aes(x = stand, y = moss, fill = stand)) +
  # geom_errorbar(data = sumSE.moss, aes(x = stand, ymin = moss - ci, ymax = moss + ci, col = stand), width = 0) +
  # geom_point(data = sumSE.moss, aes(x = stand, y = moss, col = stand), size = 3) +
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank()) +
  labs(x = 'Stand', y = "% Moss (100 m2 relevés)") +
  facet_grid(.~ site, scales = 'free')
dev.off()

# model moss
m.moss.Pisgah <- lm((moss) ~ stand, data = cover.Pisgah)
summary(m.moss.Pisgah)

m.moss.SEM <- lm(moss ~ stand, data = cover.SEM)
summary(m.moss.SEM)

# make a table of the results
tab_model(m.moss.Pisgah, m.moss.SEM, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

#### Percent native (relative cover) ####
hist(relcover.Pisgah$NATNative)
hist(relcover.SEM$NATNative)

# calculate stand means and CIs
sumSE.relcovNative <- summarySE(relcover, measurevar = "NATNative", groupvars = c('site','stand')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-releves/relcov-Native-releves.pdf"), height = 2.5, width = 4)
ggplot(data = relcover) +
  geom_boxplot(aes(x = stand, y = NATNative * 100, fill = stand)) +
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank()) +
  labs(x = 'Stand', y = "Relative % Native (100 m2 relevés)") +
  facet_grid(.~ site, scales = 'free')
dev.off()

# model relative cover native
m.relcovNative.Pisgah <- lm((NATNative) ~ stand, data = relcover.Pisgah)
summary(m.relcovNative.Pisgah)

m.relcovNative.SEM <- lm((NATNative) ~ stand, data = relcover.SEM)
summary(m.relcovNative.SEM)

# make a table of the results
tab_model(m.relcovNative.Pisgah, m.relcovNative.SEM, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)