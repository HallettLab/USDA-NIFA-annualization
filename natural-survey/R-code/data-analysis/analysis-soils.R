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
library(ggfortify)

# set working directory to cleaned data folder
setwd("C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/Natural Site Survey/Reed_USDA-NIFA_NaturalSurvey/Data-cleaned")

# set path to R-output folder within Data-analysis
Rout <- "C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/Natural Site Survey/Reed_USDA-NIFA_NaturalSurvey/R-code/Data-analysis/R-output"

########## depth and moisture field measurements ###########
# read in environmental data
env <- read.csv('environmental.csv')

# factor all the categorical variables
env[,1:7] <- lapply(env[,1:7], factor)

# subsample plots were 1m2. there were 3 plots per stand (meaning 6 plots per block; 60 total across sites).
env$block <- factor(env$block, levels = c('1','2','3','4','5',
                                          '6','7','8','9','10'))
env$block2 <- factor(rep(c(rep("a",6),rep("b",6),rep("c",6),rep("d",6),rep("e",6)),2))
env$transect <- factor(env$transect, levels = c('1a','1b','2a','2b','3a','3b','4a','4b','5a','5b',
                                                '6a','6b','7a','7b','8a','8b','9a','9b','10a','10b'))

# separate data by sites
env.Pisgah <- subset(env, site == 'Pisgah')
env.SEM <- subset(env, site == 'South Eugene Meadows')

#### soil moisture ####
hist((env.Pisgah$moist_mean))
hist((env.SEM$moist_mean))

# calculate transect means and CIs
sumSE.moist <- summarySE(env, measurevar = "moist_mean", groupvars = c('block','stand','site')) # aggregating to stand, block, site

pdf(paste0(Rout,"/environmental/moisture.pdf"), height = 2.5, width = 4)
ggplot() +
  geom_boxplot(data = env, aes(x = stand, y = moist_mean, fill = stand)) +
  geom_point(data = env, aes(x = stand, y = moist_mean, group = block2, shape = block2), alpha = 0.1, position = position_dodge(0.8), size = 2) +
  geom_line(data = env, aes(x = stand, y = moist_mean, group = transect), alpha = 0.1, position = position_dodge(0.8))+
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_shape_manual(values = rep(c(16,17,15,3,7),2)) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  labs(x = 'Stand', y = "Soil moisture (vwc)") +
  facet_grid(.~ site, scales = 'free')
dev.off()

# model moisture
m.moist.Pisgah <- lmer((moist_mean) ~ stand + (1|block), data = env.Pisgah)
m.moist.Pisgah <- lmer((moist_mean) ~ stand + (1|transect), data = env.Pisgah)
tab_model(m.moist.Pisgah, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.moist.Pisgah)

m.moist.SEM <- lmer((moist_mean) ~ stand + (1|block), data = env.SEM)
m.moist.SEM <- lmer((moist_mean) ~ stand + (1|transect), data = env.SEM)
tab_model(m.moist.SEM, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.moist.SEM)

# make a table of the results
tab_model(m.moist.Pisgah, m.moist.SEM, linebreak = FALSE, p.val = 'satterthwaite')

#### soil depth ####
hist((env.Pisgah$depth_mean))
hist((env.SEM$depth_mean))

# calculate transect means and CIs
sumSE.moist <- summarySE(env, measurevar = "depth_mean", groupvars = c('block','stand','site')) # aggregating to stand, block, site

pdf(paste0(Rout,"/environmental/depth.pdf"), height = 2.5, width = 4)
ggplot() +
  geom_boxplot(data = env, aes(x = stand, y = depth_mean, fill = stand)) +
  geom_point(data = env, aes(x = stand, y = depth_mean, group = block2, shape = block2), alpha = 0.1, position = position_dodge(0.8), size = 2) +
  geom_line(data = env, aes(x = stand, y = depth_mean, group = transect), alpha = 0.1, position = position_dodge(0.8))+
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_shape_manual(values = rep(c(16,17,15,3,7),2)) +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  labs(x = 'Stand', y = "Soil depth (cm)") +
  facet_grid(.~ site, scales = 'free')
dev.off()

# model depth
m.depth.Pisgah <- lmer((depth_mean) ~ stand + (1|block), data = env.Pisgah)
m.depth.Pisgah <- lmer((depth_mean) ~ stand + (1|transect), data = env.Pisgah)
tab_model(m.depth.Pisgah, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.depth.Pisgah)

m.depth.SEM <- lmer((depth_mean) ~ stand + (1|block), data = env.SEM)
m.depth.SEM <- lmer((depth_mean) ~ stand + (1|transect), data = env.SEM)
tab_model(m.depth.SEM, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.depth.SEM)

# make a table of the results
tab_model(m.depth.Pisgah, m.depth.SEM, linebreak = FALSE, p.val = 'satterthwaite')



#### regress soil moisture on depth ####
pdf(paste0(Rout,"/soil-moist_depth.pdf"), height = 2.5, width = 4)
ggplot() +
  geom_point(data = env, aes(x = depth_mean, y = moist_mean, col = stand, shape = block2), size = 2) +
  geom_smooth(data = env, aes(x = depth_mean, y = moist_mean), method = 'lm', se = FALSE, col = 'black') +
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_shape_manual(values = rep(c(16,17,15,3,7),2)) +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  labs(x = 'Soil depth (cm)', y = 'Soil moisture (vwc)') +
  guides(shape = 'none') +
  facet_grid(.~ site)
dev.off()

# model moisture against depth
m.moist2.Pisgah <- lmer((moist_mean) ~ depth_mean + (1|block), data = env.Pisgah)
m.moist2.Pisgah <- lmer((moist_mean) ~ depth_mean + (1|transect), data = env.Pisgah)
tab_model(m.moist2.Pisgah, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.moist2.Pisgah)

m.moist2.SEM <- lmer((moist_mean) ~ depth_mean + (1|block), data = env.SEM)
m.moist2.SEM <- lmer((moist_mean) ~ depth_mean + (1|transect), data = env.SEM)
tab_model(m.moist2.SEM, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.moist2.SEM)

# make a table of the results
tab_model(m.moist2.Pisgah, m.moist2.SEM, linebreak = FALSE, p.val = 'satterthwaite')

########## soil composition analysis ###########
# read in soils data
soils <- read.csv('soils.csv')

# factor all the categorical variables
soils[,1:8] <- lapply(soils[,1:8], factor)

# experiment was set up as a randomized complete block design with subsampling. 
# 10 blocks across 2 sites (5 blocks per site); each block had a paired annual and perennial stand.
# within each stand, I collected data from 3 subsamples (plots) along a transect.
# within each plot, I composited 3 soil samples of 0-10cm depth, and 3 samples of 10-30cm depth IF i could harvest from that depth (uncommon).

# subsample plots were 1m2. there were 3 plots per stand (meaning 6 plots per block; 60 total across sites).
soils$block <- factor(soils$block, levels = c('1','2','3','4','5',
                                          '6','7','8','9','10'))
soils$block2 <- factor(c(rep("a",8),rep("b",9),rep("c",6),rep("d",6),rep("e",9),
                         rep("a",6),rep("b",8),rep("c",8),rep("d",6),rep("e",6)))
soils$transect <- factor(soils$transect, levels = c('1a','1b','2a','2b','3a','3b','4a','4b','5a','5b',
                                                '6a','6b','7a','7b','8a','8b','9a','9b','10a','10b'))

# separate data by sites
soils.Pisgah <- subset(soils, site == 'Pisgah')
soils.SEM <- subset(soils, site == 'South Eugene Meadows')

# separate samples of different depths
# I also had texture run on a subset of the 0-10cm depth samples (one per stand per block).

# subset soils data
texture <- subset(soils, !is.na(sand_pct))
depth10 <- subset(soils, depth=='0-10cm')
depth30 <- subset(soils, depth=='10-30cm')

texture.Pisgah <- subset(texture, site == 'Pisgah')
texture.SEM <- subset(texture, site == 'South Eugene Meadows')
depth10.Pisgah <- subset(depth10, site == 'Pisgah')
depth10.SEM <- subset(depth10, site == 'South Eugene Meadows')

#### percent sand ####
hist((texture.Pisgah$sand_pct))
hist((texture.SEM$sand_pct))

# calculate transect means and CIs
sumSE.sand <- summarySE(texture, measurevar = "sand_pct", groupvars = c('stand','site'))

pdf(paste0(Rout,"/environmental/sand.pdf"), height = 2.5, width = 4)
ggplot() +
  geom_boxplot(data = texture, aes(x = stand, y = sand_pct, fill = stand)) +
  # geom_point(data = texture, aes(x = stand, y = sand_pct, group = block2, shape = block2), alpha = 0.1, position = position_dodge(0.8), size = 2) +
  # geom_line(data = texture, aes(x = stand, y = sand_pct, group = transect), alpha = 0.1, position = position_dodge(0.8))+
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_shape_manual(values = rep(c(16,17,15,3,7),2)) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank()) +
  labs(x = 'Stand', y = "Percent sand") +
  facet_grid(.~ site, scales = 'free')
dev.off()

# model sand
m.sand.Pisgah <- lm((sand_pct) ~ stand, data = texture.Pisgah)
tab_model(m.sand.Pisgah, linebreak = FALSE, p.val = 'satterthwaite')

m.sand.SEM <- lm((sand_pct) ~ stand, data = texture.SEM)
tab_model(m.sand.SEM, linebreak = FALSE, p.val = 'satterthwaite')

# make a table of the results
tab_model(m.sand.Pisgah, m.sand.SEM, linebreak = FALSE, p.val = 'satterthwaite')


#### percent silt ####
hist((texture.Pisgah$silt_pct))
hist((texture.SEM$silt_pct))

# calculate transect means and CIs
sumSE.silt <- summarySE(texture, measurevar = "silt_pct", groupvars = c('stand','site'))

pdf(paste0(Rout,"/environmental/silt.pdf"), height = 2.5, width = 4)
ggplot() +
  geom_boxplot(data = texture, aes(x = stand, y = silt_pct, fill = stand)) +
  # geom_point(data = texture, aes(x = stand, y = silt_pct, group = block2, shape = block2), alpha = 0.1, position = position_dodge(0.8), size = 2) +
  # geom_line(data = texture, aes(x = stand, y = silt_pct, group = transect), alpha = 0.1, position = position_dodge(0.8))+
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_shape_manual(values = rep(c(16,17,15,3,7),2)) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank()) +
  labs(x = 'Stand', y = "Percent silt") +
  facet_grid(.~ site, scales = 'free')
dev.off()

# model silt
m.silt.Pisgah <- lm((silt_pct) ~ stand, data = texture.Pisgah)
tab_model(m.silt.Pisgah, linebreak = FALSE, p.val = 'satterthwaite')

m.silt.SEM <- lm((silt_pct) ~ stand, data = texture.SEM)
tab_model(m.silt.SEM, linebreak = FALSE, p.val = 'satterthwaite')

# make a table of the results
tab_model(m.silt.Pisgah, m.silt.SEM, linebreak = FALSE, p.val = 'satterthwaite')


#### percent clay ####
hist((texture.Pisgah$clay_pct))
hist((texture.SEM$clay_pct))

# calculate transect means and CIs
sumSE.clay <- summarySE(texture, measurevar = "clay_pct", groupvars = c('stand','site'))

pdf(paste0(Rout,"/environmental/clay.pdf"), height = 2.5, width = 4)
ggplot() +
  geom_boxplot(data = texture, aes(x = stand, y = clay_pct, fill = stand)) +
  # geom_point(data = texture, aes(x = stand, y = clay_pct, group = block2, shape = block2), alpha = 0.1, position = position_dodge(0.8), size = 2) +
  # geom_line(data = texture, aes(x = stand, y = clay_pct, group = transect), alpha = 0.1, position = position_dodge(0.8))+
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_shape_manual(values = rep(c(16,17,15,3,7),2)) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank()) +
  labs(x = 'Stand', y = "Percent clay") +
  facet_grid(.~ site, scales = 'free')
dev.off()

# model clay
m.clay.Pisgah <- lm((clay_pct) ~ stand, data = texture.Pisgah)
tab_model(m.clay.Pisgah, linebreak = FALSE, p.val = 'satterthwaite')

m.clay.SEM <- lm((clay_pct) ~ stand, data = texture.SEM)
tab_model(m.clay.SEM, linebreak = FALSE, p.val = 'satterthwaite')

# make a table of the results
tab_model(m.clay.Pisgah, m.clay.SEM, linebreak = FALSE, p.val = 'satterthwaite')
#### percent carbon, 0-10cm ####
hist((depth10.Pisgah$C_pct))
hist((depth10.SEM$C_pct))

# calculate transect means and CIs
sumSE.C_pct <- summarySE(depth10, measurevar = "C_pct", groupvars = c('block','stand','site'), na.rm = T) # aggregating to stand, block, site

pdf(paste0(Rout,"/environmental/carbon.pdf"), height = 2.5, width = 4)
ggplot() +
  geom_boxplot(data = depth10, aes(x = stand, y = C_pct, fill = stand)) +
  geom_point(data = depth10, aes(x = stand, y = C_pct, group = block2, shape = block2), alpha = 0.1, position = position_dodge(0.8), size = 2) +
  geom_line(data = depth10, aes(x = stand, y = C_pct, group = transect), alpha = 0.1, position = position_dodge(0.8))+
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_shape_manual(values = rep(c(16,17,15,3,7),2)) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  labs(x = 'Stand', y = "Percent carbon") +
  facet_grid(.~ site, scales = 'free')
dev.off()

# model carbon
m.carbon.Pisgah <- lmer((C_pct) ~ stand + (1|block), data = depth10.Pisgah)
m.carbon.Pisgah <- lmer((C_pct) ~ stand + (1|transect), data = depth10.Pisgah)
tab_model(m.carbon.Pisgah, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.carbon.Pisgah)

m.carbon.SEM <- lmer((C_pct) ~ stand + (1|block), data = depth10.SEM)
m.carbon.SEM <- lmer((C_pct) ~ stand + (1|transect), data = depth10.SEM)
tab_model(m.carbon.SEM, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.carbon.SEM)

# make a table of the results
tab_model(m.carbon.Pisgah, m.carbon.SEM, linebreak = FALSE, p.val = 'satterthwaite')

#### percent nitrogen, 0-10cm ####
hist((depth10.Pisgah$N_pct))
hist((depth10.SEM$N_pct))

# calculate transect means and CIs
sumSE.N_pct <- summarySE(depth10, measurevar = "N_pct", groupvars = c('block','stand','site'), na.rm = T) # aggregating to stand, block, site

pdf(paste0(Rout,"/environmental/nitrogen.pdf"), height = 2.5, width = 4)
ggplot() +
  geom_boxplot(data = depth10, aes(x = stand, y = N_pct, fill = stand)) +
  geom_point(data = depth10, aes(x = stand, y = N_pct, group = block2, shape = block2), alpha = 0.1, position = position_dodge(0.8), size = 2) +
  geom_line(data = depth10, aes(x = stand, y = N_pct, group = transect), alpha = 0.1, position = position_dodge(0.8))+
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_shape_manual(values = rep(c(16,17,15,3,7),2)) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  labs(x = 'Stand', y = "Percent nitrogen") +
  facet_grid(.~ site, scales = 'free')
dev.off()

# model nitrogen
m.nitrogen.Pisgah <- lmer((N_pct) ~ stand + (1|block), data = depth10.Pisgah)
m.nitrogen.Pisgah <- lmer((N_pct) ~ stand + (1|transect), data = depth10.Pisgah)
tab_model(m.nitrogen.Pisgah, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.nitrogen.Pisgah)

m.nitrogen.SEM <- lmer((N_pct) ~ stand + (1|block), data = depth10.SEM)
m.nitrogen.SEM <- lmer((N_pct) ~ stand + (1|transect), data = depth10.SEM)
tab_model(m.nitrogen.SEM, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.nitrogen.SEM)

# make a table of the results
tab_model(m.nitrogen.Pisgah, m.nitrogen.SEM, linebreak = FALSE, p.val = 'satterthwaite')


#### percent OM, 0-10cm ####
hist((depth10.Pisgah$OM_pct))
hist((depth10.SEM$OM_pct))

# calculate transect means and CIs
sumSE.OM_pct <- summarySE(depth10, measurevar = "OM_pct", groupvars = c('block','stand','site'), na.rm = T) # aggregating to stand, block, site

pdf(paste0(Rout,"/environmental/OM.pdf"), height = 2.5, width = 4)
ggplot() +
  geom_boxplot(data = depth10, aes(x = stand, y = OM_pct, fill = stand)) +
  geom_point(data = depth10, aes(x = stand, y = OM_pct, group = block2, shape = block2), alpha = 0.1, position = position_dodge(0.8), size = 2) +
  geom_line(data = depth10, aes(x = stand, y = OM_pct, group = transect), alpha = 0.1, position = position_dodge(0.8))+
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  scale_shape_manual(values = rep(c(16,17,15,3,7),2)) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  labs(x = 'Stand', y = "Percent OM") +
  facet_grid(.~ site, scales = 'free')
dev.off()

# model OM
m.OM.Pisgah <- lmer((OM_pct) ~ stand + (1|block), data = depth10.Pisgah)
m.OM.Pisgah <- lmer((OM_pct) ~ stand + (1|transect), data = depth10.Pisgah)
tab_model(m.OM.Pisgah, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.OM.Pisgah)

m.OM.SEM <- lmer((OM_pct) ~ stand + (1|block), data = depth10.SEM)
m.OM.SEM <- lmer((OM_pct) ~ stand + (1|transect), data = depth10.SEM)
tab_model(m.OM.SEM, linebreak = FALSE, p.val = 'satterthwaite')
AICc(m.OM.SEM)

# make a table of the results
tab_model(m.OM.Pisgah, m.OM.SEM, linebreak = FALSE, p.val = 'satterthwaite')

########## PCA of all soil components ##########
# aggregate data to transect level
env.agg <- aggregate(cbind(moist_mean, depth_mean, moist_sd, depth_sd) ~ site + block + stand + transect, data = env, mean)
depth10.agg1 <- aggregate(cbind(C_pct, N_pct, OM_pct) ~ site + block + stand + transect, data = depth10[!is.na(depth10$C_pct),], mean)
depth10.agg2 <- aggregate(cbind(sand_pct, silt_pct, clay_pct) ~ site + block + stand + transect, data = depth10, mean)

soils.agg <- cbind(env.agg, depth10.agg1[,c('C_pct','N_pct','OM_pct')], depth10.agg2[,c("sand_pct",'silt_pct','clay_pct')])
soils.agg.sc <- cbind(soils.agg[,1:4], scale(soils.agg[,5:14]))

# run a PCA on scaled data
##### ordination figure #####
pca <- prcomp(~moist_mean + depth_mean + C_pct + N_pct + sand_pct + silt_pct + clay_pct, data = soils.agg.sc)
summary(pca)
rownames(pca$rotation)<-c('% Moisture','Depth (cm)','% C','% N','% Sand','% Silt', '% Clay')
soils.agg.sc$site.stand <- paste(soils.agg$site, soils.agg$stand, sep = ' - ')

pcafig <- autoplot(pca,data=soils.agg.sc,group = 'site.stand', colour='stand',shape='site',
                   frame = TRUE,frame.colour='stand',frame.fill='stand',
                   loadings=TRUE,loadings.label=TRUE,
                   loadings.colour = 'black',loadings.label.colour='black', 
                   loadings.label.size=4,scale=0,
                   size = 3) +
  scale_colour_manual(values=c("goldenrod1","darkolivegreen4")) +
  scale_fill_manual(values=c("goldenrod1","darkolivegreen4")) +
  # xlab(expression(atop('Cooler          ' %->% '          Warmer', 'PC1 (50.02%)')))+
  # ylab(expression(atop('PC2 (23.77%)','Drier          ' %->% '          Moister')))+
  theme_bw()+theme(legend.title=element_blank(),
                   legend.background=element_blank(),
                   legend.position="bottom")
pcafig
pdf(paste0(Rout,"/PCA-soils.pdf"), height=4, width=5)
png(filename = paste0(Rout,"/PCA-soils.png"), res=600, height=4, width=5, units="in")
pcafig
dev.off()
?autoplot
