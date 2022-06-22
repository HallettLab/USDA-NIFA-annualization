library(raster)
library(rgdal)
library(ggplot2)
library(tidyr)
library(dplyr)
library(PerformanceAnalytics)
library(sjPlot)
library(spdep)
library(car)
library(quantreg)
library(MuMIn)
library(jtools)
library(wesanderson)

# Is there evidence of a spatiotemporal shift in annual cover than can be related to certain landscape factors?
# What topographic and edaphic conditions are most sensitive to rising annual cover? 

# We hypothesize that Willamette Valley grasslands have shifted toward greater annual (%) cover since 1986,
# concurrent with ongoing climate change. Furthermore, we expect that annual species cover would have been historically
# associated with south-facing, steep, and/or shallower, sandier soils where heat load is high and available water storage is low.
# We hypothesize that annual cover fidelity to these conditions has decreased through time,
# as climate change and a buildup of propagules has allowed for increasing colonization of alternative microsites. 

# set working directory
setwd("C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/Other projects/RAP analysis/R_working_directory")

load('rdata_LandscapeDrivers.RData')

# set directory for figure exports
Rout <- 'C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/Other projects/RAP analysis/R_working_directory/R figure exports/allsites'

# load in anchors shapefile
anchors <- readOGR(dsn = "C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/Other projects/RAP analysis/ArcMap raster exports/WVOPC_Anchors_WGS_1984.shp", layer="WVOPC_Anchors_WGS_1984")

# load in df_allsites
df_allsites <- readRDS(file = "df_polygons.rds")

######## calculate new heatload variables using R function. ########
# df_allsites$heatload is the variable calculated using ArcMap package. In theory, this sets NE facing slopes as coolest and SW facing slopes as warmest. Not convinced it's calculated well. Distribution is hyper normal; doesn't work well for analyses.
# Calculates heat load or potential annual direct incident radiation, using the formulas published in 
# McCune & Keon (2002) based on aspect, slope and latitude.

# the below function calculates heatload index exactly as done in McCune and Keon 2002. I have two separate lines for variable "A" below:
# the first one is used to set NE as coolest, SW as warming; the second is used just for N as coolest and S as warmest.
# I will use these two versions of the function to calculate heatload2 and heatload3, repectively.
hli <- function (aspect, slope, latitude, method = 'heatload', units = 'degrees', equation = 1)
{
  if (units == 'degrees')   # convert degrees to radians
  {
    aspect <- aspect/180*pi
    slope <- slope/180*pi
    aspect[slope == 0] <- 0
    latitude <- latitude/180*pi
  }  
  # A <- if (method == 'heatload') abs (pi - abs (aspect - (5*pi/4))) else pi - abs (aspect-pi)
  A <- pi - abs (aspect-pi)
  S <- slope
  L <- if (length (latitude) == 1) rep (latitude, length (A)) else latitude
  if (equation == 1) res <- exp (-1.467 +1.582*cos(L)*cos(S) -1.500*cos(A)*sin(S)*sin(L) -0.262*sin(L)*sin(S) +0.607*sin(A)*sin(S))
  if (equation == 2) res <- exp (-1.236 +1.350*cos(L)*cos(S) -1.376*cos(A)*sin(S)*sin(L) -0.331*sin(L)*sin(S) +0.375*sin(A)*sin(S))
  if (equation == 3) res <-      +0.339 +0.808*cos(L)*cos(S)                             -0.196*sin(L)*sin(S)                       - 0.482*cos(A)*sin(S)
  return (res)
}
df_allsites$heatload2 <- hli (aspect = df_allsites$aspect, slope = df_allsites$slope, latitude = df_allsites$y, equation = 3)
df_allsites$heatload3 <- hli (aspect = df_allsites$aspect, slope = df_allsites$slope, latitude = df_allsites$y, equation = 3)
hist(df_allsites$heatload,xlim = c(0,1.1), ylim = c(0,1200000))
hist(df_allsites$heatload2,xlim = c(0,1.1), ylim = c(0,1200000))
hist(df_allsites$heatload3,xlim = c(0,1.1), ylim = c(0,1200000))
hist(log(df_allsites$heatload),xlim = c(-0.8,0.1), ylim = c(0,1200000))
hist(log(df_allsites$heatload2),xlim = c(-0.8,0.1), ylim = c(0,1200000))
hist(log(df_allsites$heatload3),xlim = c(-0.8,0.1), ylim = c(0,1200000))

#### assign quantiles to predictor variables ####
# library(gtools)
quantile(df_allsites$heatload, na.rm=T, probs = seq(0,1,(1/3)))
df_allsites$heatload_quants <- gtools::quantcut(df_allsites$heatload, q = 3)
quantile(df_allsites$heatload2, na.rm=T, probs = seq(0,1,(1/3)))
df_allsites$heatload2_quants <- gtools::quantcut(df_allsites$heatload2, q = 3)
quantile(df_allsites$heatload3, na.rm=T, probs = seq(0,1,(1/3)))
df_allsites$heatload3_quants <- gtools::quantcut(df_allsites$heatload3, q = 3)
df_allsites$aspect_cats <- NA
df_allsites[!is.na(df_allsites$aspect),]$aspect_cats <- 'North'
df_allsites[!is.na(df_allsites$aspect) & df_allsites$aspect > 45,]$aspect_cats <- 'East'
df_allsites[!is.na(df_allsites$aspect) & df_allsites$aspect > 135,]$aspect_cats <- 'South'
df_allsites[!is.na(df_allsites$aspect) & df_allsites$aspect > 225,]$aspect_cats <- 'West'
df_allsites[!is.na(df_allsites$aspect) & df_allsites$aspect > 315,]$aspect_cats <- 'North'
quantile(df_allsites$slope, na.rm=T, probs = seq(0,1,(1/3)))
df_allsites$slope_quants <- gtools::quantcut(df_allsites$slope, q = 3)
quantile(df_allsites$elev, na.rm=T, probs = seq(0,1,(1/3)))
df_allsites$elev_quants <- gtools::quantcut(df_allsites$elev, q = 3)
quantile(df_allsites$depth, na.rm=T, probs = seq(0,1,(1/5)))
df_allsites$depth_quants <- gtools::quantcut(df_allsites$depth, q = 5)
quantile(df_allsites$aws100, na.rm=T, probs = seq(0,1,(1/4)))
df_allsites$aws100_quants <- gtools::quantcut(df_allsites$aws100, q = 4)
quantile(df_allsites$sand, na.rm=T, probs = seq(0,1,(1/4)))
df_allsites$sand_quants <- gtools::quantcut(df_allsites$sand, q = 4)
quantile(df_allsites$silt, na.rm=T, probs = seq(0,1,(1/4)))
df_allsites$silt_quants <- gtools::quantcut(df_allsites$silt, q = 4)
quantile(df_allsites$clay, na.rm=T, probs = seq(0,1,(1/4)))
df_allsites$clay_quants <- gtools::quantcut(df_allsites$clay, q = 3)

# subset to remove out data with NAs for any variable
# (sand has fewest data among variables we're using, so subset by that)
df_allsites <- subset(df_allsites, !is.na(sand))

# pull out a single year of data from df_allsites to isolate the static landscape predictor variables data
xvars_df <- subset(df_allsites[,c(1:2,5,14:39)], year == '2020')
xvars_df <- xvars_df[,-c(3,16)]
xvars_df$log.slope <- log(xvars_df$slope+0.1)
xvars_df$log.heatload <- log(xvars_df$heatload)
xvars_df$log.heatload2 <- log(xvars_df$heatload2)
xvars_df$log.heatload3 <- log(xvars_df$heatload3)

# look at correlations among predictors
# for figure purposes, make new heatload variable name
xvars_df$`heat load` <- xvars_df$heatload3
load('corrchart.RData')
pdf(paste0(Rout,"/landscape_predictors.pdf"),height=5,width=5)
png(paste0(Rout,"/landscape_predictors.png"),height = 5, width = 5, units = "in", res = 600)
chart.Correlation(xvars_df[,c("heat load",'sand','depth','log.slope','northness')],
                  histogram=TRUE)
                  # pch=c(1,2,0)[factor(xvars_df$condition)],
                  # col=alpha(c("#F8766D", "#00BA38", "#619CFF")[factor(xvars_df$condition)],0.2),
                  # bg=alpha(c("#F8766D", "#00BA38", "#619CFF")[factor(xvars_df$condition)],0.2))
dev.off()

############## aggregated half-decade annual cover estimates against landscape drivers ####
df_allsites_absann_agg <- aggregate(absann ~ pixel + years, data = df_allsites, mean)
df_allsites_absann_agg <- merge(df_allsites_absann_agg, xvars_df, by = 'pixel')
df_allsites_absann_agg$years2 <- as.numeric(df_allsites_absann_agg$years)

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

#### calculate logit transformation
df_allsites_absann_agg$absann.logit <- logit(df_allsites_absann_agg$absann, adjust = 0.01)
hist(df_allsites_absann_agg$absann)
hist(df_allsites_absann_agg$absann.logit)

#### scale variables before analysis
df_allsites_absann_agg.scaled <- cbind(df_allsites_absann_agg[,c(1:4,35)],
                            scale(df_allsites_absann_agg[,c(5:18,30:34)]),
                            df_allsites_absann_agg[,19:29])

levels(as.factor(df_allsites_absann_agg$site))
######## plots faceted by site ####
df_allsites_absann_agg$site <- factor(df_allsites_absann_agg$site,
                              levels = c('Baskett Slough National Wildlife Refuge', 'Ankeny National Wildlife Refuge',
                                         'Kingston Prairie Preserve', 'Chip Ross Park', 'Bald Hill-Fitton Green Complex',
                                         'Finley National Wildlife Refuge','Rattlesnake Butte','Coburg Ridge',
                                         'Murray Hill, Upper Willow Creek','Andrew Reasoner Wildlife Preserve',
                                         'South Eugene Meadows','Suzanne Arlie, Mount Baldy','Howard Buford Recreation Area',
                                         'Thurston Hills','Native Oaks Ridge','Cerro Gordo'))
df_allsites_absann_agg_sites1 <- subset(df_allsites_absann_agg, site == 'Baskett Slough National Wildlife Refuge' |
                                          site == 'Ankeny National Wildlife Refuge' | site == 'Kingston Prairie Preserve' |
                                        site == 'Chip Ross Park' | site == 'Bald Hill-Fitton Green Complex' |
                                        site == 'Finley National Wildlife Refuge' | site == 'Rattlesnake Butte' | site == 'Coburg Ridge')
df_allsites_absann_agg_sites2 <- subset(df_allsites_absann_agg, site == 'Murray Hill, Upper Willow Creek' |
                                          site == 'Andrew Reasoner Wildlife Preserve' | site == 'South Eugene Meadows' |
                                          site == 'Suzanne Arlie, Mount Baldy' | site == 'Howard Buford Recreation Area' |
                                          site == 'Thurston Hills' | site == 'Native Oaks Ridge' | site == 'Cerro Gordo')
sites1.labs <- c("BSNWR",'ANWR','KPP','CRP','BHFG','FNWR','RB','CR')
names(sites1.labs) <- c('Baskett Slough National Wildlife Refuge', 'Ankeny National Wildlife Refuge',
                        'Kingston Prairie Preserve', 'Chip Ross Park', 'Bald Hill-Fitton Green Complex',
                        'Finley National Wildlife Refuge','Rattlesnake Butte','Coburg Ridge')
sites2.labs <- c('MHUWC','ARWP','SEM','SAMB','HBRA','TH','NOR','CG')
names(sites2.labs) <- c('Murray Hill, Upper Willow Creek','Andrew Reasoner Wildlife Preserve',
                        'South Eugene Meadows','Suzanne Arlie, Mount Baldy','Howard Buford Recreation Area',
                        'Thurston Hills','Native Oaks Ridge','Cerro Gordo')

#### heatload ####
png(paste0(Rout,"/absann_sitefacet_heat_1.png"),height = 8, width = 8, units = "in", res = 600)
ggplot(df_allsites_absann_agg_sites1, aes(x = heatload3, y = absann)) +
  # geom_hline(yintercept = 0, linetype = 2, colour = "darkgrey") +
  geom_point(aes(x = heatload3, y = absann, col = depth_quants), size = 0.2, shape = 16) +
  # geom_point(aes(x = heatload, y = absann, col = sand_quants), size = 0.2, alpha = 0.2) +
  # geom_point(aes(x = heatload, y = absann), col = 'black', size = 0.2, alpha = 0.2) +
  # scale_color_manual(values = c('#FFEB38','#FF9800','#795548')) +
  # geom_quantile(quantiles = c(0.1,0.5,0.9), col = 'black') +
  scale_x_continuous(breaks = c(0.7,0.9)) +
  scale_color_viridis_d(direction = -1,
                        name = "Soil depth", labels = c('0-41cm','41-78cm','78-104cm','104-122cm','122-202cm')) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = 'Heat load', y = "AFG (%) cover") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  facet_grid(site ~ years, scales = 'free_y', labeller = labeller(site = sites1.labs))
dev.off()

png(paste0(Rout,"/absann_sitefacet_heat_2.png"),height = 8, width = 8, units = "in", res = 600)
ggplot(df_allsites_absann_agg_sites2, aes(x = heatload3, y = absann)) +
  # geom_hline(yintercept = 0, linetype = 2, colour = "darkgrey") +
  geom_point(aes(x = heatload3, y = absann, col = depth_quants), size = 0.2, shape = 16) +
  # geom_point(aes(x = heatload, y = absann, col = sand_quants), size = 0.2, alpha = 0.2) +
  # geom_point(aes(x = heatload, y = absann), col = 'black', size = 0.2, alpha = 0.2) +
  # scale_color_manual(values = c('#FFEB38','#FF9800','#795548')) +
  # geom_smooth(method = 'lm', col = 'red', lwd = 0.5, linetype = 2) +
  # geom_quantile(quantiles = c(0.1,0.5,0.9), col = 'black') +
  scale_x_continuous(breaks = c(0.7,0.9)) +
  scale_color_viridis_d(direction = -1,
                        name = "Soil depth", labels = c('0-41cm','41-78cm','78-104cm','104-122cm','122-202cm')) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = 'Heat load', y = "AFG (%) cover") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  facet_grid(site ~ years, scales = 'free_y', labeller = labeller(site = sites2.labs))
dev.off()

 
#### northness ####
png(paste0(Rout,"/absann_sitefacet_northness_1.png"),height = 8, width = 8, units = "in", res = 600)
ggplot(df_allsites_absann_agg_sites1, aes(x = northness, y = absann)) +
  # geom_hline(yintercept = 0, linetype = 2, colour = "darkgrey") +
  geom_point(aes(x = northness, y = absann, col = depth_quants), size = 0.2, shape = 16) +
  # geom_point(aes(x = heatload, y = absann, col = sand_quants), size = 0.2, alpha = 0.2) +
  # geom_point(aes(x = heatload, y = absann), col = 'black', size = 0.2, alpha = 0.2) +
  # scale_color_manual(values = c('#FFEB38','#FF9800','#795548')) +
  # geom_smooth(method = 'lm', col = 'red', lwd = 0.5, linetype = 2) +
  # geom_quantile(quantiles = c(0.1,0.5,0.9), col = 'black') +
  scale_x_continuous(breaks = c(-1,0,1)) +
  scale_color_viridis_d(direction = -1,
                        name = "Soil depth", labels = c('0-41cm','41-78cm','78-104cm','104-122cm','122-202cm')) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = 'Northness', y = "AFG (%) cover") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  facet_grid(site ~ years, scales = 'free_y', labeller = labeller(site = sites1.labs))
dev.off()

png(paste0(Rout,"/absann_sitefacet_northness_2.png"),height = 8, width = 8, units = "in", res = 600)
ggplot(df_allsites_absann_agg_sites2, aes(x = northness, y = absann)) +
  # geom_hline(yintercept = 0, linetype = 2, colour = "darkgrey") +
  geom_point(aes(x = northness, y = absann, col = depth_quants), size = 0.2, shape = 16) +
  # geom_point(aes(x = heatload, y = absann, col = sand_quants), size = 0.2, alpha = 0.2) +
  # geom_point(aes(x = heatload, y = absann), col = 'black', size = 0.2, alpha = 0.2) +
  # scale_color_manual(values = c('#FFEB38','#FF9800','#795548')) +
  # geom_smooth(method = 'lm', col = 'red', lwd = 0.5, linetype = 2) +
  # geom_quantile(quantiles = c(0.1,0.5,0.9), col = 'black') +
  scale_x_continuous(breaks = c(-1,0,1)) +
  scale_color_viridis_d(direction = -1,
                        name = "Soil depth", labels = c('0-41cm','41-78cm','78-104cm','104-122cm','122-202cm')) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = 'Northness', y = "AFG (%) cover") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  facet_grid(site ~ years, scales = 'free_y', labeller = labeller(site = sites2.labs))
dev.off()

#### slope ####
png(paste0(Rout,"/absann_sitefacet_slope_1.png"),height = 8, width = 8, units = "in", res = 600)
ggplot(df_allsites_absann_agg_sites1, aes(x = log.slope, y = absann)) +
  # geom_hline(yintercept = 0, linetype = 2, colour = "darkgrey") +
  geom_point(aes(x = log.slope, y = absann, col = depth_quants), size = 0.2, shape = 16) +
  # geom_point(aes(x = heatload, y = absann, col = sand_quants), size = 0.2, alpha = 0.2) +
  # geom_point(aes(x = heatload, y = absann), col = 'black', size = 0.2, alpha = 0.2) +
  # scale_color_manual(values = c('#FFEB38','#FF9800','#795548')) +
  # geom_smooth(method = 'lm', col = 'red', lwd = 0.5, linetype = 2) +
  # geom_quantile(quantiles = c(0.1,0.5,0.9), col = 'black') +
  scale_color_viridis_d(direction = -1,
                        name = "Soil depth", labels = c('0-41cm','41-78cm','78-104cm','104-122cm','122-202cm')) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = 'log(Slope)', y = "AFG (%) cover") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  facet_grid(site ~ years, scales = 'free_y', labeller = labeller(site = sites1.labs))
dev.off()

png(paste0(Rout,"/absann_sitefacet_slope_2.png"),height = 8, width = 8, units = "in", res = 600)
ggplot(df_allsites_absann_agg_sites2, aes(x = log.slope, y = absann)) +
  # geom_hline(yintercept = 0, linetype = 2, colour = "darkgrey") +
  geom_point(aes(x = log.slope, y = absann, col = depth_quants), size = 0.2, shape = 16) +
  # geom_point(aes(x = heatload, y = absann, col = sand_quants), size = 0.2, alpha = 0.2) +
  # geom_point(aes(x = heatload, y = absann), col = 'black', size = 0.2, alpha = 0.2) +
  # scale_color_manual(values = c('#FFEB38','#FF9800','#795548')) +
  # geom_smooth(method = 'lm', col = 'red', lwd = 0.5, linetype = 2) +
  # geom_quantile(quantiles = c(0.1,0.5,0.9), col = 'black') +
  scale_color_viridis_d(direction = -1,
                        name = "Soil depth", labels = c('0-41cm','41-78cm','78-104cm','104-122cm','122-202cm')) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = 'log(Slope)', y = "AFG (%) cover") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  facet_grid(site ~ years, scales = 'free_y', labeller = labeller(site = sites2.labs))
dev.off()

#### aspect ####
png(paste0(Rout,"/absann_sitefacet_aspect_1.png"),height = 8, width = 8, units = "in", res = 600)
ggplot(df_allsites_absann_agg_sites1, aes(x = aspect, y = absann)) +
  # geom_hline(yintercept = 0, linetype = 2, colour = "darkgrey") +
  geom_point(aes(x = aspect, y = absann, col = depth_quants), size = 0.2, shape = 16) +
  # geom_point(aes(x = heatload, y = absann, col = sand_quants), size = 0.2, alpha = 0.2) +
  # geom_point(aes(x = heatload, y = absann), col = 'black', size = 0.2, alpha = 0.2) +
  # scale_color_manual(values = c('#FFEB38','#FF9800','#795548')) +
  # geom_smooth(method = 'lm', col = 'red', lwd = 0.5, linetype = 2) +
  # geom_quantile(quantiles = c(0.1,0.5,0.9), col = 'black') +
  scale_color_viridis_d(direction = -1,
                        name = "Soil depth", labels = c('0-41cm','41-78cm','78-104cm','104-122cm','122-202cm')) +
  scale_x_continuous(breaks = c(0,180,360)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = 'Aspect (°)', y = "AFG (%) cover") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  facet_grid(site ~ years, scales = 'free_y', labeller = labeller(site = sites1.labs))
dev.off()

png(paste0(Rout,"/absann_sitefacet_aspect_2.png"),height = 8, width = 8, units = "in", res = 600)
ggplot(df_allsites_absann_agg_sites2, aes(x = aspect, y = absann)) +
  # geom_hline(yintercept = 0, linetype = 2, colour = "darkgrey") +
  geom_point(aes(x = aspect, y = absann, col = depth_quants), size = 0.2, shape = 16) +
  # geom_point(aes(x = heatload, y = absann, col = sand_quants), size = 0.2, alpha = 0.2) +
  # geom_point(aes(x = heatload, y = absann), col = 'black', size = 0.2, alpha = 0.2) +
  # scale_color_manual(values = c('#FFEB38','#FF9800','#795548')) +
  # geom_smooth(method = 'lm', col = 'red', lwd = 0.5, linetype = 2) +
  # geom_quantile(quantiles = c(0.1,0.5,0.9), col = 'black') +
  scale_color_viridis_d(direction = -1,
                        name = "Soil depth", labels = c('0-41cm','41-78cm','78-104cm','104-122cm','122-202cm')) +
  scale_x_continuous(breaks = c(0,180,360)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = 'Aspect (°)', y = "AFG (%) cover") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  facet_grid(site ~ years, scales = 'free_y', labeller = labeller(site = sites2.labs))
dev.off()

#### depth ####
png(paste0(Rout,"/absann_sitefacet_depth_1.png"),height = 8, width = 8, units = "in", res = 600)
ggplot(df_allsites_absann_agg_sites1, aes(x = depth, y = absann)) +
  # geom_hline(yintercept = 0, linetype = 2, colour = "darkgrey") +
  geom_point(aes(x = depth, y = absann, col = heatload3_quants), size = 0.2, shape = 16) +
  # geom_point(aes(x = heatload, y = absann, col = sand_quants), size = 0.2, alpha = 0.2) +
  # geom_point(aes(x = heatload, y = absann), col = 'black', size = 0.2, alpha = 0.2) +
  # scale_color_manual(values = c('#FFEB38','#FF9800','#795548')) +
  # geom_smooth(method = 'lm', col = 'red', lwd = 0.5, linetype = 2) +
  # geom_quantile(quantiles = c(0.1,0.5,0.9), col = 'black') +
  scale_color_manual(values = c('#3B9AB2','#EBCC2A','#F21A00'), labels = c('0.59-0.91','0.91-0.93','0.93-1.02'),
                     name = "Heat load") +
  scale_x_continuous(breaks = c(0,50,150)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = 'Soil depth (cm)', y = "AFG (%) cover") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  facet_grid(site ~ years, scales = 'free_y', labeller = labeller(site = sites1.labs))
dev.off()

png(paste0(Rout,"/absann_sitefacet_depth_2.png"),height = 8, width = 8, units = "in", res = 600)
ggplot(df_allsites_absann_agg_sites2, aes(x = depth, y = absann)) +
  # geom_hline(yintercept = 0, linetype = 2, colour = "darkgrey") +
  geom_point(aes(x = depth, y = absann, col = heatload3_quants), size = 0.2, shape = 16) +
  # geom_point(aes(x = heatload, y = absann, col = sand_quants), size = 0.2, alpha = 0.2) +
  # geom_point(aes(x = heatload, y = absann), col = 'black', size = 0.2, alpha = 0.2) +
  # scale_color_manual(values = c('#FFEB38','#FF9800','#795548')) +
  # geom_smooth(method = 'lm', col = 'red', lwd = 0.5, linetype = 2) +
  # geom_quantile(quantiles = c(0.1,0.5,0.9), col = 'black') +
  scale_color_manual(values = c('#3B9AB2','#EBCC2A','#F21A00'), labels = c('0.59-0.91','0.91-0.93','0.93-1.02'),
                     name = "Heat load") +
  scale_x_continuous(breaks = c(0,50,150)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = 'Soil depth (cm)', y = "AFG (%) cover") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  facet_grid(site ~ years, scales = 'free_y', labeller = labeller(site = sites2.labs))
dev.off()

#### sand ####
png(paste0(Rout,"/absann_sitefacet_sand_1.png"),height = 8, width = 8, units = "in", res = 600)
ggplot(df_allsites_absann_agg_sites1, aes(x = sand, y = absann)) +
  # geom_hline(yintercept = 0, linetype = 2, colour = "darkgrey") +
  geom_point(aes(x = sand, y = absann, col = heatload3_quants), size = 0.2, shape = 16) +
  # geom_point(aes(x = heatload, y = absann, col = sand_quants), size = 0.2, alpha = 0.2) +
  # geom_point(aes(x = heatload, y = absann), col = 'black', size = 0.2, alpha = 0.2) +
  # scale_color_manual(values = c('#FFEB38','#FF9800','#795548')) +
  # geom_smooth(method = 'lm', col = 'red', lwd = 0.5, linetype = 2) +
  # geom_quantile(quantiles = c(0.1,0.5,0.9), col = 'black') +
  scale_color_manual(values = c('#3B9AB2','#EBCC2A','#F21A00'), labels = c('0.59-0.91','0.91-0.93','0.93-1.02'),
                     name = "Heat load") +
  # scale_x_continuous(breaks = c(0,50,150)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = '(%) Sand', y = "AFG (%) cover") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  facet_grid(site ~ years, scales = 'free_y', labeller = labeller(site = sites1.labs))
dev.off()

png(paste0(Rout,"/absann_sitefacet_sand_2.png"),height = 8, width = 8, units = "in", res = 600)
ggplot(df_allsites_absann_agg_sites2, aes(x = sand, y = absann)) +
  # geom_hline(yintercept = 0, linetype = 2, colour = "darkgrey") +
  geom_point(aes(x = sand, y = absann, col = heatload3_quants), size = 0.2, shape = 16) +
  # geom_point(aes(x = heatload, y = absann, col = sand_quants), size = 0.2, alpha = 0.2) +
  # geom_point(aes(x = heatload, y = absann), col = 'black', size = 0.2, alpha = 0.2) +
  # scale_color_manual(values = c('#FFEB38','#FF9800','#795548')) +
  # geom_smooth(method = 'lm', col = 'red', lwd = 0.5, linetype = 2) +
  # geom_quantile(quantiles = c(0.1,0.5,0.9), col = 'black') +
  scale_color_manual(values = c('#3B9AB2','#EBCC2A','#F21A00'), labels = c('0.59-0.91','0.91-0.93','0.93-1.02'),
                     name = "Heat load") +
  # scale_x_continuous(breaks = c(0,50,150)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = '(%) Sand', y = "AFG (%) cover") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  facet_grid(site ~ years, scales = 'free_y', labeller = labeller(site = sites2.labs))
dev.off()

######## heatload rq model ####
m.allsites.absann.heat <- rq(absann.logit ~ years2 * (log.heatload + depth + sand), data = df_allsites_absann_agg.scaled, tau = c(0.1,0.5,0.9))
m.allsites.absann.heat2 <- rq(absann.logit ~ years2 * (heatload2 + depth + sand), data = df_allsites_absann_agg.scaled, tau = c(0.1,0.5,0.9))
m.allsites.absann.heat3 <- rq(absann.logit ~ years2 * (heatload3 + depth + sand), data = df_allsites_absann_agg.scaled, tau = c(0.1,0.5,0.9))

######## northness + slope rq model ####
m.allsites.absann.north <- rq(absann.logit ~ years2 * (northness + log.slope + depth + sand), data = df_allsites_absann_agg.scaled, tau = c(0.1,0.5,0.9))

######## compare AIC between heatload and northness ####
AIC.heat <- AIC(m.allsites.absann.heat)
AIC.heat2 <- AIC(m.allsites.absann.heat2)
AIC.heat3 <- AIC(m.allsites.absann.heat3)
AIC.north <- AIC(m.allsites.absann.north)
# heatload3 model is generally the best! that's what I was hoping/expecting. 

# calculate R2
m.allsites.absann0 <- rq(absann.logit ~ 1, data = df_allsites_absann_agg.scaled, tau = c(0.1,0.5,0.9))
R2.allsites.absann.heat <- 1 - m.allsites.absann.heat$rho/m.allsites.absann0$rho
R2.allsites.absann.heat2 <- 1 - m.allsites.absann.heat2$rho/m.allsites.absann0$rho
R2.allsites.absann.heat3 <- 1 - m.allsites.absann.heat3$rho/m.allsites.absann0$rho
R2.allsites.absann.north <- 1 - m.allsites.absann.north$rho/m.allsites.absann0$rho

# need to fit separate quantregs to make predictions.
m.allsites.absann.1 <- rq(absann.logit ~ years2 * (heatload3 + depth + sand), data = df_allsites_absann_agg.scaled, na.action = na.fail, tau = 0.1)
m.allsites.absann.5 <- rq(absann.logit ~ years2 * (heatload3 + depth + sand), data = df_allsites_absann_agg.scaled, na.action = na.fail, tau = 0.5)
m.allsites.absann.9 <- rq(absann.logit ~ years2 * (heatload3 + depth + sand), data = df_allsites_absann_agg.scaled, na.action = na.fail, tau = 0.9)
summ(m.allsites.absann.9)

pdf(paste0(Rout,"/allsites_coeffplot_heatload_absann3.pdf"),height = 3, width = 5)
plot_summs(m.allsites.absann.9, m.allsites.absann.5, m.allsites.absann.1, model.names = c('0.9','0.5','0.1'), colors = c("grey10",'grey50','grey80'), legend.title = 'Quantile', point.size = 1)
dev.off()

# save results of model
summary(m.allsites.absann.heat3) # screenshot this for pasting into results
anova(m.allsites.absann.heat3) # screenshot this for pasting into results
write.table(R2.allsites.absann.heat3, "clipboard", sep="\t", row.names = FALSE, col.names = FALSE)


######## predictions from model (for heatload) ####
# predictions focusing on years2 and heatload, holding depth and sand at their mean (which is 0 when scaled)
nd = expand.grid(years2 = seq(min(df_allsites_absann_agg.scaled$years2), max(df_allsites_absann_agg.scaled$years2), length.out = 7),
                 heatload3 = seq(min(df_allsites_absann_agg.scaled$heatload3), max(df_allsites_absann_agg.scaled$heatload3), length.out = 20),
                 depth = 0,
                 sand = 0)

# make nd.raw to tack on for plotting purposes
nd.raw = expand.grid(years = levels(df_allsites_absann_agg$years),
                     heatload3 = seq(min(df_allsites_absann_agg$heatload3), max(df_allsites_absann_agg$heatload3), length.out = 20))

# in order to get confidence intervals around rq predictions, you must use models fitted separately for each level of tau
p.allsites.absann.1 = predict(m.allsites.absann.1,
                             newdata = nd, interval = 'confidence')
p.allsites.absann.5 = predict(m.allsites.absann.5,
                             newdata = nd, interval = 'confidence')
p.allsites.absann.9 = predict(m.allsites.absann.9,
                             newdata = nd, interval = 'confidence')

# now tack the predictions onto a dataframe with the RAW values
p.allsites.absann.fit <- data.frame(nd.raw, 'tau0.1' = p.allsites.absann.1[,1], 'tau0.5' = p.allsites.absann.5[,1], 'tau0.9' = p.allsites.absann.9[,1])
p.allsites.absann.lower <- data.frame(nd.raw, 'tau0.1' = p.allsites.absann.1[,2], 'tau0.5' = p.allsites.absann.5[,2], 'tau0.9' = p.allsites.absann.9[,2])
p.allsites.absann.upper <- data.frame(nd.raw, 'tau0.1' = p.allsites.absann.1[,3], 'tau0.5' = p.allsites.absann.5[,3], 'tau0.9' = p.allsites.absann.9[,3])
p.allsites.absann <- gather(p.allsites.absann.fit, Quantile, fit, tau0.1:tau0.9, factor_key=TRUE)
p.allsites.absann <- cbind(p.allsites.absann,
                          lower = gather(p.allsites.absann.lower, Quantile, lower, tau0.1:tau0.9, factor_key=TRUE)[,4])
p.allsites.absann <- cbind(p.allsites.absann,
                          upper = gather(p.allsites.absann.upper, Quantile, upper, tau0.1:tau0.9, factor_key=TRUE)[,4])
p.allsites.absann$Quantile <- gsub("tau","",as.character(p.allsites.absann$Quantile))

## back-transform from logit scale
p.allsites.absann$fit.back <- logit2prob(p.allsites.absann$fit) *100
p.allsites.absann$lower.back <- logit2prob(p.allsites.absann$lower) *100
p.allsites.absann$upper.back <- logit2prob(p.allsites.absann$upper) *100

p.allsites.absann.1 <- subset(p.allsites.absann, Quantile == '0.1')
p.allsites.absann.5 <- subset(p.allsites.absann, Quantile == '0.5')
p.allsites.absann.9 <- subset(p.allsites.absann, Quantile == '0.9')

p.heat.depth <- ggplot() +
                  # geom_hline(yintercept = 0, linetype = 2, colour = "darkgrey") +
                  geom_point(data = df_allsites_absann_agg, aes(x = heatload3, y = absann, col = depth_quants), size = 0.2, shape = 16) +
                  geom_ribbon(data = p.allsites.absann, aes(x = heatload3, ymin = lower.back, ymax = lower.back, group = Quantile), col = NA, alpha = 0.2) +
                  geom_line(data = p.allsites.absann, aes(x = heatload3, y = fit.back, group = Quantile)) +
                  # scale_color_brewer(palette = 'OrRd') +
                  # scale_color_manual(values = c("#619CFF","#00BA38","#F8766D")) +
                  # scale_shape_manual(values = c(3,16)) +
                  scale_x_continuous(breaks = c(0.7,0.9)) +
                  # scale_color_manual(values = rev(wes_palette('Zissou1',n = 5)),
                  #                    name = "Soil depth", labels = c('0-41cm','41-78cm','78-104cm','104-122cm','122-202cm')) +
                  scale_color_viridis_d(direction = -1,
                                     name = "Soil depth", labels = c('0-41cm','41-78cm','78-104cm','104-122cm','122-202cm')) +
                  theme_bw() +
                  theme(panel.grid = element_blank()) +
                  labs(x = 'Heat load', y = "AFG\n(%) Cover") +
                  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1)),
                         shape = guide_legend(override.aes = list(alpha = 1, size = 1))) +
                  facet_grid(. ~ years)
png(paste0(Rout,"/absann_preds_heatload_depth.png"),height = 1.6, width = 8, units = "in", res = 1200)
p.heat.depth
dev.off()

######## predictions from model (for depth) ####
# predictions focusing on years2 and depth, holding heatload and sand at their mean (which is 0 when scaled)
nd = expand.grid(years2 = seq(min(df_allsites_absann_agg.scaled$years2), max(df_allsites_absann_agg.scaled$years2), length.out = 7),
                 depth = seq(min(df_allsites_absann_agg.scaled$depth), max(df_allsites_absann_agg.scaled$depth), length.out = 20),
                 heatload3 = 0,
                 sand = 0)

# make nd.raw to tack on for plotting purposes
nd.raw = expand.grid(years = levels(df_allsites_absann_agg$years),
                     depth = seq(min(df_allsites_absann_agg$depth), max(df_allsites_absann_agg$depth), length.out = 20))

# in order to get confidence intervals around rq predictions, you must use models fitted separately for each level of tau
p.allsites.absann.1 = predict(m.allsites.absann.1,
                              newdata = nd, interval = 'confidence')
p.allsites.absann.5 = predict(m.allsites.absann.5,
                              newdata = nd, interval = 'confidence')
p.allsites.absann.9 = predict(m.allsites.absann.9,
                              newdata = nd, interval = 'confidence')

# now tack the predictions onto a dataframe with the RAW values
p.allsites.absann.fit <- data.frame(nd.raw, 'tau0.1' = p.allsites.absann.1[,1], 'tau0.5' = p.allsites.absann.5[,1], 'tau0.9' = p.allsites.absann.9[,1])
p.allsites.absann.lower <- data.frame(nd.raw, 'tau0.1' = p.allsites.absann.1[,2], 'tau0.5' = p.allsites.absann.5[,2], 'tau0.9' = p.allsites.absann.9[,2])
p.allsites.absann.upper <- data.frame(nd.raw, 'tau0.1' = p.allsites.absann.1[,3], 'tau0.5' = p.allsites.absann.5[,3], 'tau0.9' = p.allsites.absann.9[,3])
p.allsites.absann <- gather(p.allsites.absann.fit, Quantile, fit, tau0.1:tau0.9, factor_key=TRUE)
p.allsites.absann <- cbind(p.allsites.absann,
                           lower = gather(p.allsites.absann.lower, Quantile, lower, tau0.1:tau0.9, factor_key=TRUE)[,4])
p.allsites.absann <- cbind(p.allsites.absann,
                           upper = gather(p.allsites.absann.upper, Quantile, upper, tau0.1:tau0.9, factor_key=TRUE)[,4])
p.allsites.absann$Quantile <- gsub("tau","",as.character(p.allsites.absann$Quantile))

## back-transform from logit scale
p.allsites.absann$fit.back <- logit2prob(p.allsites.absann$fit) *100
p.allsites.absann$lower.back <- logit2prob(p.allsites.absann$lower) *100
p.allsites.absann$upper.back <- logit2prob(p.allsites.absann$upper) *100

p.allsites.absann.1 <- subset(p.allsites.absann, Quantile == '0.1')
p.allsites.absann.5 <- subset(p.allsites.absann, Quantile == '0.5')
p.allsites.absann.9 <- subset(p.allsites.absann, Quantile == '0.9')

p.depth.heat <- ggplot() +
                  # geom_hline(yintercept = 0, linetype = 2, colour = "darkgrey") +
                  geom_point(data = df_allsites_absann_agg, aes(x = depth, y = absann, col = heatload3_quants), size = 0.2, shape = 16) +
                  geom_ribbon(data = p.allsites.absann, aes(x = depth, ymin = lower.back, ymax = lower.back, group = Quantile), col = NA, alpha = 0.2) +
                  geom_line(data = p.allsites.absann, aes(x = depth, y = fit.back, group = Quantile)) +
                  # scale_color_brewer(palette = 'OrRd') +
                  # scale_color_manual(values = c("#619CFF","#00BA38","#F8766D")) +
                  # scale_shape_manual(values = c(3,16)) +
                  scale_x_continuous(breaks = c(0,50,150)) +
                  # scale_color_manual(values = (wes_palette('Zissou1', n = 3)),
                  #                    name = "Heat load") +
                  scale_color_manual(values = c('#3B9AB2','#EBCC2A','#F21A00'), labels = c('0.59-0.91','0.91-0.93','0.93-1.02'),
                                     name = "Heat load") +
                  theme_bw() +
                  theme(panel.grid = element_blank()) +
                  labs(x = 'Soil depth (cm)', y = "AFG\n(%) Cover") +
                  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1)),
                         shape = guide_legend(override.aes = list(alpha = 1, size = 1))) +
                  facet_grid(. ~ years)
png(paste0(Rout,"/absann_preds_depth_heatload.png"),height = 1.6, width = 8, units = "in", res = 1200)
p.depth.heat
dev.off()

######## predictions from model (for sand) ####
# predictions focusing on years2 and sand, holding northness, slope, and sand at their mean (which is 0 when scaled)
nd = expand.grid(years2 = seq(min(df_allsites_absann_agg.scaled$years2), max(df_allsites_absann_agg.scaled$years2), length.out = 7),
                 sand = seq(min(df_allsites_absann_agg.scaled$sand), max(df_allsites_absann_agg.scaled$sand), length.out = 20),
                 heatload3 = 0,
                 depth = 0)

# make nd.raw to tack on for plotting purposes
nd.raw = expand.grid(years = levels(df_allsites_absann_agg$years),
                     sand = seq(min(df_allsites_absann_agg$sand), max(df_allsites_absann_agg$sand), length.out = 20))

# in order to get confidence intervals around rq predictions, you must use models fitted separately for each level of tau
p.allsites.absann.1 = predict(m.allsites.absann.1,
                              newdata = nd, interval = 'confidence')
p.allsites.absann.5 = predict(m.allsites.absann.5,
                              newdata = nd, interval = 'confidence')
p.allsites.absann.9 = predict(m.allsites.absann.9,
                              newdata = nd, interval = 'confidence')

# now tack the predictions onto a dataframe with the RAW values
p.allsites.absann.fit <- data.frame(nd.raw, 'tau0.1' = p.allsites.absann.1[,1], 'tau0.5' = p.allsites.absann.5[,1], 'tau0.9' = p.allsites.absann.9[,1])
p.allsites.absann.lower <- data.frame(nd.raw, 'tau0.1' = p.allsites.absann.1[,2], 'tau0.5' = p.allsites.absann.5[,2], 'tau0.9' = p.allsites.absann.9[,2])
p.allsites.absann.upper <- data.frame(nd.raw, 'tau0.1' = p.allsites.absann.1[,3], 'tau0.5' = p.allsites.absann.5[,3], 'tau0.9' = p.allsites.absann.9[,3])
p.allsites.absann <- gather(p.allsites.absann.fit, Quantile, fit, tau0.1:tau0.9, factor_key=TRUE)
p.allsites.absann <- cbind(p.allsites.absann,
                           lower = gather(p.allsites.absann.lower, Quantile, lower, tau0.1:tau0.9, factor_key=TRUE)[,4])
p.allsites.absann <- cbind(p.allsites.absann,
                           upper = gather(p.allsites.absann.upper, Quantile, upper, tau0.1:tau0.9, factor_key=TRUE)[,4])
p.allsites.absann$Quantile <- gsub("tau","",as.character(p.allsites.absann$Quantile))

## back-transform from logit scale
p.allsites.absann$fit.back <- logit2prob(p.allsites.absann$fit) *100
p.allsites.absann$lower.back <- logit2prob(p.allsites.absann$lower) *100
p.allsites.absann$upper.back <- logit2prob(p.allsites.absann$upper) *100

p.sand.heat <- ggplot() +
                  # geom_hline(yintercept = 0, linetype = 2, colour = "darkgrey") +
                  geom_point(data = df_allsites_absann_agg, aes(x = sand, y = absann, col = heatload3_quants), size = 0.2, shape = 16) +
                  geom_ribbon(data = p.allsites.absann, aes(x = sand, ymin = lower.back, ymax = upper.back, group = Quantile), col = NA, alpha = 0.2) +
                  geom_line(data = p.allsites.absann, aes(x = sand, y = fit.back, group = Quantile)) +
                  # scale_color_brewer(palette = 'OrRd') +
                  # scale_color_manual(values = c("#619CFF","#00BA38","#F8766D")) +
                  # scale_shape_manual(values = c(3,16)) +
                  # scale_x_continuous(breaks = c(0,50,150)) +
                  # scale_color_manual(values = (wes_palette('Zissou1', n = 3)),
                  #                    name = "Heat load") +
                  scale_color_manual(values = c('#3B9AB2','#EBCC2A','#F21A00'), labels = c('0.59-0.91','0.91-0.93','0.93-1.02'),
                                     name = "Heat load") +
                  theme_bw() +
                  theme(panel.grid = element_blank()) +
                  labs(x = '(%) Sand', y = "AFG\n(%) Cover") +
                  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1)),
                         shape = guide_legend(override.aes = list(alpha = 1, size = 1))) +
                  facet_grid(. ~ years)
png(paste0(Rout,"/absann_preds_sand_heatload.png"),height = 1.6, width = 8, units = "in", res = 1200)
p.sand.heat
dev.off()

#### cowplot figure ####
library(cowplot)

# remake figures without legend
noleg_p.heat.depth <- p.heat.depth + theme(legend.position = 'none')
noleg_p.depth.heat <- p.depth.heat + theme(legend.position = 'none')
noleg_p.sand.heat <- p.sand.heat + theme(legend.position = 'none')

# extract the legends from original
lgd1<-get_legend(p.heat.depth)
lgd2<-get_legend(p.depth.heat)

# Arrange the plots and legends in a grid, with the 3 plots taking the left column,
# and the legends taking the right column (scale relative widths so the right column is only 10% of the left)
col1 <- plot_grid(noleg_p.heat.depth, noleg_p.depth.heat, noleg_p.sand.heat,
                labels = c("A","B","C"), align = "h",axis = "tb", ncol = 1)
col1
# row1<-plot_grid(plot_adtsurv,plot_sdlsurv,labels=c("A","B"),align="h",axis="tb")

col2 <- plot_grid(lgd1,lgd2, ncol = 1, rel_heights = c(0.5,1))
col2

fig <- plot_grid(col1, col2, ncol = 2,
                         rel_widths = c(1,0.15))

png(paste0(Rout,"/absann_preds_complete2.png"),height = 5, width = 8.5, units = "in", res = 600)
fig
dev.off()
