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


# in these analyses, we will be taking the rate of change raster layer and calculating summary statistics to report.

# set working directory
setwd("C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/Other projects/RAP analysis/R_working_directory")

# set directory for figure exports
Rout <- 'C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/Other projects/RAP analysis/R_working_directory/R figure exports/allsites'

# load in sites shapefile
sites <- readOGR(dsn = "C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/Other projects/RAP analysis/ArcMap raster exports/sites16shape.shp", layer="sites16shape")
plot(sites)

#### load in rate of change TIFF file, clip to sites ####
# load in processed rate of change layer
ratechange <- brick('ratechange_AbsAnn_processed.tif')
names(ratechange) <- c('Coeffs.','P-Val','R-squared','Signif. coeffs.')

# clip ratechange layers to sites area
ratechange <- crop(ratechange, extent(sites))
ratechange <- mask(ratechange, sites)
plot(ratechange[[1]])

# convert to dataframe
ratechange_df <- as.data.frame(ratechange, xy = TRUE)

# remove all NAs for coeffs column
ratechange_df <- subset(ratechange_df, !is.na(Coeffs.))

# create xy factor to use as a random effect
ratechange_df$pixel <- factor(paste(ratechange_df$x, ratechange_df$y, sep = ','))

# bring in df_allsites data so we can merge sitenames, xvars data
df_allsites <- readRDS(file = "df_polygons.rds")
df_allsites$heatload_quants <- gtools::quantcut(df_allsites$heatload, q = 3)
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

# narrow down to the only things we need from this df to merge
df_allsites2 <- subset(df_allsites, year == '2020')
df_allsites2 <- df_allsites2[,c(1:4,13:25,27:35)]

# for some reason, it seems to be duplicating some rows at Bald Hill... make sure to remove duplicates
df_allsites2 <- unique(df_allsites2)

# merge to ratechange_df
ratechange_df2 <- merge(ratechange_df,df_allsites2, by = 'pixel')
# some pixel IDs are not matching up perfectly... recreate a new pixelID column

# first, sort dataframes by x,y
df_allsites2 <- df_allsites2[order(df_allsites2$x, df_allsites2$y),]
ratechange_df <- ratechange_df[order(ratechange_df$x, ratechange_df$y),]

# now, make new pixelID columns
df_allsites2$pixelID <- as.factor(seq_len(nrow(df_allsites2)))
ratechange_df$pixelID <- as.factor(seq_len(nrow(ratechange_df)))

# now redo merge; should work
ratechange_df2 <- merge(ratechange_df,df_allsites2, by = 'pixelID')

# assign whether a pixel annualized, lost AFG cover, or did not change
ratechange_df2$status <- 'ns'
ratechange_df2[!is.na(ratechange_df2$Signif..coeffs.) & ratechange_df2$Signif..coeffs. < 0,]$status <- 'AFG declined'
ratechange_df2[!is.na(ratechange_df2$Signif..coeffs.) & ratechange_df2$Signif..coeffs. > 0,]$status <- 'AFG increased'
ratechange_gained <- subset(ratechange_df2, status == 'AFG increased')
ratechange_lost <- subset(ratechange_df2, status == 'AFG declined')
ratechange_ns <- subset(ratechange_df2, status == 'ns')

# get statistics for different groups
median(ratechange_df2$Coeffs.)
median(ratechange_gained$Coeffs.)
median(ratechange_lost$Coeffs.)
median(ratechange_ns$Coeffs.)

median(ratechange_df2$R.squared)
median(ratechange_gained$R.squared)
median(ratechange_lost$R.squared)
median(ratechange_ns$R.squared)

# calculate medians by status
medians <- aggregate(Coeffs. ~ status, data = ratechange_df2, FUN = 'median')
medians <- subset(medians, !status == 'ns')

# calculate medians by status and site
medians.bysite <- aggregate(Coeffs. ~ status + site, data = ratechange_df2, FUN = 'median')
medians.bysite <- subset(medians.bysite, !status == 'ns')
median.bysite <- aggregate(Coeffs. ~ site, data = ratechange_df2, FUN = 'median')

# # calculate quantiles
# quants <- quantile(ratechange_df2$Coeffs., probs = c(0.1,0.5,0.9))
# medians <- subset(medians, !status == 'ns')
# 
# # calculate quantiles by site
# quants.bysite <- ratechange_df2 %>% group_by(site) %>%
#   summarize(res=quantile(x,probs=0.5))
# library(plyr)
# quants.bysite <- do.call("rbind", tapply(ratechange_df2$Coeffs., ratechange_df2$site, quantile, c(0.1, 0.5, 0.9)))


# overall histogram
ratechange_df2$status <- factor(ratechange_df2$status, levels = c('ns','AFG declined','AFG increased'))

pdf(paste0(Rout,"/ratechange_hist.pdf"),height = 3, width = 3)
ggplot(data = ratechange_df2,aes(x = Coeffs., fill = status)) +
  geom_histogram(alpha = 0.75, position  =  'identity') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed') +
  geom_vline(data = medians, aes(xintercept = Coeffs., color = status)) +
  geom_vline(aes(xintercept = median(ratechange_df2$Coeffs.))) +
  scale_fill_manual(values = c('darkgrey','#4F2C7F','#C06728')) +
  scale_color_manual(values = c('#4F2C7F','#C06728')) +
  theme_bw() +
  guides(color = 'none') +
  theme(legend.position = 'bottom',panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(fill = "",x = 'AFG cover year-1',y = 'Number of pixels')
dev.off()

# histograms within sites
ratechange_df2$site <- factor(ratechange_df2$site,
                              levels = c('Baskett Slough National Wildlife Refuge', 'Ankeny National Wildlife Refuge',
                                         'Kingston Prairie Preserve', 'Chip Ross Park', 'Bald Hill-Fitton Green Complex',
                                         'Finley National Wildlife Refuge','Rattlesnake Butte','Coburg Ridge',
                                         'Murray Hill, Upper Willow Creek','Andrew Reasoner Wildlife Preserve',
                                         'South Eugene Meadows','Suzanne Arlie, Mount Baldy','Howard Buford Recreation Area',
                                         'Thurston Hills','Native Oaks Ridge','Cerro Gordo'))

pdf(paste0(Rout,"/ratechange_hists_bysite.pdf"),height = 7, width = 7)
ggplot(data = ratechange_df2,aes(x = Coeffs., fill = status)) +
  geom_histogram(alpha = 0.75, position  =  'identity',bins = 15) +
  geom_vline(aes(xintercept = 0), linetype = 'dashed')+
  geom_vline(data = medians.bysite, aes(xintercept = Coeffs., color = status)) +
  geom_vline(data = median.bysite, aes(xintercept = Coeffs.)) +
  scale_fill_manual(values = c('darkgrey','#4F2C7F','#C06728')) +
  scale_color_manual(values = c('#4F2C7F','#C06728')) +
  theme_bw() + 
  guides(color = 'none') +
  theme(legend.position = 'bottom',panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(fill = "",x = 'AFG cover year-1',y = 'Number of pixels') +
  facet_wrap(site~., scales = 'free')
dev.off()
