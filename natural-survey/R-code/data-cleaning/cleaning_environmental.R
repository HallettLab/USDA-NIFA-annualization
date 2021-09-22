# clean natural site survey environmental datasets
# created: 2021-09-22
library(matrixStats)

# script purpose:
# read in environmental data and calculate mean/var of soil moisture and depth across 5 subsamples per plot.
# write out clean environmental dataset

# set working directory to access .csv files
setwd("C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/Natural Site Survey/Reed_USDA-NIFA_NaturalSurvey/Data-entered")

# set path to cleaned data folder
cl <- "C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/Natural Site Survey/Reed_USDA-NIFA_NaturalSurvey/Data-cleaned"

# read in data
env <- read.csv("environmental.csv",header=TRUE)

# calculate mean and variance for soil moisture and depth across the 5 subsamples per plot
env2 <- env[,1:9]

env2$moist_mean <- rowMeans(env[,10:14])
env2$moist_var <- rowVars(as.matrix(env[,10:14]))
env2$moist_sd <- rowSds(as.matrix(env[,10:14]))

env2$depth_mean <- rowMeans(env[,15:19])
env2$depth_var <- rowVars(as.matrix(env[,15:19]))
env2$depth_sd <- rowSds(as.matrix(env[,15:19]))

# save data to cleaned file
write.csv(env2,file=paste0(cl,"/environmental.csv"),row.names = FALSE)

