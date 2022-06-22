library(raster)
library(rgdal)
# library(ggplot2)
library(tidyr)
library(dplyr)
library(car)

# set working directory
setwd("C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/RAP analysis/R_working_directory")

# prep for loading each yearly file at once into a list
years <- 1986:2020

#### load in processed RAP anchor layers ####
layers_processed_anchors <- list()
for (i in 1:length(years)) {
  filename <- paste0('RAP_processed_anchors_',years[i],'.tif')
  layer <- brick(filename)
  names(layer) <- c(paste0('Annual_',years[i]), paste0('Perennial_',years[i]), paste0('Total_',years[i]),
                    paste0('Rel_Ann_',years[i]), paste0('Rel_Peren_',years[i]), paste0('LogRatio_',years[i]))
  layers_processed_anchors <- c(layers_processed_anchors, list(layer))
  rm(filename); rm(layer) # remove temporary objects to save memory
}


########## ABSOLUTE ANNUAL COVER ##########
#### convert AbsAnn layers to matrices then dataframes; make rownames into y value ####
matrices_AbsAnn <- list()
for (i in 1:length(years)) {
  mtrx <- data.frame( as.matrix( layers_processed_anchors[[i]][[1]]) )
  mtrx$y <- rownames(mtrx)
  matrices_AbsAnn <- c(matrices_AbsAnn, list(mtrx))
  # rm(mtrx) # remove temporary object to save memory
}


#### convert wide to long ####
matrices_AbsAnn.long <- list()
for (i in 1:length(years)) {
  mtrx.long <- gather(matrices_AbsAnn[[i]], x, pixelvalue, X1:X3115, factor_key=TRUE, na.rm = T)
  mtrx.long$year <- years[i]
  matrices_AbsAnn.long <- c(matrices_AbsAnn.long, list(mtrx.long))
  # rm(mtrx.long) # remove temporary object to save memory
}



#### stack all years into one super long dataframe ####
AbsAnn.long <- bind_rows(matrices_AbsAnn.long)
AbsAnn.long$xy <- paste(AbsAnn.long$x,AbsAnn.long$y,sep=',')

#### prep for looping regressions across every pixel ####
AbsAnn.long$xy <- factor(AbsAnn.long$xy)
xy <- levels(AbsAnn.long$xy)
hist(AbsAnn.long$pixelvalue)
hist(logit(AbsAnn.long$pixelvalue, adjust = 0.01))
hist(AbsAnn.long[AbsAnn.long$xy == 'X2,1291',]$pixelvalue)


#### loop models over every pixel ####
# using by and apply functions is exponentially faster than for-loop
Mods <- with(AbsAnn.long,
             by(AbsAnn.long, xy, function(x) lm(pixelvalue ~ as.numeric(year), data = x)))
Summaries <- lapply(Mods, summary)
Coeffs <- sapply(Summaries, function(x) x$coefficients[2,1])
Ps <- sapply(Summaries, function(x) x$coefficients[2,4])
Rs <- sapply(Summaries, function(x) r_sq = x$r.squared)
# rm(Mods) # help regain some memory


#### steal dataframe structure to have something to put coefficient value into ####
rasterframe <- mtrx.long[,1:3]
rasterframe$xy <- factor(paste(rasterframe$x,rasterframe$y,sep=','))
rasterframe <- rasterframe[order(rasterframe$xy),]
rasterframe$Coeff <- Coeffs
rasterframe$Pval <- Ps
rasterframe$Rsquared <- Rs
rasterframe$Coeff_Signif <- NA

# a few pixels had 0 annual cover for the entire 35 years, so those have NaN for P-value and R-squared (model wasn't able to compute)
# Change those to 1 and 0, respectively.
rasterframe[is.nan(rasterframe$Pval),]$Pval <- 1
rasterframe[is.nan(rasterframe$Rsquared),]$Rsquared <- 0

rasterframe[rasterframe$Pval < 0.05,]$Coeff_Signif <- rasterframe[rasterframe$Pval < 0.05,]$Coeff
rasterframe$Rsquared_Signif <- NA
rasterframe[rasterframe$Pval < 0.05,]$Rsquared_Signif <- rasterframe[rasterframe$Pval < 0.05,]$Rsquared

# now create a dataframe with only NAs so you can bind it all back together
rasterframe.NA <- gather(matrices_AbsAnn[[i]], x, pixelvalue, X1:X3115, factor_key=TRUE)
rasterframe.NA <- rasterframe.NA[is.na(rasterframe.NA$pixelvalue),]
# rasterframe.NA$xy <- factor(paste(rasterframe.NA$x,rasterframe.NA$y,sep=','))
rasterframe.NA$Coeff <- NA
rasterframe.NA$Pval <- NA
rasterframe.NA$Rsquared <- NA
rasterframe.NA$Coeff_Signif <- NA
rasterframe.NA$Rsquared_Signif <- NA

# bind the two together
rasterframe <- rasterframe[,-c(3:4)] # throw away columns you don't need
rasterframe.NA <- rasterframe.NA[,-3] # throw away column you don't need
rasterframe2 <- rbind(rasterframe, rasterframe.NA)
rasterframe2 <- rasterframe2[order(rasterframe2$x, as.numeric(rasterframe2$y)),]


#### spread new dataframe back out wide; do separately for Coeff, Pval, Rsquared, Coeff_Signif, Rsquared_Signif ####
rasterframe2.Coeff <- rasterframe2[,-c(4:7)]
rasterframe2.Pval <- rasterframe2[,-c(3,5:7)]
rasterframe2.Rsquared <- rasterframe2[,-c(3:4,6:7)]
rasterframe2.Coeff_Signif <- rasterframe2[,-c(3:5,7)]
rasterframe2.Rsquared_Signif <- rasterframe2[,-c(3:6)]

rasterframe2.Coeff.wide <- spread(rasterframe2.Coeff, x, Coeff)
rasterframe2.Coeff.wide <- rasterframe2.Coeff.wide[order(as.numeric(rasterframe2.Coeff.wide$y)),]
rownames(rasterframe2.Coeff.wide) <- rasterframe2.Coeff.wide[,1]
rasterframe2.Coeff.wide <- rasterframe2.Coeff.wide[,-1]

rasterframe2.Pval.wide <- spread(rasterframe2.Pval, x, Pval)
rasterframe2.Pval.wide <- rasterframe2.Pval.wide[order(as.numeric(rasterframe2.Pval.wide$y)),]
rownames(rasterframe2.Pval.wide) <- rasterframe2.Pval.wide[,1]
rasterframe2.Pval.wide <- rasterframe2.Pval.wide[,-1]

rasterframe2.Rsquared.wide <- spread(rasterframe2.Rsquared, x, Rsquared)
rasterframe2.Rsquared.wide <- rasterframe2.Rsquared.wide[order(as.numeric(rasterframe2.Rsquared.wide$y)),]
rownames(rasterframe2.Rsquared.wide) <- rasterframe2.Rsquared.wide[,1]
rasterframe2.Rsquared.wide <- rasterframe2.Rsquared.wide[,-1]

rasterframe2.Coeff_Signif.wide <- spread(rasterframe2.Coeff_Signif, x, Coeff_Signif)
rasterframe2.Coeff_Signif.wide <- rasterframe2.Coeff_Signif.wide[order(as.numeric(rasterframe2.Coeff_Signif.wide$y)),]
rownames(rasterframe2.Coeff_Signif.wide) <- rasterframe2.Coeff_Signif.wide[,1]
rasterframe2.Coeff_Signif.wide <- rasterframe2.Coeff_Signif.wide[,-1]

rasterframe2.Rsquared_Signif.wide <- spread(rasterframe2.Rsquared_Signif, x, Rsquared_Signif)
rasterframe2.Rsquared_Signif.wide <- rasterframe2.Rsquared_Signif.wide[order(as.numeric(rasterframe2.Rsquared_Signif.wide$y)),]
rownames(rasterframe2.Rsquared_Signif.wide) <- rasterframe2.Rsquared_Signif.wide[,1]
rasterframe2.Rsquared_Signif.wide <- rasterframe2.Rsquared_Signif.wide[,-1]

#### finally, convert back to matrices and then rasters, using template from starting raster to get WGS 1984
ratechange_AbsAnn_Coeff <- as.matrix(rasterframe2.Coeff.wide)
ratechange_AbsAnn_Pval <- as.matrix(rasterframe2.Pval.wide)
ratechange_AbsAnn_Rsquared <- as.matrix(rasterframe2.Rsquared.wide)
ratechange_AbsAnn_Coeff_Signif <- as.matrix(rasterframe2.Coeff_Signif.wide)
ratechange_AbsAnn_Rsquared_Signif <- as.matrix(rasterframe2.Rsquared_Signif.wide)

writeRaster(
  stack(raster(ratechange_AbsAnn_Coeff, template = layers_processed_anchors[[35]]),
        raster(ratechange_AbsAnn_Pval, template = layers_processed_anchors[[35]]),
        raster(ratechange_AbsAnn_Rsquared, template = layers_processed_anchors[[35]]),
        raster(ratechange_AbsAnn_Coeff_Signif, template = layers_processed_anchors[[35]])),
        # raster(ratechange_AbsAnn_Rsquared_Signif, template = layers_processed_anchors[[35]])),
  filename = paste0('ratechange_AbsAnn.tif'), format="GTiff", overwrite = TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
# Barely able to execute this; nearly at max memory capacity; needed to not include significant layers into raster in order to complete.

#### load in ratechange_AbsAnn layers ####
ratechange_AbsAnn <- brick('ratechange_AbsAnn.tif')
names(ratechange_AbsAnn) <- c('Coefficients','P-value','R-squared')
# names(ratechange_AbsAnn) <- c('Coefficients','P-value','R-squared','Significant coefficients','Significant R-squared')
ratechange_AbsAnn

plot(ratechange_AbsAnn)
plot(ratechange_AbsAnn[[2]], breaks = c(0, 0.05,1), col = c('red','grey','white'))
plot(ratechange_AbsAnn[[3]], col = hsv(1, seq(0,1,length.out = 8) , 1))

# plot(AbsAnn.long$pixelvalue ~ AbsAnn.long$year)
hist(Coeffs)
