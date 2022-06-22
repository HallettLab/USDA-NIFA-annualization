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

#### calculate new raster that averages total cover across the 35 years ####
list_totherb <- sapply(layers_processed_anchors, "[[", 3) # extract the total herbaceous cover from each year in the list
totherb <- brick(list_totherb) # turn that list now into a brick

avg_totherb <- mean(totherb) # calculate average total herbaceous cover across years
plot(avg_totherb, col = 'black')
hist(avg_totherb)
writeRaster(avg_totherb, filename = 'avg_totalherbaceous.tif', format="GTiff", overwrite = TRUE,options="COMPRESS=LZW")


#### create mask to filter out any cells with <20% total herbaceous cover ####
mask_totherb20 <- mask(avg_totherb,
                       avg_totherb < 20, maskvalue = TRUE,
                  filename="mask_totherb20.tif", overwrite=TRUE)
mask_totherb20 <- raster('mask_totherb20.tif')
plot(mask_totherb20, col = 'black')

#### create new yearly RAP layers that are masked by minimum total herbaceous cover of 20% ####
for (i in 1:length(layers_processed_anchors)) {
  mask(layers_processed_anchors[[i]], mask_totherb20, maskvalue = NA,
       filename=paste0('RAP_processed_anchors_totherb20_',years[i],'.tif'), format=, overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
}

#### group data based on starting condition in 1986: <20% relative annual = perennial dominated; 20-50% relative annual = mixed; >50% relative annual = annual system ####
layer_1986_totherb20 <- mask(layers_processed_anchors[[1]], mask_totherb20, maskvalue = NA,
                             filename="mask_1986_totherb20.tif", format=, overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
layer_1986_totherb20
plot(layer_1986_totherb20[[4]], col='black')
hist(layer_1986_totherb20[[4]],breaks = 8)
hist(layers_processed_anchors[[1]][[4]],breaks = 8)

#### ratechange_AbsAnn ####
ratechange_AbsAnn <- brick('ratechange_AbsAnn.tif')
plot(ratechange_AbsAnn[[1]])

AnnDom <- mask(ratechange_AbsAnn,
               layer_1986_totherb20[[4]] < 0.5, maskvalue = TRUE)
# for some reason we need to redo the mask with the mask_totherb20 to make sure we're only including cells with >20% herbaceous cover.
AnnDom <- mask(AnnDom,mask_totherb20, maskvalue = NA,
               filename="ratechange_AbsAnn_processed_AnnDom.tif", format=, overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))

PerenDom <- mask(ratechange_AbsAnn,
                 layer_1986_totherb20[[4]] >= 0.2, maskvalue = TRUE)
PerenDom <- mask(PerenDom,mask_totherb20, maskvalue = NA,
                 filename="ratechange_AbsAnn_processed_PerenDom.tif", format=, overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))

Mixed <- mask(ratechange_AbsAnn,
              layer_1986_totherb20[[4]] < 0.2 | layer_1986_totherb20[[4]] >= 0.5, maskvalue = TRUE)
Mixed <- mask(Mixed,mask_totherb20, maskvalue = NA,
              filename="ratechange_AbsAnn_processed_Mixed.tif", format=, overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))

# finally, make a processed layer that includes cells with >20% herbaceous cover, but doesn't group by starting condition.
ratechange_AbsAnn_processed <- mask(ratechange_AbsAnn, mask_totherb20, maskvalue = NA,
                                    filename="ratechange_AbsAnn_processed.tif",
                                    format=, overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))