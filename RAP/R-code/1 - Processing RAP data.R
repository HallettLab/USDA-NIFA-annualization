library(raster)
require(rgdal)

# set wd
setwd('C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/RAP analysis/R_working_directory')

# set path to read in geotiff files downloaded from RAP
RAPpath <- 'C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/RAP analysis/RAP data downloads/drive-download-20220205T005002Z-001/'

# read in wvopc study area
wvopc <- readOGR(dsn = "C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/RAP analysis/ArcMap raster exports/WVOPC_Boundary.shp", layer="WVOPC_Boundary")
plot(wvopc)

# read in candidate anchor sites
anchors <- readOGR(dsn = "C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/RAP analysis/ArcMap raster exports/WVOPC_Anchors_WGS_1984.shp", layer="WVOPC_Anchors_WGS_1984")
plot(anchors)

# prep for loading each yearly file at once into a list
years <- 1986:2020

#### load in RAP datasets as bricks in a list ####
layers_year <- list() # each element in this list will be a yearly RAP tif file
for (i in 1:length(years)) {
  filename <- paste0(RAPpath,'RAP_VegCover_',years[i],'.tif')
  layer <- brick(filename)
  layer <- layer[[-c(2:3)]] # remove the bare ground and litter bands; now, within each brick, band 1 = annual, band 2 = perennial
  names(layer) <- c(paste0('AFG_',years[i]),paste0('PFG_',years[i]))
  layers_year <- c(layers_year, list(layer))
  rm(filename); rm(layer) # remove temporary objects to save memory
}

# plot one of these to look at it
layers_year[[1]]
plot(layers_year[[1]])


#### clip NLCD raster to WVOPC study area; save output as tif file ####
library(raster)
nlcd <- raster("C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/RAP analysis/ArcMap raster exports/NLCD_2019_Land_Cover_L48_20211_wgs1984.tif")
nlcd
plot(nlcd)

nlcd.crop <- crop(nlcd, extent(wvopc))
nlcd.crop2 <- mask(nlcd.crop, wvopc, filename="NLCD_cropped_wgs1984.tif", overwrite=TRUE)
plot(nlcd.crop)
plot(nlcd.crop2)

#### mask nlcd to exclude certain cover classes ####
# filter out the following cover classes:
# 11(open water), 12(perennial ice/snow), 21(developed,open), 22(developed,low), 23(developed,medium), 24(developed,high),
# 31(Barren land), 41(deciduous forest), 42(evergreen), 43(mixed forest), 82(cultivated crops), 90(woody wetlands), and 95(emergent herbaceous wetlands).
# This leaves us with just 52(shrub/scrub), 71(grassland/herbaceous), and 81(pasture/hay). Decided against keeping 95(emergent herbaceous wetlands) after
# carefully analyzing satellite imagery where this cover class occurs and deciding that it is a poor representative of grassland (there is a lot of tree cover,
# bare ground, true marshy swamp, actual open water, etc.)
nlcd_mask <- mask(nlcd.crop2,
                  nlcd.crop2 == 11 | nlcd.crop2 == 12 | nlcd.crop2 == 21 | nlcd.crop2 == 22 |
                    nlcd.crop2 == 23 | nlcd.crop2 == 24 | nlcd.crop2 == 31 | nlcd.crop2 == 41 |
                    nlcd.crop2 == 42 | nlcd.crop2 == 43 | nlcd.crop2 == 82 | nlcd.crop2 == 90 | nlcd.crop2 == 95, maskvalue = TRUE,
                  filename="nlcd_mask.tif", overwrite=TRUE)
plot(nlcd_mask)


#### clip RAP datasets to nlcd_mask ####
nlcd_mask <- raster('nlcd_mask.tif')
plot(nlcd_mask)

# first you need to crop one of the RAP layers to wvopc extent
c <- crop(layers_year[[i]], extent(wvopc))
plot(c)
# need to resample NLCD data to change resolution to that of layers_year
nlcd_mask_resample <- resample(nlcd_mask, c, method='bilinear', filename="nlcd_mask_resample.tif", overwrite=TRUE)
plot(nlcd_mask_resample)

# now use a loop to spit out cropped and masked layers_year files
# layers_processed <- list()
for (i in 1:length(layers_year)) {
  c <- crop(layers_year[[i]], extent(wvopc), filename = paste0('RAP_processed_',years[i],'.tif'), overwrite = TRUE)
  m <- mask(c, nlcd_mask_resample, filename = paste0('RAP_processed_',years[i],'.tif'), overwrite = TRUE)
  # layers_processed <- c(layers_processed, list(m))
}


#### calculate 4 new layers for each year: total herbaceous cover, relative annual, relative perennial, and log ratio ####
# stick those new layers into the layers_processed raster bricks to create 6-layer files for each year.

# read in layers_processed
layers_processed <- list()
for (i in 1:length(years)) {
  filename <- paste0('RAP_processed_',years[i],'.tif')
  layer <- brick(filename)
  names(layer) <- c(paste0('AFG_',years[i]),paste0('PFG_',years[i]))
  layers_processed <- c(layers_processed, list(layer))
  rm(filename); rm(layer) # remove temporary objects to save memory
}


# write new layers_processed files with 6 bands
for (i in 1:length(layers_processed)) {
  writeRaster(
    stack(layers_processed[[i]],
          overlay(layers_processed[[i]], fun=function(x,y) { x + y }),
          overlay(layers_processed[[i]], fun=function(x,y) { x/(x + y) }),
          overlay(layers_processed[[i]], fun=function(x,y) { y/(x + y) }),
          overlay(layers_processed[[i]], fun=function(x,y) { log((x+1)/(y+1)) })),
    filename = paste0('RAP_processed2_',years[i],'.tif'), overwrite = TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
}



#### load in processed RAP layers ####
layers_processed2 <- list()
for (i in 1:length(years)) {
  filename <- paste0('RAP_processed2_',years[i],'.tif')
  layer <- brick(filename)
  names(layer) <- c(paste0('Annual_',years[i]), paste0('Perennial_',years[i]), paste0('Total_',years[i]),
                    paste0('Rel_Ann_',years[i]), paste0('Rel_Peren_',years[i]), paste0('LogRatio_',years[i]))
  layers_processed2 <- c(layers_processed2, list(layer))
  rm(filename); rm(layer) # remove temporary objects to save memory
}

layers_processed2[[1]]
plot(layers_processed2[[1]])




#### clip processed RAP layers to anchor sites ####
for (i in 1:length(layers_processed2)) {
  c <- crop(layers_processed2[[i]], extent(anchors))
  m <- mask(c, anchors, filename = paste0('RAP_processed_anchors_',years[i],'.tif'), overwrite = TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
}

#### load in processed anchor layers ####
layers_processed_anchors <- list()
for (i in 1:length(years)) {
  filename <- paste0('RAP_processed_anchors_',years[i],'.tif')
  layer <- brick(filename)
  names(layer) <- c(paste0('Annual_',years[i]), paste0('Perennial_',years[i]), paste0('Total_',years[i]),
                    paste0('Rel_Ann_',years[i]), paste0('Rel_Peren_',years[i]), paste0('LogRatio_',years[i]))
  layers_processed_anchors <- c(layers_processed_anchors, list(layer))
  rm(filename); rm(layer) # remove temporary objects to save memory
}

layers_processed_anchors[[1]]
plot(layers_processed_anchors[[1]])
