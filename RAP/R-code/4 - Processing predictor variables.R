library(raster)
require(rgdal)

setwd('C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/RAP analysis/R_working_directory')

# read in wvopc study area
wvopc <- readOGR(dsn = "C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/RAP analysis/ArcMap raster exports/WVOPC_Boundary.shp", layer="WVOPC_Boundary")
plot(wvopc)

# read in anchor sites
anchors <- readOGR(dsn = "C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/RAP analysis/ArcMap raster exports/WVOPC_Anchors_WGS_1984.shp", layer="WVOPC_Anchors_WGS_1984")
plot(anchors)

# set path to read in geotiff files of ssurgo and DEM predictor variables exported from ArcGIS
path <- 'C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/RAP analysis/ArcMap raster exports'

# load in rasters of each variable
ssurgo_pH <- raster(paste0(path,'/ssurgo_WV_pH.tif'))
ssurgo_aws100 <- raster(paste0(path,'/ssurgo_WV_aws100.tif'))
ssurgo_depth <- raster(paste0(path,'/ssurgo_WV_depth.tif'))
ssurgo_om <- raster(paste0(path,'/ssurgo_WV_om.tif'))
ssurgo_sand <- raster(paste0(path,'/ssurgo_WV_sand.tif'))
ssurgo_silt <- raster(paste0(path,'/ssurgo_WV_silt.tif'))
ssurgo_clay <- raster(paste0(path,'/ssurgo_WV_clay.tif'))
dem_slope <- raster(paste0(path, '/WV_slope.tif'))
dem_aspect <- raster(paste0(path, '/WV_aspect.tif'))
dem_heatload <- raster(paste0(path, '/WV_heatload.tif')) # this is heatload calculated in ArcMap 10.5 using the Geomorphometry and Gradients Metrics toolbox: https://evansmurphy.wixsite.com/evansspatial/arcgis-gradient-metrics-toolbox
dem_elev <- raster('C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/RAP analysis/rasters_SRTMGL1.tar/output_SRTMGL1.tif')

#### throw all these rasters into a list ####
predictors <- list(ssurgo_aws100,ssurgo_depth,ssurgo_pH,ssurgo_om,ssurgo_sand,ssurgo_silt,ssurgo_clay,
                   dem_elev,dem_slope,dem_aspect,dem_heatload)


#### process predictor variables
# use the ratechange_absann_processed raster as the layer for resampling, clipping extent, and masking.
mask <- raster('ratechange_AbsAnn_processed.tif')
plot(mask)

# first you need to resample predictor layers to change resolution to 30x30
# then you need to mask each layer by the ratechange raster

# use a loop to spit out resampled and masked ssurgo files
vars <- c('aws100','depth','pH','om','sand','silt','clay','elev','slope','aspect','heatload')
for (i in 1:7) {
  c <- resample(predictors[[i]], mask, method='bilinear', filename = paste0('ssurgo_processed_',vars[i],'.tif'), overwrite = TRUE)
  m <- mask(c, mask, filename = paste0('ssurgo_processed_',vars[i],'.tif'), overwrite = TRUE)
}

# use a loop to spit out resampled and masked DEM files
for (i in 8:11) {
  c <- resample(predictors[[i]], mask, method='bilinear', filename = paste0('dem_processed_',vars[i],'.tif'), overwrite = TRUE)
  m <- mask(c, mask, filename = paste0('dem_processed_',vars[i],'.tif'), overwrite = TRUE)
}
