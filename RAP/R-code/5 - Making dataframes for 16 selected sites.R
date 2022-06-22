library(raster)
library(rgdal)
library(ggplot2)
library(tidyr)
library(dplyr)
library(car)
library(lme4)
library(sjPlot)
library(MuMIn)
library(emmeans)
library(quantreg)
library(lqmm)

# set working directory
setwd("C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/RAP analysis/R_working_directory")

# set directory for figure exports
Rout <- 'C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/RAP analysis/R_working_directory/R figure exports'

# prep for loading each yearly file at once into a list
years <-1986:2020


######## load in year-by-year RAP layers to a list ########
layers_processed_anchors <- list()
for (i in 1:length(years)) {
  filename <- paste0('RAP_processed_anchors_totherb20_',years[i],'.tif')
  layer <- brick(filename)
  names(layer) <- c(paste0('Annual_',years[i]), paste0('Perennial_',years[i]), paste0('Total_',years[i]),
                    paste0('Rel_Ann_',years[i]), paste0('Rel_Peren_',years[i]), paste0('LogRatio_',years[i]))
  layers_processed_anchors <- c(layers_processed_anchors, list(layer))
  rm(filename); rm(layer) # remove temporary objects to save memory
}
plot(layers_processed_anchors[[1]][[1]])

######## load in landscape driver predictor variables to a list ########
predictors <- c('aws100','depth','pH','om','sand','silt','clay','elev','slope','aspect','heatload')
type <- c(rep('ssurgo',7),rep('dem',4))

xvars <- list()
for (i in 1:length(predictors)) {
  filename <- paste0(type[i],'_processed_',predictors[i],'.tif')
  layer <- raster(filename)
  xvars <- c(xvars, list(layer))
  rm(filename); rm(layer) # remove temporary objects to save memory
}
plot(xvars[[1]])


###### Howard Buford Recreation Area (hbra) ######

#### load in HBRA shapefile
hbra <- readOGR(dsn = "C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/RAP analysis/RAPanalysis-anchors/HBRA/HBRA.shp", layer="HBRA")
plot(hbra)

#### clip annual and perennial layers to hbra
annual_hbra <- list()
perennial_hbra <- list()
for (i in 1:length(layers_processed_anchors)) {
  c1 <- crop(layers_processed_anchors[[i]][[1]], extent(hbra))
  m1 <- mask(c1, hbra)
  annual_hbra <- c(annual_hbra, list(m1))
  
  c2 <- crop(layers_processed_anchors[[i]][[2]], extent(hbra))
  m2 <- mask(c2, hbra)
  perennial_hbra <- c(perennial_hbra, list(m2))
  rm(c1);rm(m1);rm(c2);rm(m2)
}
plot(annual_hbra[[1]], main = '1986 Absolute annual cover'); plot(hbra, bg="transparent", add=TRUE)
# plot(perennial_hbra[[1]], main = '1986 Absolute perennial cover'); plot(hbra, bg="transparent", add=TRUE)

#### convert lists into a brick for easy conversion to an analysis- and plot-friendly dataframe
annual_hbra <- brick(annual_hbra)
perennial_hbra <- brick(perennial_hbra)

#### convert bricks to dataframe
annual_hbra <- as.data.frame(annual_hbra, xy = TRUE)
perennial_hbra <- as.data.frame(perennial_hbra, xy = TRUE)

#### convert wide to long
annual_hbra <- gather(annual_hbra, year, absann, Annual_1986:Annual_2020, factor_key=TRUE, na.rm = T)
perennial_hbra <- gather(perennial_hbra, year, absperen, Perennial_1986:Perennial_2020, factor_key=TRUE, na.rm = T)

#### turn into 1 df; adjust year column; add pixel ID column to beginning
df_hbra <- annual_hbra; df_hbra$absperen <- perennial_hbra$absperen
df_hbra$year <- as.numeric(gsub("Annual_", "", df_hbra$year))
df_hbra <- cbind(pixel = factor(paste(df_hbra$x, df_hbra$y, sep = ',')),
                 df_hbra)
df_hbra <- cbind(site = 'Howard Buford Recreation Area',df_hbra)

#### calculate additional variables (total, relative covers, logratio)
df_hbra$total <- df_hbra$absann + df_hbra$absperen
df_hbra$relann <- df_hbra$absann/df_hbra$total
df_hbra$relperen <- df_hbra$absperen/df_hbra$total
df_hbra$logratio <- log((df_hbra$absann+1) / (df_hbra$absperen+1))

#### add in variable for relann_1986 so we can categorize the data by starting condition
hbra_1986 <- subset(df_hbra, year == '1986')
hbra_1986 <- hbra_1986[,-c(2:7,9:11)]; colnames(hbra_1986)[2] <- 'relann_1986'
df_hbra <- merge(df_hbra, hbra_1986, by= 'pixel')
df_hbra$condition <- "Mixed"
df_hbra[df_hbra$relann_1986 >= 0.5,]$condition <- 'Annual Dom.'
df_hbra[df_hbra$relann_1986 < 0.2,]$condition <- 'Peren. Dom.'

#### clip xvars to HBRA area
xvars_hbra <- lapply(xvars, function(x) crop(x, extent(hbra)))
xvars_hbra <- lapply(xvars_hbra, function(x) mask(x, hbra))

#### convert list into a brick
xvars_hbra <- brick(xvars_hbra)
names(xvars_hbra) <- predictors
xvars_df <- as.data.frame(xvars_hbra, xy = TRUE)

#### prep xvars_df to merge with df_hbra
xvars_df$pixel <- factor(paste(xvars_df$x, xvars_df$y, sep = ','))
xvars_df <- xvars_df[,-c(1:2)]

#### calculate northness
deg2rad <- function(deg) {(deg * pi) / (180)}
xvars_df$northness <- deg2rad(xvars_df$aspect)
xvars_df$northness <- cos(xvars_df$northness)

#### fix aspect (-1 value)
xvars_df[!is.na(xvars_df$aspect) & xvars_df$aspect < 0,]$aspect <- 360 + xvars_df[!is.na(xvars_df$aspect) & xvars_df$aspect < 0,]$aspect

#### assign quantiles to predictor variables
# library(gtools)
xvars_df$heatload_quants <- gtools::quantcut(xvars_df$heatload, q = 3)
# xvars_df$aspect_quants <- gtools::quantcut(xvars_df$aspect, q = 5)
# xvars_df$northness_quants <- gtools::quantcut(xvars_df$northness, q = 5)
xvars_df$aspect_cats <- NA
xvars_df[!is.na(xvars_df$aspect),]$aspect_cats <- 'North'
xvars_df[!is.na(xvars_df$aspect) & xvars_df$aspect > 45,]$aspect_cats <- 'East'
xvars_df[!is.na(xvars_df$aspect) & xvars_df$aspect > 135,]$aspect_cats <- 'South'
xvars_df[!is.na(xvars_df$aspect) & xvars_df$aspect > 225,]$aspect_cats <- 'West'
xvars_df[!is.na(xvars_df$aspect) & xvars_df$aspect > 315,]$aspect_cats <- 'North'

xvars_df$slope_quants <- gtools::quantcut(xvars_df$slope, q = 3)
xvars_df$elev_quants <- gtools::quantcut(xvars_df$elev, q = 3)
xvars_df$depth_quants <- gtools::quantcut(xvars_df$depth, q = 2)
quantile(xvars_df$depth, na.rm=T, probs = seq(0,1,(1/3)))
xvars_df$aws100_quants <- gtools::quantcut(xvars_df$aws100, q = 2)
quantile(xvars_df$aws100, na.rm=T, probs = seq(0,1,(1/2)))
xvars_df$aws100_quants <- gtools::quantcut(xvars_df$aws100, q = 2)
quantile(xvars_df$aws100, na.rm=T, probs = seq(0,1,(1/2)))
xvars_df$sand_quants <- gtools::quantcut(xvars_df$sand, q = 3)
quantile(xvars_df$sand, na.rm=T, probs = seq(0,1,(1/3)))
xvars_df$silt_quants <- gtools::quantcut(xvars_df$silt, q = 3)
quantile(xvars_df$silt, na.rm=T, probs = seq(0,1,(1/4)))
xvars_df$clay_quants <- gtools::quantcut(xvars_df$clay, q = 3)
quantile(xvars_df$clay, na.rm=T, probs = seq(0,1,(1/3)))

#### merge to df_hbra by x,y
df_hbra <- merge(df_hbra, xvars_df, by = 'pixel')

#### create column for 1/2 decade
df_hbra$years <- "'86-90"
df_hbra[df_hbra$year > 1990,]$years <- "'91-95"
df_hbra[df_hbra$year > 1995,]$years <- "'96-00"
df_hbra[df_hbra$year > 2000,]$years <- "'01-05"
df_hbra[df_hbra$year > 2005,]$years <- "'06-10"
df_hbra[df_hbra$year > 2010,]$years <- "'11-15"
df_hbra[df_hbra$year > 2015,]$years <- "'16-20"
df_hbra$years <- factor(df_hbra$years, levels = c("'86-90","'91-95","'96-00","'01-05","'06-10","'11-15","'16-20"))

#### sort by year, pixel ID
df_hbra <- df_hbra[order(df_hbra$pixel, df_hbra$year),]

#### save dataframe to a file
# write.csv(df_hbra, file = 'df_hbra.csv', row.names = FALSE)
saveRDS(df_hbra, file = "df_hbra.rds")
# df_hbra <- readRDS(file = "df_hbra.rds")


###### All other candidate anchor sites ######
#### load in anchors shapefile
anchors <- readOGR(dsn = "C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/RAP analysis/ArcMap raster exports/WVOPC_Anchors_WGS_1984.shp", layer="WVOPC_Anchors_WGS_1984")
plot(anchors)

anchors$sitename
sites <- c('Andrew Reasoner Wildlife Preserve','Coburg Ridge','Murray Hill, Upper Willow Creek','South Eugene Meadows',
           'Quamash Prairie','Thurston Hills','Suzanne Arlie, Mount Baldy','Finley National Wildlife Refuge','Fitton Green County Natural Area',
           'Bald Hill Natural Area','Chip Ross Park','Kingston Prairie Preserve','Ankeny National Wildlife Refuge','Baskett Slough National Wildlife Refuge',
           'Cerro Gordo','Bald Hill Farm, Mulkey Ridge','Rattlesnake Butte','Native Oaks Ridge')
sort(sites)
#### clip annual and perennial layers
polygs_annual <- list()
polygs_perennial <- list()

for (j in 1:length(sites)) {
  site <- subset(anchors, sitename == sites[j])
  annual <- list()
  perennial <- list()
  
    for (i in 1:length(layers_processed_anchors)) {
      c1 <- crop(layers_processed_anchors[[i]][[1]], extent(site))
      m1 <- mask(c1, site)
      annual <- c(annual, list(m1))
    
      c2 <- crop(layers_processed_anchors[[i]][[2]], extent(site))
      m2 <- mask(c2, site)
      perennial <- c(perennial, list(m2))
      
      rm(c1);rm(m1);rm(c2);rm(m2)
    }
  polygs_annual <- c(polygs_annual, list(annual))
  polygs_perennial <- c(polygs_perennial, list(perennial))
}

# plot(annual_hbra[[1]], main = '1986 Absolute annual cover'); plot(hbra, bg="transparent", add=TRUE)
# plot(perennial_hbra[[1]], main = '1986 Absolute perennial cover'); plot(hbra, bg="transparent", add=TRUE)

#### convert lists into a brick for easy conversion to an analysis- and plot-friendly dataframe
polygs_annual2 <- lapply(polygs_annual, function(x) brick(x))
polygs_perennial2 <- lapply(polygs_perennial, function(x) brick(x))

#### convert bricks to dataframe
polygs_annual2 <- lapply(polygs_annual2, function(x) as.data.frame(x, xy = TRUE)) 
polygs_perennial2 <- lapply(polygs_perennial2, function(x) as.data.frame(x, xy = TRUE)) 

#### convert wide to long
polygs_annual2 <- lapply(polygs_annual2, function(x) gather(x, year, absann, Annual_1986:Annual_2020, factor_key=TRUE, na.rm = T))
polygs_perennial2 <- lapply(polygs_perennial2, function(x) gather(x, year, absperen, Perennial_1986:Perennial_2020, factor_key=TRUE, na.rm = T))

#### add site name into each dataframe
polygs_annual2 <- mapply(`[<-`, polygs_annual2, 'site', value = sites, SIMPLIFY = FALSE)
polygs_perennial2 <- mapply(`[<-`, polygs_perennial2, 'site', value = sites, SIMPLIFY = FALSE)

#### turn list of site dataframes into one long df
polygs_annual2 <- do.call("rbind", polygs_annual2)
polygs_perennial2 <- do.call("rbind", polygs_perennial2)

#### turn into 1 df; adjust year and site columns; add pixel ID column to beginning
df_polygs <- polygs_annual2; df_polygs$absperen <- polygs_perennial2$absperen
df_polygs$year <- as.numeric(gsub("Annual_", "", df_polygs$year))
df_polygs <- cbind(pixel = factor(paste(df_polygs$x, df_polygs$y, sep = ',')),
                   df_polygs)
df_polygs <- cbind(site = df_polygs$site,df_polygs[,-6])

#### rbind with df_hbra
df_polygs <- rbind(df_hbra, df_polygs)

#### calculate additional variables (total, relative covers, logratio)
df_polygs$total <- df_polygs$absann + df_polygs$absperen
df_polygs$relann <- df_polygs$absann/df_polygs$total
df_polygs$relperen <- df_polygs$absperen/df_polygs$total
df_polygs$logratio <- log((df_polygs$absann+1) / (df_polygs$absperen+1))

#### add in variable for relann_1986 so we can categorize the data by starting condition
polygs_1986 <- subset(df_polygs, year == '1986')
polygs_1986 <- polygs_1986[,c(2,9)]; colnames(polygs_1986)[2] <- 'relann_1986'
df_polygs <- merge(df_polygs, polygs_1986, by= 'pixel')
df_polygs$condition <- "Mixed"
df_polygs[!is.na(df_polygs$relann_1986) & df_polygs$relann_1986 >= 0.5,]$condition <- 'Annual Dom.'
df_polygs[!is.na(df_polygs$relann_1986) & df_polygs$relann_1986 < 0.2,]$condition <- 'Peren. Dom.'
df_polygs[is.nan(df_polygs$relann_1986),]$condition <- NA





#### clip xvars to polygs area
sites <- c('Howard Buford Recreation Area',sites)
polygs <-   anchors[anchors$sitename %in% sites,]
plot(anchors)
plot(polygs,col = 'red', bg="transparent", add=TRUE)
# xvars_polygs <- lapply(xvars, function(x) crop(x, extent(polygs)))
# xvars_polygs <- lapply(xvars_polygs, function(x) mask(x, polygs))

xvars_polygs <- list()

for (j in 1:length(sites)) {
  site <- subset(anchors, sitename == sites[j])
  xvars_poly <- list()

  for (i in 1:length(xvars)) {
    c1 <- crop(xvars[[i]], extent(site))
    m1 <- mask(c1, site)
    xvars_poly <- c(xvars_poly, list(m1))
    
    rm(c1);rm(m1)
  }
  xvars_polygs2 <- c(xvars_polygs, list(xvars_poly))
}


#### convert list into a brick
xvars_polygs <- brick(xvars_polygs)
xvars_polygs <- lapply(xvars_polygs, function(x) brick(x))

#### convert bricks to dataframe
xvars_polygs <- lapply(xvars_polygs, function(x) as.data.frame(x, xy = TRUE)) 

#### add site name into each dataframe
xvars_polygs <- mapply(`[<-`, xvars_polygs, 'site', value = sites, SIMPLIFY = FALSE)

#### turn list of site dataframes into one long df
xvars_polygs <- do.call("rbind", xvars_polygs)
colnames(xvars_polygs)[3:13] <- predictors

#### prep xvars_df to merge with df_polygs
xvars_polygs$pixel <- factor(paste(xvars_polygs$x, xvars_polygs$y, sep = ','))
xvars_polygs <- xvars_polygs[,-c(1:2,14)]

#### calculate northness
deg2rad <- function(deg) {(deg * pi) / (180)}
xvars_polygs$northness <- deg2rad(xvars_polygs$aspect)
xvars_polygs$northness <- cos(xvars_polygs$northness)

#### fix aspect (-1 value)
xvars_polygs[!is.na(xvars_polygs$aspect) & xvars_polygs$aspect < 0,]$aspect <- 360 + xvars_polygs[!is.na(xvars_polygs$aspect) & xvars_polygs$aspect < 0,]$aspect

#### assign quantiles to predictor variables (will need to do this separately for each site???)
# # library(gtools)
# xvars_df$heatload_quants <- gtools::quantcut(xvars_df$heatload, q = 3)
# # xvars_df$aspect_quants <- gtools::quantcut(xvars_df$aspect, q = 5)
# # xvars_df$northness_quants <- gtools::quantcut(xvars_df$northness, q = 5)
# xvars_df$aspect_cats <- NA
# xvars_df[!is.na(xvars_df$aspect),]$aspect_cats <- 'North'
# xvars_df[!is.na(xvars_df$aspect) & xvars_df$aspect > 45,]$aspect_cats <- 'East'
# xvars_df[!is.na(xvars_df$aspect) & xvars_df$aspect > 135,]$aspect_cats <- 'South'
# xvars_df[!is.na(xvars_df$aspect) & xvars_df$aspect > 225,]$aspect_cats <- 'West'
# xvars_df[!is.na(xvars_df$aspect) & xvars_df$aspect > 315,]$aspect_cats <- 'North'
# 
# xvars_df$slope_quants <- gtools::quantcut(xvars_df$slope, q = 3)
# xvars_df$elev_quants <- gtools::quantcut(xvars_df$elev, q = 3)
# xvars_df$depth_quants <- gtools::quantcut(xvars_df$depth, q = 2)
# quantile(xvars_df$depth, na.rm=T, probs = seq(0,1,(1/3)))
# xvars_df$aws100_quants <- gtools::quantcut(xvars_df$aws100, q = 2)
# quantile(xvars_df$aws100, na.rm=T, probs = seq(0,1,(1/2)))
# xvars_df$aws100_quants <- gtools::quantcut(xvars_df$aws100, q = 2)
# quantile(xvars_df$aws100, na.rm=T, probs = seq(0,1,(1/2)))
# xvars_df$sand_quants <- gtools::quantcut(xvars_df$sand, q = 3)
# quantile(xvars_df$sand, na.rm=T, probs = seq(0,1,(1/3)))
# xvars_df$silt_quants <- gtools::quantcut(xvars_df$silt, q = 3)
# quantile(xvars_df$silt, na.rm=T, probs = seq(0,1,(1/4)))
# xvars_df$clay_quants <- gtools::quantcut(xvars_df$clay, q = 3)
# quantile(xvars_df$clay, na.rm=T, probs = seq(0,1,(1/3)))

#### merge to df_polygs by x,y
df_polygons <- merge(df_polygs, xvars_polygs, by = 'pixel')
# after merging, there SHOULD be the same number of rows; it should just be pasting in the xvars data to the df_polygs data when pixelIDs match
# it is actually increasing the number of rows somewhat substantially (by 61,005), which divided by 35 = 1743 pixels. 
# my best guess is that there are some areas within different sites (e.g., border of Bald Hill and Mulkey) that have the same pixel ID.

#### create column for 1/2 decade
df_polygons$years <- "'86-90"
df_polygons[df_polygons$year > 1990,]$years <- "'91-95"
df_polygons[df_polygons$year > 1995,]$years <- "'96-00"
df_polygons[df_polygons$year > 2000,]$years <- "'01-05"
df_polygons[df_polygons$year > 2005,]$years <- "'06-10"
df_polygons[df_polygons$year > 2010,]$years <- "'11-15"
df_polygons[df_polygons$year > 2015,]$years <- "'16-20"
df_polygons$years <- factor(df_polygons$years, levels = c("'86-90","'91-95","'96-00","'01-05","'06-10","'11-15","'16-20"))

#### sort by year, pixel ID
df_polygons <- df_polygons[order(df_polygons$pixel, df_polygons$year),]

#### collapse Bald Hill, Fitton Green sites to a single complex (all completely adjacent)
levels(as.factor(df_polygons$site))
df_polygons[df_polygons$site == 'Bald Hill Farm, Mulkey Ridge' |
              df_polygons$site == 'Bald Hill Natural Area' |
              df_polygons$site == 'Fitton Green County Natural Area',]$site <- 'Bald Hill-Fitton Green Complex' 
  
# throw away Quamash Prairie, which has a lot of pixels within ponds and highly disturbed ground
df_polygons <- subset(df_polygons, !site == 'Quamash Prairie')


#### save dataframe to a file
# write.csv(df_polygons, file = 'df_polygons.csv', row.names = FALSE)
saveRDS(df_polygons, file = "df_polygons.rds")
# df_polygons <- readRDS(file = "df_polygons.rds")

#### load df
df_polygons <- readRDS(file = "df_polygons.rds")
