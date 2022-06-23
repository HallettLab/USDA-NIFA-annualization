# updated 2022-06-22
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

# set working directory to cleaned data folder
setwd("C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/MBF Experiment/Reed_USDA-NIFA_MBF Experiment/Data-cleaned")

# set path to R-output folder within Data-analysis
Rout <- "C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/MBF Experiment/Reed_USDA-NIFA_MBF Experiment/R-code/Data-analysis/R-output"

# read in cover data
cov.long <- read.csv('Cover_MBF_compiled_long.csv') # plots as rows
cov.wide <- read.csv('Cover_MBF_compiled.csv', header = FALSE) # plots as columns


# read in species info metadata to group species
fxgrps <- read.csv('species_info.csv')

# prep a cover dataframe with metadata; factor all the categorical variables
cover <- lapply(cov.long[,1:9], factor) %>%
  as.data.frame(.)

# prep functional group dataframe
fxgrps <- cbind(fxgrps, cov.wide[13:64,2:673])
colnames(fxgrps)[8:679] <- cov.wide[8,2:673]

# make cover measurements numeric
fxgrps[,8:679] <- sapply(fxgrps[,8:679],as.numeric)

# aggregate data by group2
agg.group2 <- aggregate(.~ group2, fxgrps[,c(4,8:679)], sum)
# transpose
agg.group2 <- setNames(data.frame(t(agg.group2[,-1])), agg.group2[,1])

# calculate native richness (still at month-level)
natives <- subset(fxgrps, group == 'Seeded native')
natives <- setNames(data.frame(t(natives[,-c(1:7)])), natives[,1])

# get the seeded forage species
forage <- subset(fxgrps, group2 == 'Forage grass' | group2 == 'Forage forb')
# aggregate by species2
forage <- aggregate(.~ species2, forage[,c(2,8:679)], sum)
forage <- setNames(data.frame(t(forage[,-1])), forage[,1])

#### add variables of interest to cover dataframe ####
cover <- cbind(cover, cov.long[,10:12]) # get bareground, litter, moss in there
cover$covTotal <- colSums(fxgrps[,8:679]) # total vegetative cover
cover <- cbind(cover, agg.group2) # summed cover for group2
cover <- cbind(cover, natives) # each individual native species
cover$natrich <- specnumber(natives) # native richness
cover <- cbind(cover, forage) # forage species

#### calculate average values across entire sampling period ####
cover.avg <- aggregate(.~ subplot.number, cover[,c(8,10:36)], mean)

#### total native richness (cumulative native species across time periods) ####
natives2 <- cbind(subplot = cov.long$subplot.number, natives[,-9])
agg.natives2 <- aggregate(.~ subplot, natives2, sum)
cover.avg$natrich <- specnumber(agg.natives2[,2:11])
cover.avg <- cbind(cover[cover$month=='2021-11',1:10], cover.avg[,-1])
cover.avg <- cover.avg[,-c(1:2)]
colnames(cover.avg)[8:9] <- c('Bareground_November','Bareground')
  

######## save processed cover data as .csv files ########
write.csv(cover, file = 'cover-processed2.csv', row.names = FALSE)
write.csv(cover.avg, file = 'cover_avg-processed2.csv', row.names = FALSE)
