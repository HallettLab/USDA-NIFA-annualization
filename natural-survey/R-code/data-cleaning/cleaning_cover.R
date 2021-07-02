# clean natural site survey cover datasets
# created: 2021-07-02 (script will be modified over project lifetime)

# script purpose:
# read in Pisgah and South Eugene Meadows plant cover files and combine into one dataframe
# build 2 types of cover dataframes:
## 1) long form (where site,block,stand,plot are columns, species are rows)
## 2) wide form (where species are columns, site,block,stand,plot are rows)
# write out clean cover datasets

# set working directory to access .csv files
setwd("C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/Natural Site Survey/Reed_USDA-NIFA_NaturalSurvey/Data-entered")

# set path to cleaned data folder
cl <- "C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/Natural Site Survey/Reed_USDA-NIFA_NaturalSurvey/Data-cleaned"

# read in both data
pisgah <- read.csv("cover_pisgah.csv",header=FALSE)
sem <- read.csv("cover_southeugenemeadows.csv",header=FALSE)

# the data are currently in long format (site/block/stand/plot as columns, species as rows).
# append sem data to pisgah data as additional columns. merge the dataframes based on similar species (i.e., column 1).
library(dplyr)
combined <- full_join(pisgah, sem, by = "V1")

# reorder the species alphabetically (leave off first 9 rows, which are not the species)
df <- combined[10:115,]
df <- df[order(df$V1),]
combined[10:115,] <- df; rm(df)

# replace NAs in species cover with zeros
combined[10:115,][is.na(combined[10:115,])] <- 0

# export long format to match it up with a functional group file I will create in Excel
colnames(combined) <- as.character(combined[1,]) # first replace the column headers with the first row
write.csv(combined[-1,],file=paste0(cl,"/cover_plots-columns.csv"),row.names = FALSE)

# transpose to wide format (plots as rows; species as columns) and export
combined.wide = setNames(data.frame(t(combined[,-1])), combined[,1])
rownames(combined.wide) = NULL
write.csv(combined.wide,file=paste0(cl,"/cover_plots-rows.csv"),row.names = FALSE)

