# updated 2022-04-08
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
setwd("C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/Natural Site Survey/Reed_USDA-NIFA_NaturalSurvey/Data-cleaned")

# set path to R-output folder within Data-analysis
Rout <- "C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/Natural Site Survey/Reed_USDA-NIFA_NaturalSurvey/R-code/Data-analysis/R-output"

# read in cover data
cov.long <- read.csv('cover_plots-rows.csv') # plots as rows
cov.wide <- read.csv('cover_plots-columns.csv', header = FALSE) # plots as columns

# subset data to subsample plots
cov.long <- subset(cov.long, !is.na(transect))
cov.wide <- cov.wide[,-c(2,6,10,14,18,22,26,30,34,38,
                        42,46,50,54,58,62,66,70,74,78)]

# subsample plots were 1m2. there were 3 plots per stand (meaning 6 plots per block; 60 total across sites).
# within plots, I estimated cover to species level, estimated litter and bareground, and measured litter depth.

# create a matrix of only the species as columns
mtrx <- cov.long[,11:116]

# read in functional group metadata to group species
fxgrps <- read.csv('functional_groups.csv')

# prep a cover dataframe with metadata; factor all the categorical variables
cover <- lapply(cov.long[,1:6], factor) %>%
  as.data.frame(.)

# prep functional group dataframe
fxgrps <- cbind(fxgrps, cov.wide[11:116,2:61])
colnames(fxgrps)[12:71] <- cov.wide[3,2:61]

# make cover measurements numeric
fxgrps[,12:71] <- sapply(fxgrps[,12:71],as.numeric)

#### calculate total and functional group covers ####
cover <- cbind(cover, cov.long[,7:10])
cover$covTotal <- colSums(fxgrps[,12:31])
cover$covNative <- colSums(fxgrps[fxgrps$nativity=="Native",12:71])
cover$covIntro <- colSums(fxgrps[fxgrps$nativity=="Introduced",12:71])
cover$covAnnual <- colSums(fxgrps[fxgrps$duration=="Annual",12:71])
cover$covPeren <- colSums(fxgrps[fxgrps$duration=="Perennial",12:71])
cover$covGram <- colSums(fxgrps[fxgrps$form2=="Graminoid",12:71])
cover$covForb <- colSums(fxgrps[fxgrps$form2=="Forb",12:71])
cover$covNatAnn <- colSums(fxgrps[fxgrps$nativity=="Native"&fxgrps$duration=="Annual",12:71])
cover$covNatPer <- colSums(fxgrps[fxgrps$nativity=="Native"&fxgrps$duration=="Perennial",12:71])
cover$covIntAnn <- colSums(fxgrps[fxgrps$nativity=="Introduced"&fxgrps$duration=="Annual",12:71])
cover$covIntPer <- colSums(fxgrps[fxgrps$nativity=="Introduced"&fxgrps$duration=="Perennial",12:71])
cover$covAG <- colSums(fxgrps[fxgrps$form2=="Graminoid"&fxgrps$duration=="Annual",12:71])
cover$covPG <- colSums(fxgrps[fxgrps$form2=="Graminoid"&fxgrps$duration=="Perennial",12:71])
cover$covAF <- colSums(fxgrps[fxgrps$form2=="Forb"&fxgrps$duration=="Annual",12:71])
cover$covPF <- colSums(fxgrps[fxgrps$form2=="Forb"&fxgrps$duration=="Perennial",12:71])
cover$covNAG <- colSums(fxgrps[fxgrps$funcgroup=="NAG",12:71])
cover$covNAF <- colSums(fxgrps[fxgrps$funcgroup=="NAF",12:71])
cover$covNPG <- colSums(fxgrps[fxgrps$funcgroup=="NPG",12:71])
cover$covNPF <- colSums(fxgrps[fxgrps$funcgroup=="NPF",12:71])
cover$covIAG <- colSums(fxgrps[fxgrps$funcgroup=="IAG",12:71])
cover$covIAF <- colSums(fxgrps[fxgrps$funcgroup=="IAF",12:71])
cover$covIPG <- colSums(fxgrps[fxgrps$funcgroup=="IPG",12:71])
cover$covIPF <- colSums(fxgrps[fxgrps$funcgroup=="IPF",12:71])
cover$covWoody <- colSums(fxgrps[fxgrps$funcgroup=="Woody",12:71])
cover$covUnk <- colSums(fxgrps[fxgrps$funcgroup=="Unk",12:71])

#### calculate relative covers ####
relcover<-cover[,1:6]
relcover$relNative<-cover$covNative/cover$covTotal
relcover$relIntro<-cover$covIntro/cover$covTotal
relcover$relAnnual<-cover$covAnnual/cover$covTotal
relcover$relPeren<-cover$covPeren/cover$covTotal
relcover$relGram<-cover$covGram/cover$covTotal
relcover$relForb<-cover$covForb/cover$covTotal
relcover$relNatAnn<-cover$covNatAnn/cover$covTotal
relcover$relNatPer<-cover$covNatPer/cover$covTotal
relcover$relIntAnn<-cover$covIntAnn/cover$covTotal
relcover$relIntPer<-cover$covIntPer/cover$covTotal
relcover$relNAG<-cover$covNAG/cover$covTotal
relcover$relNAF<-cover$covNAF/cover$covTotal
relcover$relNPG<-cover$covNPG/cover$covTotal
relcover$relNPF<-cover$covNPF/cover$covTotal
relcover$relIAG<-cover$covIAG/cover$covTotal
relcover$relIAF<-cover$covIAF/cover$covTotal
relcover$relIPG<-cover$covIPG/cover$covTotal
relcover$relIPF<-cover$covIPF/cover$covTotal
relcover$relWoody<-cover$covWoody/cover$covTotal
relcover$relUnk<-cover$covUnk/cover$covTotal



#### calculate various richness measurements, diversity, etc. ####
cover$rich <- specnumber(mtrx)

# native richness
mtrx.nat <- subset(fxgrps, nativity =='Native') %>%
  t(.) %>%
  .[12:71,] %>%
  as.data.frame(.) %>%
  sapply(.,as.numeric) %>%
  as.data.frame(.)
cover$rich.nat <- specnumber(mtrx.nat)

# forb richness
mtrx.forb <- subset(fxgrps, form =='Forb') %>%
  t(.) %>%
  .[12:71,] %>%
  as.data.frame(.) %>%
  sapply(.,as.numeric) %>%
  as.data.frame(.)
cover$rich.forb <- specnumber(mtrx.forb)

# native forb richness
mtrx.natforb <- subset(fxgrps, nativity == 'Native' & form =='Forb') %>%
  t(.) %>%
  .[12:71,] %>%
  as.data.frame(.) %>%
  sapply(.,as.numeric) %>%
  as.data.frame(.)
cover$rich.natforb <- specnumber(mtrx.natforb)

cover$diversity <- diversity(mtrx, index = "simpson")
cover$dominance <- 1 - cover$diversity # Simpson's dominance
cover$invsimp <- 1/cover$dominance
cover$evenness <- cover$invsimp/cover$rich


######## save processed cover data as .csv files ########
write.csv(cover, file = 'community_cover_plots-processed.csv')
write.csv(relcover, file = 'community_relcover_plots-processed.csv')