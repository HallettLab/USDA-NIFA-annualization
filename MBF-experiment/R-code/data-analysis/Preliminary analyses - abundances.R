#### 6/22/2022: Preliminary analyses MBF experiment
library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(MuMIn)
library(Rmisc)
library(ggplot2)
library(dplyr)
library(sjPlot)
library(tidyr)
options(contrasts = c("contr.sum", "contr.poly"))

# set working directory to cleaned data folder
setwd("C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/MBF Experiment/Reed_USDA-NIFA_MBF Experiment/Data-cleaned")

# set path to R-output folder within Data-analysis
Rout <- "C:/Users/paulr/Dropbox (University of Oregon)/OREGON/Postdoc/USDA-NIFA project/MBF Experiment/Reed_USDA-NIFA_MBF Experiment/R-code/Data-analysis/R-output"

# read in June abundance data
abundance <- read.csv('Abundances_MBF_2022-06.csv')