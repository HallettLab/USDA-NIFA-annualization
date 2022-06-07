# updated 2022-04-06
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

# subset data to large releve plots
cov.long <- subset(cov.long, is.na(transect))
cov.wide <- cov.wide[,c(1,2,6,10,14,18,22,26,30,34,38,
                        42,46,50,54,58,62,66,70,74,78)]

# large releve plots were 100m2. there was 1 releve per stand (meaning 2 releves per block; 20 total across sites).
# within releves, I estimated cover to species level, estimated litter and bareground. did not measure litter depth.

# create a matrix of only the species as columns
mtrx <- cov.long[,11:116]

# read in functional group metadata to group species
fxgrps <- read.csv('functional_groups.csv')

# prep a cover dataframe with metadata; factor all the categorical variables
cover <- lapply(cov.long[,1:6], factor) %>%
  as.data.frame(.)

# prep functional group dataframe
fxgrps <- cbind(fxgrps, cov.wide[11:116,2:21])
colnames(fxgrps)[12:31] <- cov.wide[3,2:21]

# make cover measurements numeric
fxgrps[,12:31] <- sapply(fxgrps[,12:31],as.numeric)

#### calculate total and functional group covers ####
cover <- cbind(cover, cov.long[,8:10])
cover$Total <- colSums(fxgrps[,12:31])
cover$NATNative <- colSums(fxgrps[fxgrps$nativity=="Native",12:31])
cover$NATIntroduced <- colSums(fxgrps[fxgrps$nativity=="Introduced",12:31])
cover$NATUnknown <- colSums(fxgrps[fxgrps$nativity=="Unknown",12:31])
cover$DURAnnual <- colSums(fxgrps[fxgrps$duration2=="Annual",12:31])
cover$DURPerennial <- colSums(fxgrps[fxgrps$duration2=="Perennial",12:31])
cover$DURBiennial <- colSums(fxgrps[fxgrps$duration2=="Biennial",12:31])
cover$DURUnknown <- colSums(fxgrps[fxgrps$duration2=="Unknown",12:31])
cover$FORGraminoid <- colSums(fxgrps[fxgrps$form2=="Graminoid",12:31])
cover$FORForb <- colSums(fxgrps[fxgrps$form2=="Forb",12:31])
cover$FORWoody <- colSums(fxgrps[fxgrps$form2=="Woody",12:31])
cover$NatAnn <- colSums(fxgrps[fxgrps$nativity=="Native"&fxgrps$duration2=="Annual",12:31])
cover$NatPer <- colSums(fxgrps[fxgrps$nativity=="Native"&fxgrps$duration2=="Perennial",12:31])
cover$IntAnn <- colSums(fxgrps[fxgrps$nativity=="Introduced"&fxgrps$duration2=="Annual",12:31])
cover$IntPer <- colSums(fxgrps[fxgrps$nativity=="Introduced"&fxgrps$duration2=="Perennial",12:31])
cover$AG <- colSums(fxgrps[fxgrps$form2=="Graminoid"&fxgrps$duration2=="Annual",12:31])
cover$PG <- colSums(fxgrps[fxgrps$form2=="Graminoid"&fxgrps$duration2=="Perennial",12:31])
cover$AF <- colSums(fxgrps[fxgrps$form2=="Forb"&fxgrps$duration2=="Annual",12:31])
cover$PF <- colSums(fxgrps[fxgrps$form2=="Forb"&fxgrps$duration2=="Perennial",12:31])
cover$NAG <- colSums(fxgrps[fxgrps$funcgroup=="NAG",12:31])
cover$NAF <- colSums(fxgrps[fxgrps$funcgroup=="NAF",12:31])
cover$NPG <- colSums(fxgrps[fxgrps$funcgroup=="NPG",12:31])
cover$NPF <- colSums(fxgrps[fxgrps$funcgroup=="NPF",12:31])
cover$IAG <- colSums(fxgrps[fxgrps$funcgroup=="IAG",12:31])
cover$IAF <- colSums(fxgrps[fxgrps$funcgroup=="IAF",12:31])
cover$IPG <- colSums(fxgrps[fxgrps$funcgroup=="IPG",12:31])
cover$IPF <- colSums(fxgrps[fxgrps$funcgroup=="IPF",12:31])
cover$Woody <- colSums(fxgrps[fxgrps$funcgroup=="Woody",12:31])

#### calculate relative covers ####
relcover<-cover[,1:6]
relcover$NATNative<-cover$NATNative/cover$Total
relcover$NATIntroduced<-cover$NATIntroduced/cover$Total
relcover$NATUnknown<-cover$NATUnknown/cover$Total
relcover$DURAnnual<-cover$DURAnnual/cover$Total
relcover$DURPerennial<-cover$DURPerennial/cover$Total
relcover$DURBiennial<-cover$DURBiennial/cover$Total
relcover$DURUnknown<-cover$DURUnknown/cover$Total
relcover$FORGraminoid<-cover$FORGraminoid/cover$Total
relcover$FORForb<-cover$FORForb/cover$Total
relcover$FORWoody<-cover$FORWoody/cover$Total
relcover$NatAnn<-cover$NatAnn/cover$Total
relcover$NatPer<-cover$NatPer/cover$Total
relcover$IntAnn<-cover$IntAnn/cover$Total
relcover$IntPer<-cover$IntPer/cover$Total
relcover$NAG<-cover$NAG/cover$Total
relcover$NAF<-cover$NAF/cover$Total
relcover$NPG<-cover$NPG/cover$Total
relcover$NPF<-cover$NPF/cover$Total
relcover$IAG<-cover$IAG/cover$Total
relcover$IAF<-cover$IAF/cover$Total
relcover$IPG<-cover$IPG/cover$Total
relcover$IPF<-cover$IPF/cover$Total
relcover$Woody<-cover$Woody/cover$Total



#### calculate various richness measurements, diversity, etc. ####
cover$rich <- specnumber(mtrx)

# native richness
mtrx.nat <- subset(fxgrps, nativity =='Native') %>%
  t(.) %>%
  .[12:31,] %>%
  as.data.frame(.) %>%
  sapply(.,as.numeric) %>%
  as.data.frame(.)
cover$rich.nat <- specnumber(mtrx.nat)

# forb richness
mtrx.forb <- subset(fxgrps, form =='Forb') %>%
  t(.) %>%
  .[12:31,] %>%
  as.data.frame(.) %>%
  sapply(.,as.numeric) %>%
  as.data.frame(.)
cover$rich.forb <- specnumber(mtrx.forb)

# native forb richness
mtrx.natforb <- subset(fxgrps, nativity == 'Native' & form =='Forb') %>%
  t(.) %>%
  .[12:31,] %>%
  as.data.frame(.) %>%
  sapply(.,as.numeric) %>%
  as.data.frame(.)
cover$rich.natforb <- specnumber(mtrx.natforb)

cover$diversity <- diversity(mtrx, index = "simpson")
cover$dominance <- 1 - cover$diversity # Simpson's dominance
cover$invsimp <- 1/cover$dominance
cover$evenness <- cover$invsimp/cover$rich


######## save processed cover data as .csv files ########
write.csv(cover, file = 'community_cover_releves-processed.csv')
write.csv(relcover, file = 'community_relcover_releves-processed.csv')


#### Annual:Perennial ratio ####
relcover.DUR <- gather(relcover, Duration, value, DURAnnual:DURUnknown, factor_key = TRUE)
relcover.DUR$Duration <- gsub( "DUR", "", as.character(relcover.DUR$Duration))
relcover.DUR$Duration <- factor(relcover.DUR$Duration, levels = c('Annual', 'Perennial', 'Biennial', 'Unknown'))

# sumSE.DUR <- summarySE(relcover.DUR, measurevar="value", groupvars = c('site','stand','Duration'))
# sumSE.DUR$Duration <- factor(sumSE.DUR$Duration, levels = c('Annual', 'Perennial', 'Biennial', 'Unknown'))
# 
# ggplot(sumSE.DUR, aes(x = stand, y = value, fill = Duration))+
#   scale_fill_manual(values = c("goldenrod1","darkolivegreen4",'blue','grey')) +
#   geom_bar(stat = "identity")+
#   scale_y_continuous(expand=c(0,0),limits=c(0,1),breaks=c(0,0.25,0.5,0.75,1))+
#   labs(x="Stand",y="Relative cover",fill="")+
#   theme_bw() +
#   theme(legend.position = 'top') +
#   facet_grid(.~ site,scales = 'free')
relcover.DUR$stand2 <- relcover.DUR$stand
levels(relcover.DUR$stand2) <- c("Annual stand", "Perennial stand")

pdf(paste0(Rout,"/cover-releves/relcover-releves.pdf"), height = 4, width = 4)
ggplot(data = relcover.DUR, aes(x = block, y = value, fill = Duration, group = stand2)) +
  geom_bar(aes(x = block, y = value, group = stand2, fill = Duration), stat = "identity")+
  scale_fill_manual(values = c("goldenrod1","darkolivegreen4",'blue','grey')) +
  scale_y_continuous(expand=c(0,0),limits=c(0,1),breaks=c(0,0.25,0.5,0.75)) +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  labs(x = 'Block', y = 'Relative cover (100 m2 relevés)') +
  guides(shape = FALSE) +
  facet_grid(stand2 ~ site, scales = 'free')
dev.off()

#### cover broken down by species ####
cov.long.spp <- cbind(site=cov.long$site,cov.long[,11:116])
cov.long.spp <- aggregate(. ~ site, data = cov.long.spp, FUN = sum, na.rm = TRUE)

# go from wide to long
cov.spp <- gather(cov.long.spp, species, cover, Achillea.millefolium:Wyethia.angustifolia, factor_key=TRUE)

# sort by site and species
cov.spp <- cov.spp[order(cov.spp$site, cov.spp$species),]

# append fxgrps data to it
cov.spp <- cbind(cov.spp[,-2], rbind(fxgrps[,1:11],fxgrps[,1:11]))

# give any species with cover = 1 a value of cover = 1.5 so it plots better on log scale
cov.spp[cov.spp$cover=='1',]$cover <- 1.5
# plot species ordered by cover
# cov.spp$funcgroup <- factor(cov.spp$funcgroup, 
#                             levels = c('NPF','IPF','NAF','IBF','IAF','NPG','IPG','IAG','Woody','Unk'))
cov.spp$funcgroup <- factor(cov.spp$funcgroup, 
                            levels = c('NPF','NAF','IPF','IAF','IBF','??F','NPG','IPG','I?G','IAG','Woody'))

lgd <- ggplot(data = cov.spp[cov.spp$site=='Pisgah' & !cov.spp$cover==0,], aes(x = reorder(species_complete, cover),y = cover, fill = funcgroup)) +
  geom_bar(stat='identity',width = 0.5) +
  scale_fill_manual(values = c('goldenrod3','goldenrod1','mediumpurple4','mediumpurple3','mediumpurple1','mediumorchid1','darkgreen','palegreen3','darkolivegreen1','palegreen1','chocolate4'))+
  theme_bw()+theme(legend.position='top',
                   axis.text.y = element_text(face='italic'),
                   legend.direction = 'horizontal')+
  labs(x='',fill='',y='Total Cover')+
  scale_y_log10() +
  facet_grid(.~site) +
  guides(fill = guide_legend(ncol = 11)) +
  coord_flip()
lgd<-get_legend(lgd)
lgd.grid <- plot_grid(lgd)
lgd.grid

sp.pisgah <- ggplot(data = cov.spp[cov.spp$site=='Pisgah' & !cov.spp$cover==0,], aes(x = reorder(species_complete, cover),y = cover, fill = funcgroup)) +
  geom_bar(stat='identity',width = 0.5) +
  scale_fill_manual(values = c('goldenrod3','goldenrod1','mediumpurple4','mediumpurple3','mediumpurple1','mediumorchid1','darkgreen','palegreen3','darkolivegreen1','palegreen1','chocolate4'))+
  theme_bw()+theme(legend.position='none',
                   axis.text.y = element_text(face='italic'))+
  labs(x='',fill='',y='Total Cover')+
  scale_y_log10() +
  facet_grid(.~site) +
  coord_flip()
sp.pisgah

sp.SEM <- ggplot(data = cov.spp[cov.spp$site=='South Eugene Meadows' & !cov.spp$cover==0,], aes(x = reorder(species_complete, cover),y = cover, fill = funcgroup)) +
  geom_bar(stat='identity',width = 0.5) +
  scale_fill_manual(values = c('goldenrod3','goldenrod1','mediumpurple4','mediumpurple3','mediumpurple1','darkgreen','palegreen3','darkolivegreen1','palegreen1','chocolate4'))+
  theme_bw()+theme(legend.position='none',
                   axis.text.y = element_text(face='italic'))+
  labs(x='',fill='',y='Total Cover')+
  scale_y_log10() +
  facet_grid(.~site) +
  coord_flip()
sp.SEM

sp.plots <- plot_grid(sp.pisgah, sp.SEM)
sp.plots <- plot_grid(lgd.grid, sp.plots, nrow = 2,
                      rel_heights = c(0.05,1))

ggsave(sp.plots, file=paste0(Rout,"/cover-releves/sp-releves.pdf"),height=10,width=10)
