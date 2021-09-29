# updated 2021-09-29
# preliminary analyses for releve cover/richness
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
# despite estimating cover to species in releves, I do not think we should analyze cover at the releve level.
# these estimates were very challenging over such a large area, and are going to be more reliable at the subsampled
# plot level. Instead, releve cover should be primarily used for richness.

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
cover$Unk <- colSums(fxgrps[fxgrps$funcgroup=="Unk",12:31])

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
relcover$Unk<-cover$Unk/cover$Total



#### Richness ####
cover$rich <- specnumber(mtrx)

# histogram of richness
hist((cover$rich))

# calculate stand means and CIs
sumSE.rich <- summarySE(cover, measurevar = "rich", groupvars = c('site','stand')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-releves/richness-releve_sampledata.pdf"), height = 2.5, width = 4)
ggplot(data = cover) +
  geom_line(aes(x = stand, y = rich, group = block), alpha = 0.2)+
  geom_point(aes(x = stand, y = rich, shape = site), alpha = 0.2, size = 1.5) +
  geom_errorbar(data = sumSE.rich, aes(x = stand, ymin = rich - ci, ymax = rich + ci, col = stand), width = 0) +
  geom_point(data = sumSE.rich, aes(x = stand, y = rich, col = stand), size = 3) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank()) +
  labs(x = 'Stand', y = 'Richness (relevés)') +
  guides(shape = FALSE) +
  facet_grid(.~ site,scales = 'free')
dev.off()

# model richness
m.rich1 <- lmer(rich ~ stand + (1|site/block), data = cover)
summary(m.rich1)
AICc(m.rich1)
plot(m.rich1)
# singular fit issue coming from block... remove from model

m.rich1 <- lmer(rich ~ stand + (1|site), data = cover)
summary(m.rich1)

# calculate modeled means and confidence intervals
emm.rich1 <- data.frame(emmeans(m.rich1, ~ stand, type = 'response'))
emm.rich1

# make a table of the results
tab_model(m.rich1, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

# plot the modeled estimates! 
pdf(paste0(Rout,"/cover-releves/richness-releve_modeled.pdf"), height = 2.5, width = 2.5)
ggplot(data = cover) +
  geom_line(aes(x = stand, y = rich, group = block), alpha = 0.2)+
  geom_point(aes(x = stand, y = rich, shape = site), alpha = 0.2, size = 1.5) +
  geom_errorbar(data = emm.rich1, aes(x = stand, ymin = lower.CL, ymax = upper.CL, col = stand), width = 0) +
  geom_point(data = emm.rich1, aes(x = stand, y = emmean, col = stand), size = 3.5) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = '', y = 'Richness (relevés)') +
  guides(col = FALSE)
dev.off()

#### Simpson's Diversity ####
# (probably not reliable for releves since it depends on cover estimates) 
cover$diversity <- diversity(mtrx, index = "simpson")

# histogram of diversity
hist((cover$diversity))
hist(logit(cover$diversity))

# calculate stand means and CIs
sumSE.diversity <- summarySE(cover, measurevar = "diversity", groupvars = c('site','stand')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-releves/diversity-releve_sampledata.pdf"), height = 2.5, width = 4)
ggplot(data = cover) +
  geom_line(aes(x = stand, y = diversity, group = block), alpha = 0.2)+
  geom_point(aes(x = stand, y = diversity, shape = site), alpha = 0.2, size = 1.5) +
  geom_errorbar(data = sumSE.diversity, aes(x = stand, ymin = diversity - ci, ymax = diversity + ci, col = stand), width = 0) +
  geom_point(data = sumSE.diversity, aes(x = stand, y = diversity, col = stand), size = 3) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank()) +
  labs(x = 'Stand', y = "Simpson's Diversity (relevés)") +
  guides(shape = FALSE) +
  facet_grid(.~ site,scales = 'free')
dev.off()

# model diversity
m.diversity1 <- lmer(logit(diversity) ~ stand + (1|site/block), data = cover)
summary(m.diversity1)
AICc(m.diversity1)
plot(m.diversity1)
# singular fit issue coming from site... remove from model

m.diversity1 <- lmer(logit(diversity) ~ stand + (1|block), data = cover)
summary(m.diversity1)

# calculate modeled means and confidence intervals
emm.diversity1 <- data.frame(emmeans(m.diversity1, ~ stand, type = 'response'))
emm.diversity1

# make a table of the results
tab_model(m.diversity1, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

# plot the modeled estimates! 
pdf(paste0(Rout,"/cover-releves/diversity-releve_modeled.pdf"), height = 2.5, width = 2.5)
ggplot(data = cover) +
  geom_line(aes(x = stand, y = diversity, group = block), alpha = 0.2)+
  geom_point(aes(x = stand, y = diversity, shape = site), alpha = 0.2, size = 1.5) +
  geom_errorbar(data = emm.diversity1, aes(x = stand, ymin = lower.CL, ymax = upper.CL, col = stand), width = 0) +
  geom_point(data = emm.diversity1, aes(x = stand, y = response, col = stand), size = 3.5) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = '', y = 'Diversity (relevés)') +
  guides(col = FALSE)
dev.off()

#### Evenness ####
cover$dominance <- 1 - cover$diversity # Simpson's dominance
cover$invsimp <- 1/cover$dominance
cover$evenness <- cover$invsimp/cover$rich

# histogram of evenness
hist((cover$evenness))
hist(logit(cover$evenness))

# calculate stand means and CIs
sumSE.evenness <- summarySE(cover, measurevar = "evenness", groupvars = c('site','stand')) # aggregating to stand, block, site

pdf(paste0(Rout,"/cover-releves/evenness-releve_sampledata.pdf"), height = 2.5, width = 4)
ggplot(data = cover) +
  geom_line(aes(x = stand, y = evenness, group = block), alpha = 0.2)+
  geom_point(aes(x = stand, y = evenness, shape = site), alpha = 0.2, size = 1.5) +
  geom_errorbar(data = sumSE.evenness, aes(x = stand, ymin = evenness - ci, ymax = evenness + ci, col = stand), width = 0) +
  geom_point(data = sumSE.evenness, aes(x = stand, y = evenness, col = stand), size = 3) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank()) +
  labs(x = 'Stand', y = "Simpson's evenness (relevés)") +
  guides(shape = FALSE) +
  facet_grid(.~ site,scales = 'free')
dev.off()

# model evenness
m.evenness1 <- lmer(logit(evenness) ~ stand + (1|site/block), data = cover)
summary(m.evenness1)
AICc(m.evenness1)
plot(m.evenness1)

# calculate modeled means and confidence intervals
emm.evenness1 <- data.frame(emmeans(m.evenness1, ~ stand, type = 'response'))
emm.evenness1

# make a table of the results
tab_model(m.evenness1, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)

# plot the modeled estimates! 
pdf(paste0(Rout,"/cover-releves/evenness-releve_modeled.pdf"), height = 2.5, width = 2.5)
ggplot(data = cover) +
  geom_line(aes(x = stand, y = evenness, group = block), alpha = 0.2)+
  geom_point(aes(x = stand, y = evenness, shape = site), alpha = 0.2, size = 1.5) +
  geom_errorbar(data = emm.evenness1, aes(x = stand, ymin = lower.CL, ymax = upper.CL, col = stand), width = 0) +
  geom_point(data = emm.evenness1, aes(x = stand, y = response, col = stand), size = 3.5) +
  # geom_hline(yintercept = 0.1732621, col = 'goldenrod1') +
  # geom_hline(yintercept = 0.2751764, col = 'darkolivegreen4') +
  scale_color_manual(values = c("goldenrod1","darkolivegreen4")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = '', y = 'Evenness (relevés)') +
  guides(col = FALSE)
dev.off()

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

# pdf(paste0(Rout,"/cover-releves/relcoverDUR-releve_sampledata.pdf"), height = 4, width = 4)
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
  labs(x = 'Block', y = 'Relative cover (relevés)') +
  guides(shape = FALSE) +
  facet_grid(stand2 ~ site, scales = 'free')
# dev.off()

#### cover broken down by species ####
cov.long.spp <- cbind(site=cov.long$site,cov.long[,11:116])
cov.long.spp <- aggregate(. ~ site, data = cov.long.spp, FUN = sum, na.rm = TRUE)

# go from wide to long
cov.spp <- gather(cov.long.spp, species, cover, Achillea.millefolium:Wyethia.angustifolia, factor_key=TRUE)

# sort by site and species
cov.spp <- cov.spp[order(cov.spp$site, cov.spp$species),]

# append fxgrps data to it
cov.spp <- cbind(cov.spp[,-2], rbind(fxgrps[,1:11],fxgrps[,1:11]))

# plot species ordered by cover
# cov.spp$funcgroup <- factor(cov.spp$funcgroup, 
#                             levels = c('NPF','IPF','NAF','IBF','IAF','NPG','IPG','IAG','Woody','Unk'))
cov.spp$funcgroup <- factor(cov.spp$funcgroup, 
                            levels = c('NPF','NPG','NAF','IPG','IPF','IAG','IBF','Woody','IAF','Unk'))

lgd <- ggplot(data = cov.spp[cov.spp$site=='Pisgah' & !cov.spp$cover==0,], aes(x = reorder(species_complete, cover),y = cover, fill = funcgroup)) +
  geom_bar(stat='identity',width = 0.5) +
  scale_fill_manual(values = c('goldenrod3','darkgreen','goldenrod1','palegreen3','mediumpurple4','palegreen1','mediumpurple3','chocolate4','mediumpurple1','gray60'))+
  theme_bw()+theme(legend.position='top',
                   axis.text.y = element_text(face='italic'))+
  labs(x='',fill='',y='Total Cover')+
  scale_y_log10() +
  facet_grid(.~site) +
  coord_flip()
lgd<-get_legend(lgd)
lgd.grid <- plot_grid(lgd)
lgd.grid

sp.pisgah <- ggplot(data = cov.spp[cov.spp$site=='Pisgah' & !cov.spp$cover==0,], aes(x = reorder(species_complete, cover),y = cover, fill = funcgroup)) +
  geom_bar(stat='identity',width = 0.5) +
  scale_fill_manual(values = c('goldenrod3','darkgreen','goldenrod1','palegreen3','mediumpurple4','palegreen1','mediumpurple3','chocolate4','mediumpurple1','gray60'))+
  theme_bw()+theme(legend.position='none',
                   axis.text.y = element_text(face='italic'))+
  labs(x='',fill='',y='Total Cover')+
  scale_y_log10() +
  facet_grid(.~site) +
  coord_flip()
sp.pisgah

sp.SEM <- ggplot(data = cov.spp[cov.spp$site=='South Eugene Meadows' & !cov.spp$cover==0,], aes(x = reorder(species_complete, cover),y = cover, fill = funcgroup)) +
  geom_bar(stat='identity',width = 0.5) +
  scale_fill_manual(values = c('goldenrod3','darkgreen','goldenrod1','palegreen3','mediumpurple4','palegreen1','mediumpurple3','chocolate4','mediumpurple1','gray60'))+
  theme_bw()+theme(legend.position='none',
                   axis.text.y = element_text(face='italic'))+
  labs(x='',fill='',y='Total Cover')+
  scale_y_log10() +
  facet_grid(.~site) +
  coord_flip()
sp.SEM

sp.plots <- plot_grid(sp.pisgah, sp.SEM)
sp.plots <- plot_grid(lgd.grid, sp.plots, nrow = 2,
                      rel_heights = c(0.1,1))

ggsave(sp.plots, file=paste0(Rout,"/cover-releves/sp-plots.pdf"),height=10,width=10)