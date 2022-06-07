library(vegan)
library(ggplot2)

# run cover-releves script through line 30

# make a matrix for functional groups
mtrx <- cover[,31:38]

# across both sites
bray <- vegdist(mtrx, method = "bray", binary = FALSE)

# #Perform a non-linear ordination (Non-Metric Multidimensional Scaling)
ord <- metaMDS(mtrx, distance = "bray", binary= FALSE, k=2, trymax=100)
ord

# You can investigate the species which may be driving the site distribution pattern, referred to as intrinsic variables.
spp.fit <- envfit(ord, mtrx, permutations = 999)
head(spp.fit)

ordiplot(ord, type = "n", main = "Functional group ordination", xlim = c(-1,1), ylim = c(-0.5,0.5))
points(ord, display = "site", pch = c(16, 8) [as.factor(cover$site)], col = c("goldenrod1","darkolivegreen4") [as.factor(cover$stand)])
ordiellipse(ord, conf = 0.95, groups = cover$stand, label = TRUE, draw = "polygon", cex=0.6)
plot(spp.fit, p.max = 1, col = "grey", cex = 0.7) # change the significance level of species shown with p.max
plot(spp.fit, p.max = 0.05, col = "black", cex = 0.7) # change the significance level of species shown with p.max
legend("topright", legend = c(levels(as.factor(cover$site)), levels(as.factor(cover$stand))), pch = c(16, 8, 16, 16), col = c('black','black',"goldenrod1","darkolivegreen4"), bty = "n", cex = 0.8) # displays symbol and colour legend

#### now do with all species
# run 'processing-cover-releves.R' through line 30
mtrx2 <- cov.long[,8:116]
colnames(mtrx2)[1] <- 'litter'

# across both sites
bray2 <- vegdist(mtrx2, method = "bray", binary = FALSE)

# #Perform a non-linear ordination (Non-Metric Multidimensional Scaling)
ord2 <- metaMDS(mtrx2, distance = "bray", binary= FALSE, k=2, trymax=100)
ord2

# You can investigate the species which may be driving the site distribution pattern, referred to as intrinsic variables.
spp.fit2 <- envfit(ord2, mtrx2, permutations = 999)
head(spp.fit2)

ordiplot(ord2, type = "n", main = "Species ordination (p < 0.01)")
points(ord2, display = "site", pch = c(16, 8) [as.factor(cover$site)], col = c("goldenrod1","darkolivegreen4") [as.factor(cover$stand)])
ordiellipse(ord2, conf = 0.95, groups = cover$stand, label = TRUE, draw = "polygon", cex=0.6)
# plot(spp.fit2, p.max = 1, col = "grey", cex = 0.7) # change the significance level of species shown with p.max
plot(spp.fit2, p.max = 0.01, col = "black", cex = 0.7) # change the significance level of species shown with p.max
legend("topright", legend = c(levels(as.factor(cover$site)), levels(as.factor(cover$stand))), pch = c(16, 8, 16, 16), col = c('black','black',"goldenrod1","darkolivegreen4"), bty = "n", cex = 0.8) # displays symbol and colour legend
