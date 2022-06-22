# gSSURGO database table processing
library(ggplot2)
library(plyr)
library(dplyr)

# load in muaggatt table; this is what we're striving for: variables of interest aggregated to the mapunit level
muaggatt <- read.csv('muaggatt.csv')

# load in chorizon table
chorizon <- read.csv('chorizon.csv')
chorizon2 <- chorizon[,c(1:11,171)]

# load in component table
component <- read.csv('component.csv')

# load in chtexturegrp table
chtexturegrp <- read.csv('chtexturegrp.csv')

#### soil depth ####
# for mapunits with brockdepmin records, use that variable from muaggatt table.
# for other mapunits, determine depth by choosing deepest nonrock horizon of each component (i.e., max hzdept_r from chorizon), and then averaging these numbers weighted by component percent
hzdept_r_max <- aggregate(hzdept_r ~ cokey, chorizon, FUN = max)

hzdept_r_max <- merge(component, hzdept_r_max, by= 'cokey')

depth <- ddply(hzdept_r_max, .(mukey),
      function(x) data.frame(deepest = weighted.mean(x$hzdept_r, x$comppct_r)))

depth2 <- muaggatt[,c(41,7)]

depth <- merge(depth, depth2, by= 'mukey')

ggplot(depth, aes(x = brockdepmin, y = deepest)) +
  geom_point() +
  geom_smooth(method = 'lm')
summary(lm(deepest ~ brockdepmin, data = depth))

depth$depth <- depth$brockdepmin
depth[is.na(depth$depth),]$depth <- depth[is.na(depth$depth),]$deepest

#### albedo ####
albedo <- ddply(component, .(mukey),
                function(x) data.frame(albedo = weighted.mean(x$albedodry_r, x$comppct_r)))

#### aspect ####
aspect <- ddply(component, .(mukey),
                function(x) data.frame(aspect = weighted.mean(x$aspectrep, x$comppct_r)))

#### elevation ####
elevation <- ddply(component, .(mukey),
                function(x) data.frame(elevation = weighted.mean(x$elev_r, x$comppct_r)))

#### OM ####
# average attribute of topmost horizons (where hzdept_r = 0), weighted by component percent
chorizon_top <- subset(chorizon, hzdept_r == 0)

OM <- data.frame(cokey = chorizon_top$cokey,
                 om_r = cbind(chorizon_top$om_r))
OM <- merge(component, OM, by= 'cokey')

OM <- ddply(OM, .(mukey),
                function(x) data.frame(OM = weighted.mean(x$om_r, x$comppct_r)))

#### pH ####
# average attribute of topmost horizons (where hzdept_r = 0), weighted by component percent
pH <- data.frame(cokey = chorizon_top$cokey,
                 ph1to1h2o_r = cbind(chorizon_top$ph1to1h2o_r))
pH <- merge(component, pH, by= 'cokey')

pH <- ddply(pH, .(mukey),
            function(x) data.frame(pH = weighted.mean(x$ph1to1h2o_r, x$comppct_r)))

#### caco3 ####
# average attribute of topmost horizons (where hzdept_r = 0), weighted by component percent
caco3 <- data.frame(cokey = chorizon_top$cokey,
                 caco3_r = cbind(chorizon_top$caco3_r))
caco3 <- merge(component, caco3, by= 'cokey')

caco3 <- ddply(caco3, .(mukey),
            function(x) data.frame(caco3 = weighted.mean(x$caco3_r, x$comppct_r)))

#### phosphorus ####
# # average attribute of topmost horizons (where hzdept_r = 0), weighted by component percent
# pbray1 <- data.frame(cokey = chorizon_top$cokey,
#                     pbray1_r = cbind(chorizon_top$pbray1_r))
# pbray1 <- merge(component, pbray1, by= 'cokey')
# 
# pbray1 <- ddply(pbray1, .(mukey),
#                function(x) data.frame(pbray1 = weighted.mean(x$pbray1_r, x$comppct_r)))

#### percent sand, silt, and clay ####
# use sandtotal_r for topmost horizon in the dominant component of the mapunit.
textures <- data.frame(cokey = chorizon_top$cokey,
                   chkey = chorizon_top$chkey,
                   sandtotal_r = cbind(chorizon_top$sandtotal_r),
                   silttotal_r = cbind(chorizon_top$silttotal_r),
                   claytotal_r = cbind(chorizon_top$claytotal_r))

textures <- merge(component, textures, by = 'cokey')

textures <- textures %>% 
  group_by(mukey) %>% 
  slice(which.max(comppct_r))

textures <- data.frame(cokey = textures$cokey,
                       chkey = textures$chkey,
                       mukey = textures$mukey,
                       sandtotal_r = textures$sandtotal_r,
                       silttotal_r = textures$silttotal_r,
                       claytotal_r = textures$claytotal_r)

#### mapunit_agg compiled table ####
mapunit_agg <- muaggatt[,c(1:3,41,5:6,13:16,7)]

mapunit_agg$depth <- depth$deepest[match(mapunit_agg$mukey, depth$mukey)]
mapunit_agg$depth2 <- mapunit_agg$depth
mapunit_agg[!is.na(mapunit_agg$brockdepmin),]$depth2 <- mapunit_agg[!is.na(mapunit_agg$brockdepmin),]$brockdepmin
mapunit_agg$albedo <- albedo$albedo[match(mapunit_agg$mukey, albedo$mukey)]
mapunit_agg$aspect <- aspect$aspect[match(mapunit_agg$mukey, aspect$mukey)]
mapunit_agg$elevation <- elevation$elevation[match(mapunit_agg$mukey, elevation$mukey)]
mapunit_agg$OM <- OM$OM[match(mapunit_agg$mukey, OM$mukey)]
mapunit_agg$pH <- pH$pH[match(mapunit_agg$mukey, pH$mukey)]
mapunit_agg$caco3 <- caco3$caco3[match(mapunit_agg$mukey, caco3$mukey)]
# mapunit_agg$pbray1 <- pbray1$pbray1[match(mapunit_agg$mukey, pbray1$mukey)]
mapunit_agg$sand <- textures$sandtotal_r[match(mapunit_agg$mukey, textures$mukey)]
mapunit_agg$silt <- textures$silttotal_r[match(mapunit_agg$mukey, textures$mukey)]
mapunit_agg$clay <- textures$claytotal_r[match(mapunit_agg$mukey, textures$mukey)]

write.csv(mapunit_agg,file = 'mapunit_variables.csv')

