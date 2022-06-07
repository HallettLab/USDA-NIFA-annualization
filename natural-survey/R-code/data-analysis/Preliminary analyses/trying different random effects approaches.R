pisgah <- subset(cover, site == 'Pisgah')
pisgah$block <- factor(pisgah$block)
pisgah$transect <- factor(pisgah$transect)


sumSE.richness <- summarySE(pisgah, measurevar = "rich", groupvars = 'stand') # just aggregating to stand
sumSE.richness2 <- summarySE(pisgah, measurevar = "rich", groupvars = c('stand','block')) # aggregating to stand, block, site



mod <- lmer(log(rich) ~ stand + (1|block) + (1|transect), data = pisgah)
mod2 <- lmer(log(rich) ~ stand + (1|block), data = pisgah)
mod3 <- lmer(log(rich) ~ stand + (1|transect), data = pisgah)
mod4 <- lm(log(rich) ~ stand, data = pisgah)
mod5 <- lm(log(rich) ~ stand, data = sumSE.richness2)
mod6 <- lmer(log(rich) ~ stand + (1|block), data = sumSE.richness2)


tab_model(mod, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)
tab_model(mod2, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)
tab_model(mod3, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)
tab_model(mod4, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)
tab_model(mod5, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)
tab_model(mod6, linebreak = FALSE, p.val = 'satterthwaite', transform = NULL)
AICc(mod)
AICc(mod2)
AICc(mod3)
AICc(mod4)
AICc(mod5)
AICc(mod6)