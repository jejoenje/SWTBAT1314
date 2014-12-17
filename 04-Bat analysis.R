###
###
### ANALYSE BAT DATA 2013+2014
### Avoid all edits to data file, please do in previous script.

library(lme4)
library(arm)
library(glmmADMB)
library(gstat)
library(sp)

source('../../../000_R/RSimExamples/helpjm.r')

### Load bat master data:
bats <- read.csv('data/BAT DATA 2013 and 2014 MASTER.csv',header=T)

### Check factors, etc:
is.factor(bats$SITE)
is.factor(bats$TRSCT)
is.factor(bats$fSECTION); bats$fSECTION <- factor(bats$fSECTION)
is.factor(bats$NOTURB)

### Proportion of zero's per site:
tapply(bats$ALL_PIPS, bats$SITE, function(x) sum(x==0)/length(x))
### Average proportion of zero's per site:
mean(tapply(bats$ALL_PIPS, bats$SITE, function(x) sum(x==0)/length(x)))

### Try fitting intercept-only NB GLM per site:
sitedat <- bats[bats$SITE==levels(bats$SITE)[2],]
mod_site <- glm.nb(ALL_PIPS ~ fSECTION+WINDS, data=sitedat)
dispZuur(mod_site)
plot(fitted(mod_site), resid(mod_site))


system.time(mod1 <- glmer(PASSES ~ fSECTION + (1|SITE/TRSCT) + offset(log(AREA_ha)), 
              family='poisson', data=bats))

system.time(mod2 <- glmmadmb(PASSES ~ fSECTION + (1|SITE/TRSCT) + offset(log(AREA_ha)), 
                 family='nbinom1', data=bats))

system.time(mod3 <- glmmadmb(PASSES ~ fSECTION + (1|SITE/TRSCT) + offset(log(AREA_ha)), 
                             family='nbinom1', zeroInfl=T, data=bats))

system.time(mod4 <- glmmadmb(PASSES ~ fSECTION + (1|SITE/TRSCT) + offset(log(AREA_ha)), 
                             family='poisson', data=bats))

system.time(mod5 <- glmmadmb(PASSES ~ fSECTION + (1|SITE/TRSCT) + offset(log(AREA_ha)), 
                             family='poisson', zeroInfl=T, data=bats))

AIC(mod1, mod2, mod3, mod4, mod5)
plot(fitted(mod2),resid(mod2))
plot(bats$fSECTION,resid(mod2))
plot(bats$SITE,resid(mod2))
plot(fitted(mod2),bats$PASSES); abline(a=0, b=1, col='red')

system.time(mod2a <- glmmadmb(PASSES ~ fSECTION + NOTURB + TTMIDN + (1|SITE/TRSCT) + offset(log(AREA_ha)), 
                             family='nbinom', data=bats))


bats_nona <- bats[!is.na(bats$WINDS),]; bats_nona <- droplevels(bats_nona)
system.time(mod10 <- glmmadmb(ALL_PIPS ~ fSECTION*NOTURB + 
                   WINDS + DAYNO + TTMIDN + I(TTMIDN^2) + pTREE + D_LIN + 
                   (1|SITE/TRSCT) + offset(log(AREA_ha)), 
                 family='nbinom', data=bats_nona))
plot(fitted(mod10),resid(mod10))
plot(bats_nona$fSECTION,resid(mod10))
plot(bats_nona$NOTURB,resid(mod10))
plot(bats_nona$WINDS,resid(mod10))
plot(bats_nona$DAYNO,resid(mod10))
plot(bats_nona$TTMIDN,resid(mod10))
plot(bats_nona$pTREE,resid(mod10))
plot(bats_nona$D_LIN,resid(mod10))
plot(bats_nona$SITE,resid(mod10))

bats_nona$resid10 <- resid(mod10)
plotbubble(bats_nona$MID_X, bats_nona$MID_Y, bats_nona$resid10)

par(mfrow=c(5,7))
par(mar=c(0,0,0,0))
for(i in 1:nlevels(bats_nona$SITE)) {
  sitedat <- bats_nona[bats_nona$SITE==levels(bats_nona$SITE)[i],]
  plotbubble(sitedat$MID_X, sitedat$MID_Y, sitedat$resid10, axt='n')
}
for(i in 1:nlevels(bats_nona$SITE)) {
  sitedat <- bats_nona[bats_nona$SITE==levels(bats_nona$SITE)[i],]
  coordinates(sitedat) <- c('MID_X','MID_Y')
  vario_site <- variogram(resid10 ~ 1, data=sitedat)
  plotvario(vario_site, axes=F)
}


bats_nona$sWINDS <- scale(bats_nona$WINDS)
bats_nona$sDAYNO <- scale(bats_nona$DAYNO)
bats_nona$sTTMIDN <- scale(bats_nona$TTMIDN)
bats_nona$spTREE <- scale(bats_nona$pTREE)
bats_nona$sD_LIN <- scale(bats_nona$D_LIN)
bats_nona$sTTMIDN2 <- bats_nona$sTTMIDN^2
bats_nona$sD_LIN2 <- bats_nona$sD_LIN^2
bats_nona$sWINDS2 <- bats_nona$sWINDS^2
bats_nona$OCC_PIPS <- as.numeric(as.logical(bats_nona$ALL_PIPS))
mod90 <- glmer(OCC_PIPS ~ fSECTION*NOTURB + 
                 sWINDS + sDAYNO + sTTMIDN + spTREE + sD_LIN + 
                 (1|SITE/TRSCT), 
               family=binomial, data=bats_nona)
dispZuur(mod90)
binnedplot(fitted(mod90), resid(mod90, type='pearson'))


system.time(mod11 <- update(mod10, .~. +I(D_LIN^2) +I(WINDS^2)) )
plot(fitted(mod11),resid(mod11))
plot(bats_nona$fSECTION,resid(mod11))
plot(bats_nona$NOTURB,resid(mod11))
plot(bats_nona$WINDS,resid(mod11))
plot(bats_nona$DAYNO,resid(mod11))
plot(bats_nona$TTMIDN,resid(mod11))
plot(bats_nona$pTREE,resid(mod11))
plot(bats_nona$D_LIN,resid(mod11))
plot(bats_nona$SITE,resid(mod11))

AIC(mod10, mod11)



system.time(mod12 <- glmmadmb(ALL_PIPS ~ fSECTION*NOTURB + 
                                sWINDS + sWINDS2 + sDAYNO + sTTMIDN + sTTMIDN2 + 
                                spTREE + sD_LIN + sD_LIN2 + 
                                (1|SITE/TRSCT) + offset(log(AREA_ha)), 
                              family='nbinom', data=bats_nona))
plot(fitted(mod12),resid(mod12))
plot(bats_nona$fSECTION,resid(mod12))
plot(bats_nona$NOTURB,resid(mod12))
plot(bats_nona$WINDS,resid(mod12))
plot(bats_nona$DAYNO,resid(mod12))
plot(bats_nona$TTMIDN,resid(mod12))
plot(bats_nona$pTREE,resid(mod12))
plot(bats_nona$D_LIN,resid(mod12))
plot(bats_nona$SITE,resid(mod12))

AIC(mod10, mod11, mod12)

par(mfrow=c(1,2))
dispZuur(mod1)
plot(fitted(mod1),resid(mod1))
dispZuur(mod2)
plot(fitted(mod2),resid(mod2))

mod3 <- glmer(PASSES ~ fSECTION*NOTURB + 
                WINDS + I(WINDS^2) + DAYNO + TTMIDN + I(TTMIDN^2) + pTREE + D_LIN + 
                (1|SITE/TRSCT) + offset(log(AREA_ha)), 
              family='poisson', data=bats)
mod3a <- update(mod3, .~. -I(WINDS^2))
anova(mod3, mod3a)
# Drop wind speed poly
mod3 <- mod3a
mod3a <- update(mod3, .~. -I(TTMIDN^2))
anova(mod3, mod3a)
# Retain time to midnight poly.
dispZuur(mod3)


display(mod4)
dispZuur(mod4)
mod4a <- update(mod4, .~., data=bats_nona)
plot(fitted(mod4a), resid(mod4a))
plot(bats_nona$fSECTION, resid(mod4a))
plot(bats_nona$NOTURB, resid(mod4a))  # limited variation in 4 turbine sites
plot(bats_nona$WINDS, resid(mod4a))  # Terrible resid pattern with wind speed!
plot(bats_nona$DAYNO, resid(mod4a))
plot(bats_nona$TTMIDN, resid(mod4a))
plot(bats_nona$pTREE, resid(mod4a)) # Not so good in residuals with prop. trees
plot(bats_nona$D_LIN, resid(mod4a)) # Terrible resid pattern with distance to lin. features

plot(bats$WINDS, bats$ALL_PIPS)



library(mgcv)
mod3 <- gamm()
>>>>>>> 5c2af5fe0d56033a23d248deebcb496845bb47b4

par(mfrow=c(1,2))
plot(predict(mod1),resid(mod1))
hist(resid(mod1))
sum(residuals(mod1)^2)/(nrow(bats)-7)

library(glmmADMB)

mod2 <- glmmadmb(PASSES ~ fSECTION + (1|SITE/TRSCT) + offset(log(AREA_ha)), family='poisson', data=bats)
plot(predict(mod2),resid(mod2))
hist(resid(mod2))
sum(residuals(mod2)^2)/(nrow(bats)-7)

mod2a <- glmmadmb(PASSES ~ fSECTION + (1|SITE/TRSCT) + offset(log(AREA_ha)), 
                  family='poisson', zeroInflation=T, data=bats)
plot(predict(mod2a),resid(mod2a))
hist(resid(mod2a))
sum(residuals(mod2a)^2)/(nrow(bats)-7)

mod3 <- glmmadmb(PASSES ~ fSECTION + (1|SITE/TRSCT) + offset(log(AREA_ha)), family='nbinom', data=bats)
plot(predict(mod3),resid(mod3))
hist(resid(mod3))
sum(residuals(mod3)^2)/(nrow(bats)-7)

mod4 <- glmmadmb(PASSES ~ fSECTION + (1|SITE/TRSCT) + offset(log(AREA_ha)), family='nbinom1', data=bats)
plot(predict(mod4),resid(mod4))
hist(resid(mod4))
sum(residuals(mod4)^2)/(nrow(bats)-7)

