###
###
### ANALYSE BAT DATA 2013+2014
### Avoid all edits to data file, please do in previous script.

library(lme4)
library(arm)
library(glmmADMB)

dispZuur <- function(mod, restype='pearson') {
  if(restype=='pearson' | restype=='deviance') {
    if(restype=='pearson') {
      return(sum(resid(mod, type='pearson')^2)/(nrow(model.matrix(mod))-attr(logLik(mod),'df')))  
    }
    if(restype=='deviance') {
      return(sum(resid(mod, type='deviance')^2)/(nrow(model.matrix(mod))-attr(logLik(mod),'df')))  
    }
  } else {
    print('Error - residual type should be either \'pearson\' or \'deviance\'')
  }
}

### Load bat master data:
bats <- read.csv('data/BAT DATA 2013 and 2014 MASTER.csv',header=T)

### Check factors, etc:
is.factor(bats$SITE)
is.factor(bats$TRSCT)
is.factor(bats$fSECTION); bats$fSECTION <- factor(bats$fSECTION)
is.factor(bats$NOTURB)


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

mod4 <- glmer(ALL_PIPS ~ fSECTION*NOTURB + 
                WINDS + DAYNO + TTMIDN + I(TTMIDN^2) + pTREE + D_LIN + 
                (1|SITE/TRSCT) + offset(log(AREA_ha)), 
              family='poisson', data=bats)
display(mod4)
dispZuur(mod4)
bats_nona <- bats[!is.na(bats$WINDS),]; bats_nona <- droplevels(bats_nona)
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

