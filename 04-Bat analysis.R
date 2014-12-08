###
###
### ANALYSE BAT DATA 2013+2014
### Avoid all edits to data file, please do in previous script.

library(lme4)
library(glmmADMB)

### Load bat master data:
bats <- read.csv('data/BAT DATA 2013 and 2014 MASTER.csv',header=T)

### Check factors, etc:
is.factor(bats$SITE)
is.factor(bats$TRSCT)
is.factor(bats$fSECTION); bats$fSECTION <- factor(bats$fSECTION)
is.factor(bats$NOTURB)


mod1 <- glmer(PASSES ~ fSECTION + (1|SITE/TRSCT) + offset(log(AREA_ha)), 
              family='poisson', data=bats)
mod2 <- glmmadmb(PASSES ~ fSECTION + (1|SITE/TRSCT) + offset(log(AREA_ha)), 
                 family='poisson', data=bats)

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

