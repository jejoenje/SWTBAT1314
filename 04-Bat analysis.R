###
###
### ANALYSE BAT DATA 2013+2014
### Avoid all edits to data file, please do in previous script.

library(lme4)

### Load bat master data:
bats <- read.csv('data/BAT DATA 2013 and 2014 MASTER.csv',header=T)

### Load pre-processed habitat data:
hab <- read.csv('data/TRANSECT SECTION HABITAT DATA 2013 and 2014 PROCESSED.csv', header=T)

### Ignore some specific transects/sections that weren't done consistently/only once:
bats <- bats[!(bats$SITE=='Redlands' & bats$TRSCT=='E'),]
bats <- bats[!(bats$SITE=='Letham Far North T2' & bats$TRSCT=='NW' & bats$SECTION=='5'),]
bats <- bats[!(bats$SITE=='Mid Cambushinnie' & bats$TRSCT=='S' & bats$SECTION=='5'),]
bats <- bats[!(bats$SITE=='Nisbet Hill' & bats$TRSCT=='N' & bats$SECTION=='5'),]
bats <- bats[!(bats$SITE=='Wester Essendy' & bats$TRSCT=='W' & bats$SECTION=='4'),]
bats <- droplevels(bats)

### Remove spurious AREA_M2 and AREA_ha in prep for match from habitat data.
bats$AREA_M2 <- NULL
bats$AREA_ha <- NULL

### Merge pre-processed habitat data to bat data:
bats$id <- paste(bats$SITE, bats$TRSCT, bats$SECTION, sep='-')
bats <- merge(bats, hab, 'id')

mod1 <- glmer(PASSES ~ fSECTION + (1|SITE/TRSCT) + offset(log(AREA_ha)), 
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

