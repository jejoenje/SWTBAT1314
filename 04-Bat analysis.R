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
bats$fSECTION <- factor(bats$fSECTION)
is.factor(bats$fSECTION)
is.factor(bats$NOTURB)

### Proportion of zero's per site:
tapply(bats$ALL_PIPS, bats$SITE, function(x) sum(x==0)/length(x))
### Average proportion of zero's per site:
mean(tapply(bats$ALL_PIPS, bats$SITE, function(x) sum(x==0)/length(x)))

### First date:
min(as.Date(as.vector(bats$DATE),format='%Y-%m-%d'))

### Last date:
max(as.Date(as.vector(bats$DATE),format='%Y-%m-%d'))

### Binary responses and standardise predictors:
bats$OCC_PIPS <- as.numeric(as.logical(bats$ALL_PIPS))
bats$OCC_OTHR <- as.numeric(as.logical(bats$ALL_OTHER))
bats$sMINTEMP <- scale(bats$min_air_temp)
bats$sMINTEMP2 <- bats$sMINTEMP^2
bats$sWINDS <- scale(bats$WINDS)
bats$sWINDS2 <- bats$sWINDS^2
bats$sDAYNO <- scale(bats$DAYNO)
bats$sDAYNO2 <- bats$sDAYNO^2
bats$sTTMIDN <- scale(bats$TTMIDN)
bats$sTTMIDN2 <- bats$sTTMIDN^2
bats$spBUILD <- scale(bats$pBUILD)
bats$sD_TRE <- scale(bats$D_TRE)

### Create TURB indicator var for "single" or "multiple" (>1) turbines:
bats$TURB <- "single"
bats$TURB[bats$NOTURB>1] <- "multiple"
bats$TURB <- factor(bats$TURB)
bats$TURB <- relevel(bats$TURB, ref="single")

bats_nona <- bats[!is.na(bats$WINDS),]
bats_nona <- droplevels(bats_nona)

###
### OCCURRENCE MODELS (PROBABILITY OF A PASS)
###
### Basic model: GLMM with nested RE:

mod1 <- glmer(OCC_PIPS ~ SECTION*TURB + sMINTEMP +
                sWINDS + 
                sDAYNO + 
                sTTMIDN + 
                spBUILD + 
                sD_TRE + 
                (1|SITE/TRSCT), data=bats_nona, family=binomial)
summary(mod1)
dispZuur(mod1)

### To test, what happens when we remove the TRSCT nesting, and the entire RE:
mod1a <- glmer(OCC_PIPS ~ SECTION*TURB + sMINTEMP +
                 sWINDS + 
                 sDAYNO + 
                 sTTMIDN + 
                 spBUILD + 
                 sD_TRE + 
                 (1|SITE), data=bats_nona, family=binomial)
mod0 <- glm(OCC_PIPS ~ SECTION*TURB + sMINTEMP +
              sWINDS + 
              sDAYNO + 
              sTTMIDN + 
              spBUILD + 
              sD_TRE, data=bats_nona, family=binomial)
dispZuur(mod1)
dispZuur(mod1a)
dispZuur(mod0)
AIC(mod1, mod1a, mod0)

### Check SAC in binomial GLM:
library(ncf)
correlog0 <- correlog(bats_nona$MID_X, bats_nona$MID_Y, residuals(mod0),
                        na.rm=T, increment=1, resamp=0)
plot(correlog0$correlation[1:20], pch=16, type='o'); abline(h=0, lty='dashed')
# make a map of the residuals
plot(bats_nona$MID_X, bats_nona$MID_Y, 
     col=c("blue","red")[sign(resid(mod0))/2+1.5], 
     pch=19,
     cex=abs(resid(mod0))/max(resid(mod0))*2, 
     xlab="geographical xcoordinates",
     ylab="geographical y-coordinates")
# Basically, it doesn't appear there is much SAC!
# Plot per site:
bats_nona$rMod0 <- resid(mod0)
par(mfrow=c(6,6))
par(mar=c(0.5,0.5,0.5,0.5))
for(i in 1:nlevels(bats_nona$SITE)) {
  sitedat <- bats_nona[bats_nona$SITE==levels(bats_nona$SITE)[i],]
  sitecl <-  correlog(sitedat$MID_X, sitedat$MID_Y, sitedat$rMod0,
                      na.rm=T, increment=1, resamp=0)
  plot(sitecl$correlation[1:20], pch=16, type='o'); abline(h=0, lty='dashed')
}
# Nothing convincing per site.

### Check SAC in binomial GLMM with only site as RE:
correlog1a <- correlog(bats_nona$MID_X, bats_nona$MID_Y, residuals(mod1a),
                       na.rm=T, increment=1, resamp=0)
plot(correlog1a$correlation[1:20], pch=16, type='o'); abline(h=0, lty='dashed')
### Plot per site:
bats_nona$rMod1a <- resid(mod1a)
par(mfrow=c(6,6))
par(mar=c(0.5,0.5,0.5,0.5))
for(i in 1:nlevels(bats_nona$SITE)) {
  sitedat <- bats_nona[bats_nona$SITE==levels(bats_nona$SITE)[i],]
  sitecl <-  correlog(sitedat$MID_X, sitedat$MID_Y, sitedat$rMod1a,
                      na.rm=T, increment=1, resamp=0)
  plot(sitecl$correlation[1:20], pch=16, type='o'); abline(h=0, lty='dashed')
}
# Nothing convincing per site...


### So it still appears that the model w/o SAC structure and a nested RE is the best:
AIC(mod1, mod1a, mod0)
plot(mod1)
# "Error rate":
sum(bats_nona$OCC_PIPS!=as.numeric(fitted(mod1)>0.5))/nrow(bats_nona)
# R2?
r2mm(mod1)
# Binned residuals against fitted:
binnedplot(fitted(mod1), resid(mod1, type='pearson'))


### Try dropping interaction:
mod2 <- update(mod1, .~. -SECTION:TURB)
AIC(mod1, mod2); anova(mod1, mod2)
summary(mod2)


###
### COUNT MODELS (NEG. BIN.)
###
### Basic model: GLMM with nested RE:

mod10 <- glmmadmb(ALL_PIPS ~ SECTION*TURB + sMINTEMP +
                    sWINDS + 
                    sDAYNO + 
                    sTTMIDN + 
                    spBUILD + 
                    sD_TRE + 
                    (1|SITE/TRSCT), data=bats_nona, family="nbinom")
mod10a <- glmmadmb(ALL_PIPS ~ SECTION*TURB + sMINTEMP +
                    sWINDS + 
                    sDAYNO + 
                    sTTMIDN + 
                    spBUILD + 
                    sD_TRE + 
                    (1|SITE), data=bats_nona, family="nbinom")
dispZuur(mod10, mod10a)
AIC(mod10, mod10a)
# So the 'nested' model is quite a bit better.

# Drop interaction
mod11 <- update(mod10, .~. -SECTION:TURB)
AIC(mod10, mod10a, mod11)
anova(mod11, mod10)
summary(mod11)


###
### OCCURRENCE MODEL - GLMM NESTED - ALL PIPS - 
###  - QUAD TERMS WITH TEMP, WIND SPEED, TTMIDN, 
### MODEL SELECTION

full_mod <- formula(OCC_PIPS ~ SECTION*TURB +
                   sMINTEMP + I(sMINTEMP^2) +
                   sWINDS + I(sWINDS^2) + 
                   sTTMIDN + I(sTTMIDN^2) + 
                   sDAYNO +
                   spBUILD +
                   sD_TRE +
                   (1|SITE/TRSCT))
mod100 <- glmer(formula=full_mod, data=bats_nona, family=binomial, na.action='na.fail')
modset <- dredge(mod100, fixed="SECTION", 
                 subset=dc(sMINTEMP, `I(sMINTEMP^2)`) && 
                        dc(sWINDS, `I(sWINDS^2)`) &&
                        dc(sTTMIDN, `I(sTTMIDN^2)`), evaluate=F)
length(modset)
modset <- dredge(mod100, fixed="SECTION", 
                 subset=dc(sMINTEMP, `I(sMINTEMP^2)`) && 
                   dc(sWINDS, `I(sWINDS^2)`) &&
                   dc(sTTMIDN, `I(sTTMIDN^2)`), trace=T)

subset(modset, delta < 4)
subset(modset, cumsum(weight) <= .95)

save(modset, file='modset.Rdata')

### Example for complicated exclusion:
data(Cement)
fm <- lm(y ~ X1 + X2 + X3 + X4 + X5, Cement, na.action = na.fail)
dredge(fm, fixed="X1", subset=!(X2 && (X3 | X4 | X5)) &&  
                              !(X3 && (X2 | X4 | X5)) && 
                              !(X4 && (X2 | X3 | X5)) && 
                              !(X5 && (X2|X3|X4)))



mod1 <- glmer(OCC_PIPS ~ fSECTION*NOTURB + 
                sMINTEMP +
                sWINDS + 
                sDAYNO + 
                sTTMIDN + 
                spBUILD + 
                sD_TRE + 
                (1|SITE/TRSCT), data=bats_nona, family=binomial)
bats_nona$res1 <- resid(mod1, type='pearson')
# Diagnostic plots for basic binomial GLMM:
binnedplot(fitted(mod1), bats_nona$res1)
plot(as.vector(bats_nona$fSECTION), bats_nona$res1)
plot(bats_nona$NOTURB, bats_nona$res1)
plot(bats_nona$sMINTEMP, bats_nona$res1)
plot(bats_nona$sWINDS, bats_nona$res1)
plot(bats_nona$sDAYNO, bats_nona$res1)
plot(bats_nona$sTTMIDN, bats_nona$res1)
plot(bats_nona$spBUILD, bats_nona$res1)
plot(bats_nona$sD_TRE, bats_nona$res1)
plot(predict(mod1, type='response'),jitter(bats_nona$OCC_PIPS, 0.05))
dispZuur(mod1)
# 'Error rate':
sum(as.numeric(predict(mod1, type='response')>0.5)!=bats_nona$OCC_PIPS)/nrow(bats_nona)

### Same model as above but with SECTION as continuous:
mod2 <- glmer(OCC_PIPS ~ SECTION*NOTURB + 
                sMINTEMP +
                sWINDS + 
                sDAYNO + 
                sTTMIDN + 
                spBUILD + 
                sD_TRE + 
                (1|SITE/TRSCT), data=bats_nona, family=binomial)

### Same again but with TURB (single or multiple turbine?) as factor:
mod3 <- glmer(OCC_PIPS ~ SECTION*TURB + 
                sMINTEMP +
                sWINDS + 
                sDAYNO + 
                sTTMIDN + 
                spBUILD + 
                sD_TRE + 
                (1|SITE/TRSCT), data=bats_nona, family=binomial)
mod3a <- update(mod3, .~. -SECTION:TURB)
summary(mod3a)

AIC(mod1, mod2, mod3, mod3a)

### Now try 3a as a neg bin model:
mod4 <- glmmadmb(ALL_PIPS ~ SECTION*TURB + 
                   sMINTEMP +
                   sWINDS + 
                   sDAYNO + 
                   sTTMIDN + 
                   spBUILD + 
                   sD_TRE + 
                   (1|SITE/TRSCT), data=bats_nona, family="nbinom", verbose=T)
summary(mod4)

### How about a gamma distribution fit on no/ha?
bats_nona$ALL_PIPS_ha <- bats_nona$ALL_PIPS/bats_nona$AREA_ha
mod5 <- glmmadmb(ALL_PIPS_ha ~ SECTION*TURB + 
                   sMINTEMP +
                   sWINDS + 
                   sDAYNO + 
                   sTTMIDN + 
                   spBUILD + 
                   sD_TRE + 
                   (1|SITE/TRSCT), data=bats_nona, family="gamma", verbose=T)
### Total fail... massive failure in algorithm.

### So how about trying mod4 with an offset?
mod4a <- glmmadmb(ALL_PIPS ~ SECTION*TURB + 
                   sMINTEMP +
                   sWINDS + 
                   sDAYNO + 
                   sTTMIDN + 
                   spBUILD + 
                   sD_TRE + 
                   (1|SITE/TRSCT) + offset(log(AREA_ha)), data=bats_nona, family="nbinom", verbose=T)
### Drop the interaction:
mod4a1 <- update(mod4a, .~. -SECTION:TURB)
AIC(mod4a1, mod4a)



# Need quad term for wind?
mod1a <- update(mod1, .~. +sWINDS2)
anova(mod1, mod1a) # No.

# Quad term for TTMIDN?
mod1a <- update(mod1, .~. +sTTMIDN2)
anova(mod1, mod1a) # Yes!

# Quad term for DAYNO?
mod1b <- update(mod1a, .~. +sDAYNO2)
anova(mod1a, mod1b) # Yes?

# New 'full' model:
mod2 <- mod1b
bats_nona$res2 <- resid(mod2, type='pearson')
binnedplot(fitted(mod2), bats_nona$res2)
plot(as.vector(bats_nona$fSECTION), bats_nona$res2)
plot(bats_nona$NOTURB, bats_nona$res2)
plot(bats_nona$sMINTEMP, bats_nona$res2)
plot(bats_nona$sWINDS, bats_nona$res2)
plot(bats_nona$sDAYNO, bats_nona$res2)
plot(bats_nona$sTTMIDN, bats_nona$res2)
plot(bats_nona$spBUILD, bats_nona$res2)
plot(bats_nona$sD_TRE, bats_nona$res2)
plot(predict(mod2, type='response'),jitter(bats_nona$OCC_PIPS, 0.05))
dispZuur(mod2)
sum(as.numeric(predict(mod2, type='response')>0.5)!=bats_nona$OCC_PIPS)/nrow(bats_nona)
anova(mod1, mod2)

# Do we need SECTION*NOTURBS:
mod2a <- update(mod2, .~. -fSECTION:NOTURB)
anova(mod2a, mod2) # Borderline! Keep for now...

# Do we need sMINTEMP?
mod2a <- update(mod2, .~. -sMINTEMP)
anova(mod2a, mod2) # No, try dropping this.

mod3 <- mod2a
summary(mod3)
# Diagnostics:
bats_nona$res3 <- resid(mod3, type='pearson')
binnedplot(fitted(mod3), bats_nona$res3)
plot(as.vector(bats_nona$fSECTION), bats_nona$res3)
plot(bats_nona$NOTURB, bats_nona$res3)
plot(bats_nona$sWINDS, bats_nona$res3)
plot(bats_nona$sDAYNO, bats_nona$res3)
plot(bats_nona$sTTMIDN, bats_nona$res3)
plot(bats_nona$spBUILD, bats_nona$res3)
plot(bats_nona$sD_TRE, bats_nona$res3)
plot(predict(mod3, type='response'),jitter(bats_nona$OCC_PIPS, 0.05))
dispZuur(mod3)
sum(as.numeric(predict(mod3, type='response')>0.5)!=bats_nona$OCC_PIPS)/nrow(bats_nona)
round(summary(mod3)$coef,3)
r2mm(mod3)

# Point predictions:
p_s1 <- data.frame(fSECTION=factor("1", levels=levels(bats_nona$fSECTION)),
                   NOTURB=1,
                   sWINDS=median(bats_nona$sWINDS),
                   sDAYNO=median(bats_nona$sDAYNO),
                   sTTMIDN=median(bats_nona$sTTMIDN),
                   spBUILD=median(bats_nona$spBUILD),
                   sD_TRE=median(bats_nona$sD_TRE),
                   sTTMIDN2=median(bats_nona$sTTMIDN2),
                   sDAYNO2=median(bats_nona$sDAYNO2)
                   )
p_s2 <- p_s1; p_s2$fSECTION <- factor("2", levels=levels(bats_nona$fSECTION))
p_s3 <- p_s1; p_s3$fSECTION <- factor("3", levels=levels(bats_nona$fSECTION))
p_s4 <- p_s1; p_s4$fSECTION <- factor("4", levels=levels(bats_nona$fSECTION))
p_s5 <- p_s1; p_s5$fSECTION <- factor("5", levels=levels(bats_nona$fSECTION))

predict(mod3, newdata=p_s1, type='response', re.form=NA)
predict(mod3, newdata=p_s2, type='response', re.form=NA)
predict(mod3, newdata=p_s3, type='response', re.form=NA)
predict(mod3, newdata=p_s4, type='response', re.form=NA)
predict(mod3, newdata=p_s5, type='response', re.form=NA)

# Replicate the above manually:
mp_s1 <- c(1,0,0,0,0,
           2,
           median(bats_nona$sWINDS),
           median(bats_nona$sDAYNO),
           median(bats_nona$sTTMIDN),
           median(bats_nona$spBUILD),
           median(bats_nona$sD_TRE),
           median(bats_nona$sTTMIDN2),
           median(bats_nona$sDAYNO2),
           2,0,0,0
           )
mp_s2 <- mp_s1
mp_s2[2] <- 1; mp_s2[14] <- 1;
mp_s3 <- mp_s1
mp_s3[3] <- 1; mp_s3[15] <- 1;
mp_s4 <- mp_s1
mp_s4[4] <- 1; mp_s4[16] <- 1;
mp_s5 <- mp_s1
mp_s5[5] <- 1; mp_s5[17] <- 1;

plogis(fixef(mod3) %*% mp_s1)
plogis(fixef(mod3) %*% mp_s2)
plogis(fixef(mod3) %*% mp_s3)
plogis(fixef(mod3) %*% mp_s4)
plogis(fixef(mod3) %*% mp_s5)
# Good, manually works.

# Simulate coefficients:
mod3_sim <- sim(mod3, n.sim=1000)
mod3_simfix <- attr(mod3_sim, 'fixef')

# Predicted means and lower/upper quantiles:
# per section:
predict_s1 <- apply(mod3_simfix, 1, function(x) plogis(x%*%mp_s1))
predict_s2 <- apply(mod3_simfix, 1, function(x) plogis(x%*%mp_s2))
predict_s3 <- apply(mod3_simfix, 1, function(x) plogis(x%*%mp_s3))
predict_s4 <- apply(mod3_simfix, 1, function(x) plogis(x%*%mp_s4))
predict_s5 <- apply(mod3_simfix, 1, function(x) plogis(x%*%mp_s5))
predict_all <- data.frame(section=factor(rep(c("1","2","3","4","5"),each=1000)),
                          pred_y=c(predict_s1, predict_s2, predict_s3, predict_s4, predict_s5))
boxplot(predict_all$pred_y ~ predict_all$section)




# Point predictions w/o interaction between no. turbines/section AND w/o main effect of no. turbines:
mod4 <- update(mod3, .~. -fSECTION:NOTURB)
anova(mod4, mod3)
mod4a <- update(mod4, .~. -NOTURB)
anova(mod4a, mod4)
mod4 <- mod4a
p_s1a <- data.frame(fSECTION=factor("1", levels=levels(bats_nona$fSECTION)),
                   sWINDS=median(bats_nona$sWINDS),
                   sDAYNO=median(bats_nona$sDAYNO),
                   sTTMIDN=median(bats_nona$sTTMIDN),
                   spBUILD=median(bats_nona$spBUILD),
                   sD_TRE=median(bats_nona$sD_TRE),
                   sTTMIDN2=median(bats_nona$sTTMIDN2),
                   sDAYNO2=median(bats_nona$sDAYNO2)
)
p_s2a <- p_s1a; p_s2a$fSECTION <- factor("2", levels=levels(bats_nona$fSECTION))
p_s3a <- p_s1a; p_s3a$fSECTION <- factor("3", levels=levels(bats_nona$fSECTION))
p_s4a <- p_s1a; p_s4a$fSECTION <- factor("4", levels=levels(bats_nona$fSECTION))
p_s5a <- p_s1a; p_s5a$fSECTION <- factor("5", levels=levels(bats_nona$fSECTION))
predict(mod4, newdata=p_s1a, type='response', re.form=NA)
predict(mod4, newdata=p_s2a, type='response', re.form=NA)
predict(mod4, newdata=p_s3a, type='response', re.form=NA)
predict(mod4, newdata=p_s4a, type='response', re.form=NA)
predict(mod4, newdata=p_s5a, type='response', re.form=NA)
