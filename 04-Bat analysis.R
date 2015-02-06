###
###
### ANALYSE BAT DATA 2013+2014
### Avoid all edits to data file, please do in previous script.

library(lme4)
library(arm)
library(glmmADMB)
library(gstat)
library(sp)
library(MuMIn)
library(snow)
library(corrplot)
library(plyr)

source('../../../000_R/RSimExamples/helpjm.r')

### Load bat master data:
bats <- read.csv('data/BAT DATA 2013 and 2014 MASTER.csv',header=T)

### Check factors, etc:
is.factor(bats$SITE)
is.factor(bats$TRSCT)
bats$fSECTION <- factor(bats$fSECTION)
is.factor(bats$fSECTION)
is.factor(bats$NOTURB)

### Check transect section inventory - only for housekeeping/checking purposes:
temp <- ddply(bats, .(SITE, TRSCT, fSECTION), summarise, noDates=length(unique(DATE)), noObs=length(DATE))
write.csv(temp, 'temp.csv', row.names=F)

### Proportion of zero's per site:
tapply(bats$ALL_PIPS, bats$SITE, function(x) sum(x==0)/length(x))
### Average proportion of zero's per site:
mean(tapply(bats$ALL_PIPS, bats$SITE, function(x) sum(x==0)/length(x)))

### First date:
min(as.Date(as.vector(bats$DATE),format='%Y-%m-%d'))

### Last date:
max(as.Date(as.vector(bats$DATE),format='%Y-%m-%d'))

### Binary responses and some housekeeping:
bats$OCC_PIPS <- as.numeric(as.logical(bats$ALL_PIPS))
bats$OCC_OTHR <- as.numeric(as.logical(bats$ALL_OTHER))
bats$MINTEMP <- bats$min_air_temp

### Create TURB indicator var for "single" or "multiple" (>1) turbines:
bats$TURB <- "single"
bats$TURB[bats$NOTURB>1] <- "multiple"
bats$TURB <- factor(bats$TURB)
bats$TURB <- relevel(bats$TURB, ref="single")

bats_nona <- bats[!is.na(bats$WINDS),]
bats_nona <- droplevels(bats_nona)


###
###
### 1. SELECT HABITAT DESCRIPTORS USING SEPARATE MODEL SELECTION.
### MODEL FOR PIP OCCURRENCE WITH ONLY HAB VARIABLES:

# Tests if I can use cloglog link appropriately:
# mod_hab1.1 <- glmer(OCC_PIPS ~ 1 + (1|SITE/TRSCT), data=bats_nona, family='binomial'(link='logit'))
# mod_hab1.2 <- glmer(OCC_PIPS ~ 1 + (1|SITE/TRSCT), data=bats_nona, family='binomial'(link='cloglog'))
# plogis(fixef(mod_hab1.1))
# make.link('logit')$linkinv(fixef(mod_hab1.1))
# make.link('cloglog')$linkinv(fixef(mod_hab1.2))
# # Now with offset:
# mod_hab1.3 <- glmer(OCC_PIPS ~ 1 + (1|SITE/TRSCT), offset=log(AREA_ha), 
#                     data=bats_nona, family='binomial'(link='cloglog'))
# make.link('cloglog')$linkinv(fixef(mod_hab1.3))

# First fit non-standardised models:
mod_hab1 <- glmer(OCC_PIPS ~ pBUILD + pTREE + pRDTRK + pROADS + pROADS + pRGRAS + 
                    EDGED + D_LIN + D_BUI + D_WAT + D_TRE + (1|SITE/TRSCT), offset=log(AREA_ha), 
                  data=bats_nona, family='binomial'(link='cloglog'), na.action='na.fail')
mod_habz1 <- standardize(mod_hab1)

modset_habz1 <- dredge(mod_habz1, 
                   subset=
                     !(z.pBUILD && (z.pTREE | z.pRDTRK | z.pROADS | z.pRGRAS )) &&
                     !(z.pTREE && (z.pBUILD | z.pRDTRK | z.pROADS | z.pRGRAS )) &&
                     !(z.pRDTRK && (z.pTREE  | z.pBUILD | z.pROADS | z.pRGRAS )) &&
                     !(z.pROADS && (z.pTREE | z.pRDTRK | z.pBUILD | z.pRGRAS )) &&
                     !(z.pRGRAS && (z.pTREE | z.pRDTRK | z.pROADS | z.pBUILD )) &&
                     !(z.pBUILD && z.D_BUI) &&
                     !(z.pTREE && z.D_TRE) &&
                     !(z.pRDTRK && z.EDGED) &&
                     !(z.D_LIN && z.EDGED)
                   , evaluate=F)
length(modset_habz1)

# system.time({
#   modset_habz1 <- dredge(mod_habz1, 
#                          subset=
#                            !(z.pBUILD && (z.pTREE | z.pRDTRK | z.pROADS | z.pRGRAS )) &&
#                            !(z.pTREE && (z.pBUILD | z.pRDTRK | z.pROADS | z.pRGRAS )) &&
#                            !(z.pRDTRK && (z.pTREE  | z.pBUILD | z.pROADS | z.pRGRAS )) &&
#                            !(z.pROADS && (z.pTREE | z.pRDTRK | z.pBUILD | z.pRGRAS )) &&
#                            !(z.pRGRAS && (z.pTREE | z.pRDTRK | z.pROADS | z.pBUILD )) &&
#                            !(z.pBUILD && z.D_BUI) &&
#                            !(z.pTREE && z.D_TRE) &&
#                            !(z.pRDTRK && z.pROADS) &&
#                            !(z.pRDTRK && z.EDGED) &&
#                            !(z.D_LIN && z.EDGED)
#                          , trace=T)
# })
# user  system elapsed 
# 153.543   0.551 195.446 
# save(modset_habz1, file='modset_habz1.RData')
# subset(modset_habz1, delta<4)

# Repeat above on all clusters:
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 8), type = clusterType))
clusterExport(clust, "bats_nona")
clusterExport(clust, "glmer")
clusterExport(clust, "fixef")
clusterExport(clust, "vcov")
clusterExport(clust, "forceSymmetric")

system.time({
  modset_habz1 <- pdredge(mod_habz1, cluster=clust,
                         subset=
                           !(z.pBUILD && (z.pTREE | z.pRDTRK | z.pROADS | z.pRGRAS )) &&
                           !(z.pTREE && (z.pBUILD | z.pRDTRK | z.pROADS | z.pRGRAS )) &&
                           !(z.pRDTRK && (z.pTREE  | z.pBUILD | z.pROADS | z.pRGRAS )) &&
                           !(z.pROADS && (z.pTREE | z.pRDTRK | z.pBUILD | z.pRGRAS )) &&
                           !(z.pRGRAS && (z.pTREE | z.pRDTRK | z.pROADS | z.pBUILD )) &&
                           !(z.pBUILD && z.D_BUI) &&
                           !(z.pTREE && z.D_TRE) &&
                           !(z.pRDTRK && z.pROADS) &&
                           !(z.pRDTRK && z.EDGED) &&
                           !(z.D_LIN && z.EDGED)
                         , trace=T)
})
# user  system elapsed 
# 1.581   0.673  43.787 Nils

# save(modset_habz1, file='modset_habz1.Rdata')
load('modset_habz1.Rdata')
subset(modset_habz1, delta<4)
subset(modset_habz1, delta<10)

# So on the basis of AIC selection, let's consider D_BUI, D_WAT, EDGED and pTREE as our hab variables.

# # Second AIC based selection, exclude edge density:
# mod_hab2 <- glmer(OCC_PIPS ~ pBUILD + pTREE + pRDTRK + pROADS + pROADS + pRGRAS + 
#                     D_LIN + D_BUI + D_WAT + D_TRE + (1|SITE/TRSCT), offset=log(AREA_ha), 
#                   data=bats_nona, family='binomial'(link='cloglog'), na.action='na.fail')
# mod_habz2 <- standardize(mod_hab2)
# 
# modset_habz2 <- dredge(mod_habz2, 
#                        subset=
#                          !(z.pBUILD && (z.pTREE | z.pRDTRK | z.pROADS | z.pRGRAS )) &&
#                          !(z.pTREE && (z.pBUILD | z.pRDTRK | z.pROADS | z.pRGRAS )) &&
#                          !(z.pRDTRK && (z.pTREE  | z.pBUILD | z.pROADS | z.pRGRAS )) &&
#                          !(z.pROADS && (z.pTREE | z.pRDTRK | z.pBUILD | z.pRGRAS )) &&
#                          !(z.pRGRAS && (z.pTREE | z.pRDTRK | z.pROADS | z.pBUILD )) &&
#                          !(z.pBUILD && z.D_BUI) &&
#                          !(z.pTREE && z.D_TRE) &&
#                          !(z.pRDTRK && z.pROADS)
#                        , evaluate=F)
# length(modset_habz2)
# modset_habz2 <- pdredge(mod_habz2, cluster=clust,
#                        subset=
#                          !(z.pBUILD && (z.pTREE | z.pRDTRK | z.pROADS | z.pRGRAS )) &&
#                          !(z.pTREE && (z.pBUILD | z.pRDTRK | z.pROADS | z.pRGRAS )) &&
#                          !(z.pRDTRK && (z.pTREE  | z.pBUILD | z.pROADS | z.pRGRAS )) &&
#                          !(z.pROADS && (z.pTREE | z.pRDTRK | z.pBUILD | z.pRGRAS )) &&
#                          !(z.pRGRAS && (z.pTREE | z.pRDTRK | z.pROADS | z.pBUILD )) &&
#                          !(z.pBUILD && z.D_BUI) &&
#                          !(z.pTREE && z.D_TRE) &&
#                          !(z.pRDTRK && z.pROADS)
#                        , trace=T)
# save(modset_habz2, file='modset_habz2.Rdata')
# stopCluster(clust)
# subset(modset_habz2, delta<4)

# More hab variables retained - use original hab set (D_BUI, D_WAT, EDGED and pTREE)


# Check BIC selection:
# clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
# clust <- try(makeCluster(getOption("cl.cores", 8), type = clusterType))
# clusterExport(clust, "bats_nona")
# clusterExport(clust, "glmer")
# clusterExport(clust, "fixef")
# system.time({
#   modset_habz1_BIC <- pdredge(mod_habz1, cluster=clust,
#                           subset=
#                             !(z.pBUILD && (z.pTREE | z.pRDTRK | z.pROADS | z.pRGRAS )) &&
#                             !(z.pTREE && (z.pBUILD | z.pRDTRK | z.pROADS | z.pRGRAS )) &&
#                             !(z.pRDTRK && (z.pTREE  | z.pBUILD | z.pROADS | z.pRGRAS )) &&
#                             !(z.pROADS && (z.pTREE | z.pRDTRK | z.pBUILD | z.pRGRAS )) &&
#                             !(z.pRGRAS && (z.pTREE | z.pRDTRK | z.pROADS | z.pBUILD )) &&
#                             !(z.pBUILD && z.D_BUI) &&
#                             !(z.pTREE && z.D_TRE) &&
#                             !(z.pRDTRK && z.pROADS) &&
#                             !(z.pRDTRK && z.EDGED) &&
#                             !(z.D_LIN && z.EDGED)
#                           , rank='BIC', trace=T)
# })
# # user  system elapsed 
# # 1.651   0.700  43.930 
# stopCluster(clust)
# save(modset_habz1_BIC, file='modset_habz1_BIC.Rdata')
# subset(modset_habz1_BIC, delta<4)
# subset(modset_habz1_BIC, delta<10)
# 
# # Only EDGED and pTREE following BIC selection.
# 
# # Check DIC selection:
# clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
# clust <- try(makeCluster(getOption("cl.cores", 8), type = clusterType))
# clusterExport(clust, "bats_nona")
# clusterExport(clust, "glmer")
# clusterExport(clust, "fixef")
# system.time({
#   modset_habz1_DIC <- pdredge(mod_habz1, cluster=clust,
#                               subset=
#                                 !(z.pBUILD && (z.pTREE | z.pRDTRK | z.pROADS | z.pRGRAS )) &&
#                                 !(z.pTREE && (z.pBUILD | z.pRDTRK | z.pROADS | z.pRGRAS )) &&
#                                 !(z.pRDTRK && (z.pTREE  | z.pBUILD | z.pROADS | z.pRGRAS )) &&
#                                 !(z.pROADS && (z.pTREE | z.pRDTRK | z.pBUILD | z.pRGRAS )) &&
#                                 !(z.pRGRAS && (z.pTREE | z.pRDTRK | z.pROADS | z.pBUILD )) &&
#                                 !(z.pBUILD && z.D_BUI) &&
#                                 !(z.pTREE && z.D_TRE) &&
#                                 !(z.pRDTRK && z.pROADS) &&
#                                 !(z.pRDTRK && z.EDGED) &&
#                                 !(z.D_LIN && z.EDGED)
#                               , rank='DIC', trace=T)
# })
# # user  system elapsed 
# # 1.635   0.682  44.355
# stopCluster(clust)
# save(modset_habz1_DIC, file='modset_habz1_DIC.Rdata')
# subset(modset_habz1_DIC, delta<4)
# subset(modset_habz1_DIC, delta<10)
# 
# # For DIC, use D_BUI, D_TRE, D_WAT, EDGED, pROADS, and pTREE.

###
### OCCURRENCE MODELS (PROBABILITY OF A PASS)
###
### Basic model: GLMM with nested RE. Habitat vars based on AIC selection, and use AIC selection.
### (so that's D_BUI, D_WAT, EDGED and pTREE)

# First fit unstandardised model:
m1 <- glmer(OCC_PIPS  ~ fSECTION*TURB + 
                        MINTEMP + 
                        DAYNO +
                        TTMIDN +
                        I(TTMIDN^2) + 
                        WINDS +
                        D_BUI + 
                        D_WAT +
                        EDGED +
                        pTREE + 
                        (1|SITE/TRSCT), offset=log(AREA_ha), 
            data=bats_nona, family='binomial'(link='cloglog'), na.action='na.fail')
# Standardise predictors:
m1z <- standardize(m1)

m1z_set1 <- dredge(m1z, subset=dc(z.TTMIDN, `I(z.TTMIDN^2)`), evaluate=F)
length(m1z_set1)

# Set up cluster:
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 8), type = clusterType))
clusterExport(clust, "bats_nona")
clusterExport(clust, "glmer")
clusterExport(clust, "fixef")

# system.time({
#   m1z_set1 <- pdredge(m1z, cluster=clust, 
#                       subset=dc(z.TTMIDN, `I(z.TTMIDN^2)`), trace=T)  
# })
# user   system  elapsed 
# 12.258    7.935 2497.640 
# save(m1z_set1, file='m1z_set1.Rdata')
# stopCluster(clust)
load('m1z_set1.Rdata')
subset(m1z_set1, delta<4)
subset(m1z_set1, delta<2)
# m1z_av <- model.avg(m1z_set1, delta<4, fit=T)
# save(m1z_av, file='m1z_av.Rdata')
load('m1z_av.Rdata')
summary(m1z_av)

### "Simpler habitat" set: GLMM with nested RE. Only those hab vars that are present in all four d<4 hab models.
### (so that's EDGED and pTREE)

m2 <- glmer(OCC_PIPS  ~ fSECTION*TURB + 
              MINTEMP + 
              DAYNO +
              TTMIDN +
              I(TTMIDN^2) + 
              WINDS +
              EDGED +
              pTREE + 
              (1|SITE/TRSCT), offset=log(AREA_ha), 
            data=bats_nona, family='binomial'(link='cloglog'), na.action='na.fail')
# Standardise predictors:
m2z <- standardize(m2)

m2z_set1 <- dredge(m2z, subset=dc(z.TTMIDN, `I(z.TTMIDN^2)`), evaluate=F)
length(m2z_set1)

# Set up cluster:
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 8), type = clusterType))
clusterExport(clust, "bats_nona")
clusterExport(clust, "glmer")
clusterExport(clust, "fixef")

system.time({
  m2z_set1 <- pdredge(m2z, cluster=clust, 
                      subset=dc(z.TTMIDN, `I(z.TTMIDN^2)`), trace=T)  
})
save(m2z_set1, file='m2z_set1.Rdata')
subset(m2z_set1, delta<4)
load('m2z_set1.Rdata')

m2z_av <- model.avg(m2z_set1, delta<4, fit=T, trace=T)
save(m2z_av, file='m2z_av.Rdata')
load('m2z_av.Rdata')
#vcov(m2z_av)
summary(m2z_av)
summary(m2z_av)$avg.model



###
### PLOT PREDICTIONS: prediction interval by sim():
###

linkinv <- m2z@resp$family$linkinv

# Standardised values for single and multiple turb
c.turb <- unique(m2z@frame$c.TURB)[order(unique(m2z@frame$c.TURB))]

p1 <- c(1,0,0,0,0,   # fSECTION 1
        0,           # z.DAYNO
        0,           # z.EDGED
        0,           # z.TTMIDN,
        0,           # I(z.TTMIDN^2)
        0,           # z.WINDS
        0,           # z.pTREE
        c.turb[1],   # SINGLE TURB
        0,           # z.MINTEMP
        0,           # TURB*section 2
        0,           # TURB*section 3
        0,           # TURB*section 4
        0            # TURB*section 5
        )

linkinv(p1 %*% summary(m2z_av)$avg.model[,1])

p_single <- cbind(rep(1,5),          # Intercept (fSECTION 1)
                  c(0,1,0,0,0),      # fSECTION 2
                  c(0,0,1,0,0),      # fSECTION 3
                  c(0,0,0,1,0),      # fSECTION 4
                  c(0,0,0,0,1),      # fSECTION 5
                  rep(0,5),          # z.DAYNO
                  rep(0,5),          # z.EDGED
                  rep(0,5),          # z.TTMIDN,
                  rep(0,5),          # I(z.TTMIDN^2),
                  rep(0,5),          # z.WINDS
                  rep(0,5),          # z.pTREE
                  rep(c.turb[1],5),  # SINGLE TURB = -0.5639652
                  rep(0,5),           # z.MINTEMP                  
                  cbind(c(0,1,0,0,0),# Interaction, fsection 2
                        c(0,0,1,0,0),# Interaction, fsection 3
                        c(0,0,0,1,0),# Interaction, fsection 4
                        c(0,0,0,0,1))*c.turb[1] # fsection 5 - SINGLE TURB
                  )

p_multp <- cbind(rep(1,5),          # Intercept (fSECTION 1)
                  c(0,1,0,0,0),      # fSECTION 2
                  c(0,0,1,0,0),      # fSECTION 3
                  c(0,0,0,1,0),      # fSECTION 4
                  c(0,0,0,0,1),      # fSECTION 5
                  rep(0,5),          # z.DAYNO
                  rep(0,5),          # z.EDGED
                  rep(0,5),          # z.TTMIDN,
                  rep(0,5),          # I(z.TTMIDN^2),
                  rep(0,5),          # z.WINDS
                  rep(0,5),          # z.pTREE
                  rep(0,5),          # z.MINTEMP
                  rep(c.turb[2],5),  # MULT TURB
                  cbind(c(0,1,0,0,0),# Interaction, fsection 2
                        c(0,0,1,0,0),# Interaction, fsection 3
                        c(0,0,0,1,0),# Interaction, fsection 4
                        c(0,0,0,0,1))*c.turb[2] # fsection 5 - MULT TURB
)

# Point predictions, SINGLE turbine:
pts1_single <- linkinv(p_single %*% summary(m2z_av)$avg.model[,1])

# Point predictions, MULTIPLE turbines:
pts1_mult <- linkinv(p_multp %*% summary(m2z_av)$avg.model[,1])

# Now simulate from posteriors of betas to get prediction intervals:

# Custom sim() function to work with averaging object:
sim_man <- function(avmod, S) {
  sig.hat <- 1            # FIX RESIDUAL SE TO 1 (binomial model)
  V.beta <- vcov(avmod)   # Model estimated variance-covariance matrix
  n <- nrow(as.data.frame(model.matrix(avmod)))
  k <- nrow(summary(m2z_av)$avg.model) + length(ranef(m2z))
  sigma <- as.vector(NULL)
  beta <- as.data.frame(NULL)
  for (s in 1:S) {
    sigma <- c(sigma, sig.hat*sqrt((n-k)/rchisq(1,n-k)))
    beta <- rbind(beta, mvrnorm(1, summary(m2z_av)$avg.model[,1], V.beta*sigma[s]^2))
  }
  names(beta) <- names(summary(m2z_av)$avg.model[,1])
  return(list(fixef=beta, sigma=sigma))
}

m2z_av_sim <- sim_man(m2z_av, 1000)

# MEAN prediction, SINGLE turbine:
pts2_single <- linkinv(apply(apply(m2z_av_sim$fixef, 1, function(x) p_single %*% x ), 
                             1, 
                             mean))

# Lower and upper quantiles for point predictions, SINGLE turbine:
qtls_single <- linkinv(apply(apply(m2z_av_sim$fixef, 1, function(x) p_single %*% x ), 
                               1, 
                               function(x) quantile(x, probs=c(0.025, 0.975))))

# PREDICTIONS FOR SINGLE TURBINE:
preds_single <- data.frame(pts1_single, pts2_single, t(qtls_single))


# MEAN prediction, MULT turbine:
pts2_mult <- linkinv(apply(apply(m2z_av_sim$fixef, 1, function(x) p_multp %*% x ), 
                           1, 
                           mean))

# Lower and upper quantiles for point predictions, SINGLE turbine:
qtls_mult <- linkinv(apply(apply(m2z_av_sim$fixef, 1, function(x) p_multp %*% x ), 
                             1, 
                             function(x) quantile(x, probs=c(0.025, 0.975))))

# PREDICTIONS FOR MULTIPLE TURBINE:
preds_mult <- data.frame(pts1_mult, pts2_mult, t(qtls_mult))

# Calculate observered mean proportions 0/1 for single and multiple turbines:
bats_nona_single <- bats_nona[bats_nona$TURB=='single',]
bats_nona_multp <- bats_nona[bats_nona$TURB=='multiple',]
bats_nona_single <- droplevels(bats_nona_single); bats_nona_multp <- droplevels(bats_nona_multp)

obs_single <- tapply(bats_nona_single$OCC_PIPS, bats_nona_single$fSECTION, sum)/
  tapply(bats_nona_single$OCC_PIPS, bats_nona_single$fSECTION, length)

obs_multp <- tapply(bats_nona_multp$OCC_PIPS, bats_nona_multp$fSECTION, sum)/
  tapply(bats_nona_multp$OCC_PIPS, bats_nona_multp$fSECTION, length)


### PREDICTION intervals as per Bolker:

p_single_se <- sqrt(diag(p_single %*% vcov(m2z_av) %*% t(p_single)))
p_mult_se <- sqrt(diag(p_multp %*% vcov(m2z_av) %*% t(p_multp)))

preds_single$lwr <- linkinv(p_single %*% summary(m2z_av)$avg.model[,1] - 1.96*p_single_se)
preds_single$upr <- linkinv(p_single %*% summary(m2z_av)$avg.model[,1] + 1.96*p_single_se)

preds_mult$lwr <- linkinv(p_multp %*% summary(m2z_av)$avg.model[,1] - 1.96*p_mult_se)
preds_mult$upr <- linkinv(p_multp %*% summary(m2z_av)$avg.model[,1] + 1.96*p_mult_se)


# PLOT:

par(mfrow=c(2,2))
bars_single <- barplot(preds_single$pts2_single, 
                       ylim=c(0, max(c(preds_single$X97.5.,preds_mult$X97.5.))),
                       names.arg=c('0-100','100-200','200-300','300-400','400-500'),
                       main='Single turbine',
                       xlab='Distance band (m)',
                       ylab='Probability of a bat pass / ha'
)
arrows(bars_single, preds_single$X2.5., 
       bars_single, preds_single$X97.5., 
       code=3, length=0.1, angle=90)
# Add Bolker prediction interval:
# arrows(bars_single+1/50, preds_single$lwr, 
#        bars_single+1/50, preds_single$upr, 
#        code=0, length=0.1, angle=90, col='red',lty='dashed')

points(bars_single, obs_single, cex=1.5, pch=16, col='black')

bars_mult <- barplot(preds_mult$pts2_mult, 
                     ylim=c(0, max(c(preds_single$X97.5.,preds_mult$X97.5.))),
                     names.arg=c('0-100','100-200','200-300','300-400','400-500'),
                     main='Multiple turbines',
                     xlab='Distance band (m)',
                     ylab='Probability of a bat pass / ha'
)
arrows(bars_mult, preds_mult$X2.5., 
       bars_mult, preds_mult$X97.5., 
       code=3, length=0.1, angle=90)
# Add Bolker prediction interval:
# arrows(bars_mult+1/50, preds_mult$lwr, 
#        bars_mult+1/50, preds_mult$upr, 
#        code=0, length=0.1, angle=90, col='red',lty='dashed')

points(bars_mult, obs_multp, cex=1.5, pch=16, col='black')



###
### Same prediction but with 'full' model:
### 
p_single_full <- cbind(rep(1,5),          # Intercept (fSECTION 1)
                       c(0,1,0,0,0),      # fSECTION 2
                       c(0,0,1,0,0),      # fSECTION 3
                       c(0,0,0,1,0),      # fSECTION 4
                       c(0,0,0,0,1),      # fSECTION 5
                       rep(-0.5639652,5), # SINGLE TURB = -0.5639652
                       rep(0,5),           # z.MINTEMP                
                       rep(0,5),          # z.DAYNO
                       rep(0,5),          # z.TTMIDN,
                       rep(0,5),          # I(z.TTMIDN^2),
                       rep(0,5),          # z.WINDS                       
                       rep(0,5),          # z.EDGED
                       rep(0,5),          # z.pTREE
                       cbind(c(0,1,0,0,0),# Interaction, fsection 2
                             c(0,0,1,0,0),# Interaction, fsection 3
                             c(0,0,0,1,0),# Interaction, fsection 4
                             c(0,0,0,0,1))*-0.5639652 # fsection 5 - SINGLE TURB
                      )
p_mult_full   <- cbind(rep(1,5),          # Intercept (fSECTION 1)
                       c(0,1,0,0,0),      # fSECTION 2
                       c(0,0,1,0,0),      # fSECTION 3
                       c(0,0,0,1,0),      # fSECTION 4
                       c(0,0,0,0,1),      # fSECTION 5
                       rep(0.4360348,5), # SINGLE TURB = 0.4360348
                       rep(0,5),           # z.MINTEMP                
                       rep(0,5),          # z.DAYNO
                       rep(0,5),          # z.TTMIDN,
                       rep(0,5),          # I(z.TTMIDN^2),
                       rep(0,5),          # z.WINDS                       
                       rep(0,5),          # z.EDGED
                       rep(0,5),          # z.pTREE
                       cbind(c(0,1,0,0,0),# Interaction, fsection 2
                             c(0,0,1,0,0),# Interaction, fsection 3
                             c(0,0,0,1,0),# Interaction, fsection 4
                             c(0,0,0,0,1))*0.4360348 # fsection 5 - SINGLE TURB
)

m2z_sim <- sim(m2z, 1000)@fixef
temp <- apply(m2z_sim, 1, function(x) linkinv(p_single_full %*% x))
pts_single_full <- data.frame(mean=apply(temp, 1, mean))
pts_single_full <- cbind(pts_single_full, t(apply(temp, 1, function(x) quantile(x, probs=c(0.025, 0.975)))))
temp <- apply(m2z_sim, 1, function(x) linkinv(p_mult_full %*% x))
pts_mult_full <- data.frame(mean=apply(temp, 1, mean))
pts_mult_full <- cbind(pts_mult_full, t(apply(temp, 1, function(x) quantile(x, probs=c(0.025, 0.975)))))

bars_single <- barplot(pts_single_full$mean, 
                       ylim=c(0, max(c(pts_single_full[,3], pts_mult_full[,3]))),
                       names.arg=c('0-100','100-200','200-300','300-400','400-500'),
                       main='Single turbine',
                       xlab='Distance band (m)',
                       ylab='Probability of a bat pass / ha'
)
arrows(bars_single, pts_single_full[,2], 
       bars_single, pts_single_full[,3], 
       code=3, length=0.1, angle=90)

bars_mult <- barplot(pts_mult_full$mean, 
                       ylim=c(0, max(c(pts_single_full[,3], pts_mult_full[,3]))),
                       names.arg=c('0-100','100-200','200-300','300-400','400-500'),
                       main='Multiple turbines',
                       xlab='Distance band (m)',
                       ylab='Probability of a bat pass / ha'
)
arrows(bars_mult, pts_mult_full[,2], 
       bars_mult, pts_mult_full[,3], 
       code=3, length=0.1, angle=90)




###
### 'COUNT' MODEL USING FULL STRUCTURE FROM BINOMIAL MODEL ABOVE
### Also use standardised predictors from binomial model.






######
###### WEIRD ATTEMPT AT TRYING TO INCORPORATE RE UNCERTAINTY IN PREDICTIONS:
sim_ALL <- function(mod) {
  out <- as.data.frame(NULL)
  for(i in 1:50) {
    simdat <- attr(mod, 'frame')
    simdat$AREA_ha <- exp(m2z@resp$offset)
    simdat$OCC_PIPS <- as.vector(unlist(simulate(mod)))
    sim_fit <- update(mod, .~., data=simdat)
    out <- rbind(out, 
      c(
        predict(sim_fit, type='response', re.form=NA, newdata=data.frame(
          fSECTION=as.factor(1:5),
          c.TURB=rep(-0.5639652, 5),
          z.MINTEMP=rep(0, 5),
          z.DAYNO=rep(0, 5),
          z.TTMIDN=rep(0, 5),
          z.WINDS=rep(0, 5),
          z.EDGED=rep(0, 5),
          z.pTREE=rep(0, 5))
          ),
        predict(sim_fit, type='response', re.form=NA, newdata=data.frame(
          fSECTION=as.factor(1:5),
          c.TURB=rep(0.4360348, 5),
          z.MINTEMP=rep(0, 5),
          z.DAYNO=rep(0, 5),
          z.TTMIDN=rep(0, 5),
          z.WINDS=rep(0, 5),
          z.EDGED=rep(0, 5),
          z.pTREE=rep(0, 5))
          )
        )
    )
  }
  return(out)
}

### test
#sim_ALL(m2z)

# Set up multicluster:
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 24), type = clusterType))
clusterExport(clust, "bats_nona")
clusterExport(clust, "m2z")
clusterExport(clust, "glmer")
clusterExport(clust, "fixef")

system.time({ temp <- clusterCall(clust, sim_ALL, m2z) })

preds_RE <- matrix(unlist(temp), byrow=T, ncol=10)

preds_RE_single <- preds_RE[,1:5]
preds_RE_mult <- preds_RE[,6:10]

par(mfrow=c(1,2))
bars_single <- barplot(preds_single$pts1_single, 
                       ylim=c(0, max(preds_RE)),
                       names.arg=c('0-100','100-200','200-300','300-400','400-500'),
                       main='Single turbine',
                       xlab='Distance band (m)',
                       ylab='Probability of a bat pass / ha'
)
arrows(bars_single, preds_single$X2.5., 
       bars_single, preds_single$X97.5., 
       code=3, length=0.1, angle=90)
arrows(bars_single+1/25, apply(preds_RE_single, 2, function(x) quantile(x, probs=c(0.025))), 
       bars_single+1/25, apply(preds_RE_single, 2, function(x) quantile(x, probs=c(0.975))), 
       code=3, length=0.1, angle=90, col='darkred', lwd=1.5)

bars_mult <- barplot(preds_mult$pts1_mult, 
                     ylim=c(0, max(preds_RE)),
                     names.arg=c('0-100','100-200','200-300','300-400','400-500'),
                     main='Single turbine',
                     xlab='Distance band (m)',
                     ylab='Probability of a bat pass / ha'
)
arrows(bars_mult, preds_mult$X2.5., 
       bars_mult, preds_mult$X97.5., 
       code=3, length=0.1, angle=90)
arrows(bars_mult+1/25, apply(preds_RE_mult, 2, function(x) quantile(x, probs=c(0.025))), 
       bars_mult+1/25, apply(preds_RE_mult, 2, function(x) quantile(x, probs=c(0.975))), 
       code=3, length=0.1, angle=90, col='darkred', lwd=1.5)
