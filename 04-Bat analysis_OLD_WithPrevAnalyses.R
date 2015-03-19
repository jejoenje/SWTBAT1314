### Restored from Github repository 18/03/2015 (from commit on 17/03/2015)
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
#library(corrplot)
library(plyr)
library(reshape)
library(Hmisc)

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
#temp <- ddply(bats, .(SITE, TRSCT, fSECTION), summarise, noDates=length(unique(DATE)), noObs=length(DATE))
#write.csv(temp, 'temp.csv', row.names=F)

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
save(bats_nona, file='bats_nona.Rdata')
load('bats_nona.Rdata')

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

### First fit non-standardised models:
mod_hab1 <- glmer(OCC_PIPS ~ pBUILD + pTREE + pRDTRK + pROADS + pROADS + pRGRAS + 
                    EDGED + D_LIN + D_BUI + D_WAT + D_TRE + (1|SITE/TRSCT), offset=log(AREA_ha), 
                  data=bats_nona, family='binomial'(link='cloglog'), na.action='na.fail')
### Standardize:
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

### Model selection procedure for full habitat set, distributed across a bunch of clusters:

# clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
# clust <- try(makeCluster(getOption("cl.cores", 4), type = clusterType))
# clusterExport(clust, "bats_nona")
# clusterExport(clust, "glmer")
# clusterExport(clust, "fixef")
# clusterExport(clust, "vcov")
# clusterExport(clust, "forceSymmetric")
# 
# system.time({
#   modset_habz1 <- pdredge(mod_habz1, cluster=clust,
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
#                           , trace=2)
# })
# save(modset_habz1, file='modset_habz1.Rdata')

### Load the result from full model selection on hab variables:
load('modset_habz1.Rdata')
subset(modset_habz1, delta<4)
subset(modset_habz1, delta<10)

### So on the basis of AIC selection, let's consider D_BUI, D_WAT, EDGED and pTREE as our hab variables.

### Correlation between these four variables:
cor(subset(bats_nona, select=c('D_BUI','D_WAT','EDGED','pTREE')))
#corrplot(cor(subset(bats_nona, select=c('D_BUI','D_WAT','EDGED','pTREE'))))

### Summary table of ALL unstandardised habitat variables:
TAB_all_hab <- 
  as.data.frame(
  rbind(round(apply(bats_nona[,gsub("z.", "",names(data.frame(modset_habz1)[,2:11]))],2,mean),2),
        round(apply(bats_nona[,gsub("z.", "",names(data.frame(modset_habz1)[,2:11]))],2,median),2),
        round(apply(bats_nona[,gsub("z.", "",names(data.frame(modset_habz1)[,2:11]))],2,sd),2),
        round(apply(bats_nona[,gsub("z.", "",names(data.frame(modset_habz1)[,2:11]))],2,min),2),
        round(apply(bats_nona[,gsub("z.", "",names(data.frame(modset_habz1)[,2:11]))],2,max),2)
        ))
row.names(TAB_all_hab) <- c('Mean','Median','Std. dev.','Min.','Max.')
TAB_all_hab <- t(TAB_all_hab)


###
### OCCURRENCE MODELS (PROBABILITY OF A PASS)
###
### Basic model: GLMM with nested RE. Habitat vars based on AIC selection, and use AIC selection.
### (so that's D_BUI, D_WAT, EDGED and pTREE)

### First fit unstandardised model:
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
            data=bats_nona, family='binomial'(link='cloglog'), na.action='na.fail', 
            control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))

### Standardise predictors:
m1z <- standardize(m1)

### Build full model set to check length:
m1z_set1 <- dredge(m1z, subset=dc(z.TTMIDN, `I(z.TTMIDN^2)`), evaluate=F)
length(m1z_set1)

### Set up cluster and run model selection procedure:

# clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
# clust <- try(makeCluster(getOption("cl.cores", 4), type = clusterType))
# clusterExport(clust, "bats_nona")
# clusterExport(clust, "glmer")
# clusterExport(clust, "fixef")
# clusterExport(clust, "glmerControl")
# clusterExport(clust, "forceSymmetric")
# system.time({
#   m1z_set1 <- pdredge(m1z, cluster=clust, 
#                       subset=dc(z.TTMIDN, `I(z.TTMIDN^2)`), trace=2)  
# })
# save(m1z_set1, file='m1z_set1.Rdata')
# stopCluster(clust)

### Load the result from the full model selection procedure:
load('m1z_set1.Rdata')
subset(m1z_set1, delta<4)
subset(m1z_set1, delta<2)

### Run model averageing on delta<4 AIC set, and store the fits of each:
# m1z_av <- model.avg(m1z_set1, delta<4, fit=T)
# save(m1z_av, file='m1z_av.Rdata')

### Load results from model averaging above:
load('m1z_av.Rdata')
summary(m1z_av)



### "Simpler habitat" set: GLMM with nested RE. Only those hab vars that are present in all four d<4 hab models.
### (so that's EDGED and pTREE)

### Fit 'full model:
system.time({
  m2 <- glmer(OCC_PIPS  ~ fSECTION*TURB + 
                MINTEMP + 
                DAYNO +
                TTMIDN +
                I(TTMIDN^2) + 
                WINDS +
                EDGED +
                pTREE + 
                (1|SITE/TRSCT), offset=log(AREA_ha), 
              data=bats_nona, family='binomial'(link='cloglog'), na.action='na.fail',
              control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))
  ) })

### Standardise predictors:

### Build model set to check length:
m2z <- standardize(m2)
### Save this standardised model:
save(m2z, file='m2z.Rdata')

m2z_set1 <- dredge(m2z, subset=dc(z.TTMIDN, `I(z.TTMIDN^2)`), evaluate=F)
length(m2z_set1)

### Run model selection on 'reduced habitat set' models via cluster:
# clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
# clust <- try(makeCluster(getOption("cl.cores", 4), type = clusterType))
# clusterExport(clust, "bats_nona")
# clusterExport(clust, "glmer")
# clusterExport(clust, "fixef")
# clusterExport(clust, "glmerControl")
# clusterExport(clust, "forceSymmetric")
# 
# system.time({
#   m2z_set1 <- pdredge(m2z, cluster=clust, 
#                       subset=dc(z.TTMIDN, `I(z.TTMIDN^2)`), trace=2)  
# })
# save(m2z_set1, file='m2z_set1.Rdata')

### Load results from model selection above:
load('m2z_set1.Rdata')
subset(m2z_set1, delta<4)

### Summary model selection subset table:
TAB_mainmod_subset <- data.frame(subset(m2z_set1, delta<4))
TAB_mainmod_subset <- rbind(TAB_mainmod_subset, data.frame(m2z_set1)[which(row.names(m2z_set1)==0),])
TAB_mainmod_subset <- cbind(round(TAB_mainmod_subset[,!is.factor.df(TAB_mainmod_subset)],3), TAB_mainmod_subset[,is.factor.df(TAB_mainmod_subset)])
TAB_mainmod_subset <- TAB_mainmod_subset[,names(data.frame(subset(m2z_set1, delta<4)))]
row.names(TAB_mainmod_subset) <- NULL

### Run model averaging on above model set:
# m2z_av <- model.avg(m2z_set1, delta<4, fit=T, trace=T)
# save(m2z_av, file='m2z_av.Rdata')

### Load results from averaging above:
load('m2z_av.Rdata')
vcov(m2z_av)
summary(m2z_av)
coefTable(m2z_av, full=TRUE) # "Fully" averaged parameters; i.e. with shrinkage (zero-method)
coefTable(m2z_av, full=FALSE) # Averaged parameters WITHOUT shrinkage (natural average)

### Store list of d<4 models:
#m2z_set1_mods <- get.models(m2z_set1, subset=delta<4, cluster=clust)
#save(m2z_set1_mods, file='m2z_set1_mods.Rdata')
load('m2z_set1_mods.Rdata')

### Table of summary statistics of UNSTANDARDISED inputs to model m2z (so from model m2):
m2_dat <- attr(m2,'frame')[,2:10]
m2_dat_covs <- m2_dat[,!is.factor.df(m2_dat)]
TAB_main_model_inputs_summary <- as.data.frame(rbind(
  apply(m2_dat_covs, 2, mean),
  apply(m2_dat_covs, 2, median),
  apply(m2_dat_covs, 2, sd),
  apply(m2_dat_covs, 2, min),
  apply(m2_dat_covs, 2, max)
  ))
row.names(TAB_main_model_inputs_summary) <- c('Mean','Median','Std. dev','Min.','Max.')
TAB_main_model_inputs_summary <- t(round(TAB_main_model_inputs_summary,2))

### Summary tables of 'model averaged' parameter estimates and parameter importance ('weight') for presentation.
TAB_coefs_av <- round(as.data.frame(coefTable(m2z_av, full=TRUE)),3)
TAB_imp <- summary(m2z_av)$importance

###
### PLOT PREDICTIONS: prediction interval by sim():
###

linkinv <- m2z@resp$family$linkinv

# Standardised values for single and multiple turb
c.turb <- unique(m2z@frame$c.TURB)[order(unique(m2z@frame$c.TURB))]

wmean <- function(x, w=NULL) {
  if(is.null(w)) { w <- rep(1/length(x), length(x)) }
  return(as.vector(as.numeric(sum(x*w))))
}

# Use manual prediction averaging to calculate both point predictions and simmed credible intervals:

# Prediction data frame for SINGLE:
nd_single <- data.frame(
  fSECTION=factor(1:5),
  c.TURB=rep(c.turb[1],5),
  z.MINTEMP=rep(0,5),
  z.DAYNO=rep(0,5),
  z.TTMIDN=rep(0,5),
  z.WINDS=rep(0,5),
  z.EDGED=rep(0,5),
  z.pTREE=rep(0,5),
  AREA_ha=rep(1,5)
  )
# Prediction data frame for MULTIPLE:
nd_multiple <- data.frame(
  fSECTION=factor(1:5),
  c.TURB=rep(c.turb[2],5),
  z.MINTEMP=rep(0,5),
  z.DAYNO=rep(0,5),
  z.TTMIDN=rep(0,5),
  z.WINDS=rep(0,5),
  z.EDGED=rep(0,5),
  z.pTREE=rep(0,5),
  AREA_ha=rep(1,5)
)
### Test predictions with simmed parameters from a single model:
# m <- m2z_set1_mods[[1]]
# m_sim <- sim(m, 1000)
# m_sim_fe <- attr(m_sim, 'fixef')
# predict(m, newdata=nd_single, newparams=list(beta=m_sim_fe[1,], theta=getME(m, 'theta')), type='response', re.form=NA)
# predict(m, newdata=nd_single, newparams=list(beta=m_sim_fe[2,], theta=getME(m, 'theta')), type='response', re.form=NA)

### Now, loop through each of the best models.
### For each, simulate K parameters.
### For each, predict response based on K parameters.
### Store all K parameters.

### Convenience function to calculate mean and 5% quantiles:
mnqtl <- function(x) {
  c(mean(x), quantile(x, probs=c(0.025,0.975)))
}
mn_single <- as.data.frame(NULL)
lo_single <- as.data.frame(NULL)
hi_single <- as.data.frame(NULL)
mn_multiple <- as.data.frame(NULL)
lo_multiple <- as.data.frame(NULL)
hi_multiple <- as.data.frame(NULL)
for (i in 1:length(m2z_set1_mods)) {
  m <- m2z_set1_mods[[i]]
  m_sim <- sim(m, 1000)
  m_sim_fe <- attr(m_sim, 'fixef')
  m_preds_single <- apply(m_sim_fe, 1, function(x) predict(m, 
                                                    newdata=nd_single, 
                                                    newparams=list(beta=x, theta=getME(m,'theta')), 
                                                    type='response', 
                                                    re.form=NA))
  m_preds_multiple <- apply(m_sim_fe, 1, function(x) predict(m, 
                                                           newdata=nd_multiple, 
                                                           newparams=list(beta=x, theta=getME(m,'theta')), 
                                                           type='response', 
                                                           re.form=NA))
  mn_single <- rbind(mn_single, apply(m_preds_single, 1, mean))
  lo_single <- rbind(lo_single, apply(m_preds_single, 1, function(x) quantile(x, probs=c(0.025))))
  hi_single <- rbind(hi_single, apply(m_preds_single, 1, function(x) quantile(x, probs=c(0.975))))
  mn_multiple <- rbind(mn_multiple, apply(m_preds_multiple, 1, mean))
  lo_multiple <- rbind(lo_multiple, apply(m_preds_multiple, 1, function(x) quantile(x, probs=c(0.025))))
  hi_multiple <- rbind(hi_multiple, apply(m_preds_multiple, 1, function(x) quantile(x, probs=c(0.975))))
}
### So for want of a more elegant way of summarising, calculate the result as means of means and means of
###  quantiles:
wtd_mn_single <- apply(mn_single, 2, function(x) wtd.mean(x, weights=subset(m2z_set1,delta<4)$weight))
wtd_lo_single <- apply(lo_single, 2, function(x) wtd.mean(x, weights=subset(m2z_set1,delta<4)$weight))
wtd_hi_single <- apply(hi_single, 2, function(x) wtd.mean(x, weights=subset(m2z_set1,delta<4)$weight))
wtd_mn_multiple <- apply(mn_multiple, 2, function(x) wtd.mean(x, weights=subset(m2z_set1,delta<4)$weight))
wtd_lo_multiple <- apply(lo_multiple, 2, function(x) wtd.mean(x, weights=subset(m2z_set1,delta<4)$weight))
wtd_hi_multiple <- apply(hi_multiple, 2, function(x) wtd.mean(x, weights=subset(m2z_set1,delta<4)$weight))

FIG1 <- function() {
  plot(1:5, c(0,0,0,1,1), xlim=c(0.5,5.5), ylim=c(0,.6), type='n', xaxt='n',
       xlab="Distance from turbine(s)",
       ylab="Predicted probability of a bat pass / ha")
  axis(1, at=1:5, labels=c('0-100m','100-200m','200-300m','300-400m','400-500m'))
  offset <- 0.1
  for(i in 1:5){
    lines(c(i-offset,i-offset),c(wtd_lo_single[i], wtd_hi_single[i]))
    lines(c(i+offset,i+offset),c(wtd_lo_multiple[i], wtd_hi_multiple[i]))
  }
  points(1:5-offset, wtd_mn_single, cex=1.5, pch=21, col='black', bg='white')
  points(1:5+offset, wtd_mn_multiple, cex=1.5, pch=16)
  #Compare predictions above with automated predictions (predict.averaging()):
  #points(1:5-offset, pred_av_single, pch=16, col='blue')
  #points(1:5+offset, pred_av_multiple, pch=16, col='red')
  
}; FIG1()





# Prediction intervals cf. Bolker:

p_full_single_se <- sqrt(diag(p_full_single %*% vcov(m2z) %*% t(p_full_single)))
p_full_mult_se <- sqrt(diag(p_full_mult %*% vcov(m2z) %*% t(p_full_mult)))
p_single_se <- sqrt(diag(p_single %*% vcov(m2z_av) %*% t(p_single)))
p_mult_se <- sqrt(diag(p_multp %*% vcov(m2z_av) %*% t(p_multp)))

### OBSERVED proportions:
# Calculate observered mean proportions 0/1 for single and multiple turbines:
bats_nona_single <- bats_nona[bats_nona$TURB=='single',]
bats_nona_multp <- bats_nona[bats_nona$TURB=='multiple',]
bats_nona_single <- droplevels(bats_nona_single); bats_nona_multp <- droplevels(bats_nona_multp)

obs_single <- tapply(bats_nona_single$OCC_PIPS, bats_nona_single$fSECTION, sum)/
  tapply(bats_nona_single$OCC_PIPS, bats_nona_single$fSECTION, length)

obs_multp <- tapply(bats_nona_multp$OCC_PIPS, bats_nona_multp$fSECTION, sum)/
  tapply(bats_nona_multp$OCC_PIPS, bats_nona_multp$fSECTION, length)

par(mfrow=c(2,1))
bplot_single <- barplot(preds_single, beside=T, ylim=c(0, max(c(preds_single_upr, preds_mult_upr))), 
                        col=rep('lightgrey',5),
                        main='Single turbine', 
                        names.arg=c('Full','Averaged (NA)', 'Averaged (Shrunk)', 'Averaged (predict())'))
mtext(text=1:5, side=1, line=0, at=bplot_single)
for(i in 1:nrow(preds_single_lwr)) {
  for(j in 1:ncol(preds_single_lwr)) {
    arrows(bplot_single[i,j],preds_single_lwr[i,j],bplot_single[i,j],preds_single_upr[i,j],
           code=3, length=0.1, angle=90)
  } 
}
points(as.vector(bplot_single),rep(obs_single,4),pch=16)

bplot_mult <- barplot(preds_mult, beside=T, ylim=c(0, max(c(preds_single_upr, preds_mult_upr))), 
                      col=rep('lightgrey',5),
                      main='Multiple turbines', 
                      names.arg=c('Full','Averaged (NA)', 'Averaged (Shrunk)', 'Averaged (predict())'))
mtext(text=1:5, side=1, line=0, at=bplot_mult)
for(i in 1:nrow(preds_mult_lwr)) {
  for(j in 1:ncol(preds_mult_lwr)) {
    arrows(bplot_mult[i,j],preds_mult_lwr[i,j],bplot_mult[i,j],preds_mult_upr[i,j],
           code=3, length=0.1, angle=90)
  } 
}
points(as.vector(bplot_mult),rep(obs_multp,4),pch=16)


av_bars <- barplot(matrix(av_wind$fit,2,2), beside=T, col=c(shd1, shd2), 
                   ylim=c(0, max(hi_wind$upr)), names.arg=c('Control','Turbine'), 
                   main='Average wind speed (5.9 m/s)',
                   ylab='Predicted number of S. Pip passes per night')

legend(1, 500, legend=c('Before','After'), bty='n', 
       xjust=0, col=c(shd1, shd2), fill=c(shd1, shd2), x.intersp=0.25, text.width=0.25)





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
