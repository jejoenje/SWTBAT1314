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

### Convenience function to calculate mean and 5% quantiles:
mnqtl <- function(x) {
  c(mean(x), quantile(x, probs=c(0.025,0.975)))
}

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
# cor(subset(bats_nona, select=c('D_BUI','D_WAT','EDGED','pTREE')))
# corrplot(cor(subset(bats_nona, select=c('D_BUI','D_WAT','EDGED','pTREE'))))

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
save(TAB_all_hab, file='TAB_all_hab.Rdata')

###
### OCCURRENCE MODELS (PROBABILITY OF A PASS)
###
### Basic model: GLMM with nested RE. Habitat vars based on AIC selection, and use AIC selection.
### (so that's D_BUI, D_WAT, EDGED and pTREE)

### First fit unstandardised model:
# m1 <- glmer(OCC_PIPS  ~ fSECTION*TURB + 
#               MINTEMP + 
#               DAYNO +
#               TTMIDN +
#               I(TTMIDN^2) + 
#               WINDS +
#               D_BUI + 
#               D_WAT +
#               EDGED +
#               pTREE + 
#               (1|SITE/TRSCT), offset=log(AREA_ha), 
#             data=bats_nona, family='binomial'(link='cloglog'), na.action='na.fail', 
#             control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))

### Standardise predictors:
# m1z <- standardize(m1)

### Build full model set to check length:
# m1z_set1 <- dredge(m1z, subset=dc(z.TTMIDN, `I(z.TTMIDN^2)`), evaluate=F)
# length(m1z_set1)

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
# load('m1z_set1.Rdata')
# subset(m1z_set1, delta<4)
# subset(m1z_set1, delta<2)

### Run model averageing on delta<4 AIC set, and store the fits of each:
# m1z_av <- model.avg(m1z_set1, delta<4, fit=T)
# save(m1z_av, file='m1z_av.Rdata')

### Load results from model averaging above:
# load('m1z_av.Rdata')
# summary(m1z_av)



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
save(m2, file='m2.Rdata')

### Standardise predictors:

### Build model set to check length:
m2z <- standardize(m2)
### Save this standardised model:
save(m2z, file='m2z.Rdata')

### Null model (intercept only) - need this for easy reference later:
m2_null <- glmer(OCC_PIPS  ~ 1 + 
                   (1|SITE/TRSCT), offset=log(AREA_ha), 
                 data=bats_nona, family='binomial'(link='cloglog'), na.action='na.fail',
                 control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
save(m2_null, file='m2_null.Rdata')

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
### Add marginal and conditional R2:
TAB_mainmod_subset$mR2 <- NA
TAB_mainmod_subset$cR2 <- NA
for(i in 1: length(m2z_set1_mods)) {
  TAB_mainmod_subset$mR2[i] <- r2mm(m2z_set1_mods[[i]])$R2[1]
  TAB_mainmod_subset$cR2[i] <- r2mm(m2z_set1_mods[[i]])$R2[2]
}
TAB_mainmod_subset <- rbind(TAB_mainmod_subset,cbind(data.frame(m2z_set1)[which(row.names(m2z_set1)==0),],mR2=r2mm(m2_null)$R2[1],cR2=r2mm(m2_null)$R2[2]))
TAB_mainmod_subset <- cbind(round(TAB_mainmod_subset[,!is.factor.df(TAB_mainmod_subset)],3), TAB_mainmod_subset[,is.factor.df(TAB_mainmod_subset)])
TAB_mainmod_subset <- TAB_mainmod_subset[,c(names(data.frame(subset(m2z_set1, delta<4))),'mR2','cR2')]
row.names(TAB_mainmod_subset) <- NULL
save(TAB_mainmod_subset, file='TAB_mainmod_subset.Rdata')

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
save(TAB_main_model_inputs_summary, file='TAB_main_model_intputs_summary.Rdata')

### Summary tables of 'model averaged' parameter estimates and parameter importance ('weight') for presentation.
TAB_coefs_av <- round(as.data.frame(coefTable(m2z_av, full=TRUE)),3)
TAB_imp <- summary(m2z_av)$importance
save(TAB_coefs_av, file='TAB_coefs_av.Rdata')
save(TAB_imp, file='TAB_imp.Rdata')

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

# mn_single <- as.data.frame(NULL)
# lo_single <- as.data.frame(NULL)
# hi_single <- as.data.frame(NULL)
# mn_multiple <- as.data.frame(NULL)
# lo_multiple <- as.data.frame(NULL)
# hi_multiple <- as.data.frame(NULL)
# for (i in 1:length(m2z_set1_mods)) {
#   m <- m2z_set1_mods[[i]]
#   m_sim <- sim(m, 1000)
#   m_sim_fe <- attr(m_sim, 'fixef')
#   m_preds_single <- apply(m_sim_fe, 1, function(x) predict(m, 
#                                                     newdata=nd_single, 
#                                                     newparams=list(beta=x, theta=getME(m,'theta')), 
#                                                     type='response', 
#                                                     re.form=NA))
#   m_preds_multiple <- apply(m_sim_fe, 1, function(x) predict(m, 
#                                                            newdata=nd_multiple, 
#                                                            newparams=list(beta=x, theta=getME(m,'theta')), 
#                                                            type='response', 
#                                                            re.form=NA))
#   mn_single <- rbind(mn_single, apply(m_preds_single, 1, mean))
#   lo_single <- rbind(lo_single, apply(m_preds_single, 1, function(x) quantile(x, probs=c(0.025))))
#   hi_single <- rbind(hi_single, apply(m_preds_single, 1, function(x) quantile(x, probs=c(0.975))))
#   mn_multiple <- rbind(mn_multiple, apply(m_preds_multiple, 1, mean))
#   lo_multiple <- rbind(lo_multiple, apply(m_preds_multiple, 1, function(x) quantile(x, probs=c(0.025))))
#   hi_multiple <- rbind(hi_multiple, apply(m_preds_multiple, 1, function(x) quantile(x, probs=c(0.975))))
# }

### Load the data from above simulations: 
load('FIG_predictions1.Rdata')
load('m2z_set1.Rdata')

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

obsdat <- subset(bats_nona, select=c('OCC_PIPS','fSECTION','TURB'))

FIG2 <- function() {
  bars <- barplot(tapply(obsdat$OCC_PIPS, list(obsdat$TURB, obsdat$fSECTION), mean),beside=T,ylim=c(0,0.5), 
                  col=c('white','grey'), xaxt='n', xlab='Distance band', ylab='Probability of a Pipistrelle bat pass / ha')
  axis(1, at=apply(bars,2,mean), labels=c('0-100m','100-200m','200-300m','300-400m','400-500m'))
  for(i in 1:length(bars[1,])) {
#     lines(c(bars[1,][i],bars[1,][i]), c(wtd_lo_single[i],wtd_hi_single[i]))
#     lines(c(bars[2,][i],bars[2,][i]), c(wtd_lo_multiple[i],wtd_hi_multiple[i]))
      arrows(bars[1,][i], wtd_lo_single[i] ,bars[1,][i], wtd_hi_single[i], angle=90, code=3, length=0.1)
      arrows(bars[2,][i], wtd_lo_multiple[i],bars[2,][i], wtd_hi_multiple[i], angle=90, code=3, length=0.1)    
  }
  points(bars[1,], wtd_mn_single, cex=1.5, pch=21, col='black', bg='lightgrey')
  points(bars[2,], wtd_mn_multiple, cex=1.5, pch=21, col='black', bg='black')
}; FIG2()
