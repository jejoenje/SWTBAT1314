###
###
### ANALYSE HABITAT DATA PER TRANSECT SECTION
###
### 

rm(list=ls())
hab <- read.csv('data/TRANSECT SECTION HABITAT DATA 2013 and 2014 MASTER.csv',header=T)

### First express all m2 habitat cover to proportions:
hab$pBUILD <- hab$BUILD/hab$AREA_M2
hab$pCTREE <- hab$CTREE/hab$AREA_M2
hab$pGREYS <- hab$GREYS/hab$AREA_M2
hab$pLANDF <- hab$LANDF/hab$AREA_M2
hab$pNTREE <- hab$NTREE/hab$AREA_M2
hab$pRDTRK <- hab$RDTRK/hab$AREA_M2
hab$pROADS <- hab$ROADS/hab$AREA_M2
hab$pRGRAS <- hab$RGRAS/hab$AREA_M2
hab$pSCRUB <- hab$SCRUB/hab$AREA_M2
hab$pWATER <- hab$WATER/hab$AREA_M2

### Second, can I 'weed out' some habitat cover classifiers on the basis of too
###  little information content?

# Visualise individual land cover data:
par(mfrow=c(3,4))
hist(hab$pBUILD, breaks=100)
hist(hab$pCTREE, breaks=100)
hist(hab$pGREYS, breaks=100)
hist(hab$pLANDF, breaks=100)
hist(hab$pNTREE, breaks=100)
hist(hab$pRDTRK, breaks=100)
hist(hab$pROADS, breaks=100)
hist(hab$pRGRAS, breaks=100)
hist(hab$pSCRUB, breaks=100)
hist(hab$pWATER, breaks=100)

# Visualise same as above but w/o the zero's:
par(mfrow=c(3,4))
hist(hab$pBUILD[hab$pBUILD!=0], breaks=100)
hist(hab$pCTREE[hab$pCTREE!=0], breaks=100)
hist(hab$pGREYS[hab$pGREYS!=0], breaks=100)
hist(hab$pLANDF[hab$pLANDF!=0], breaks=100)
hist(hab$pNTREE[hab$pNTREE!=0], breaks=100)
hist(hab$pRDTRK[hab$pRDTRK!=0], breaks=100)
hist(hab$pROADS[hab$pROADS!=0], breaks=100)
hist(hab$pRGRAS[hab$pRGRAS!=0], breaks=100)
hist(hab$pSCRUB[hab$pSCRUB!=0], breaks=100)
hist(hab$pWATER[hab$pWATER!=0], breaks=100)

# Percentage of zero's for each of the above classifiers:
round(sum(hab$pBUILD!=0)/length(hab$pBUILD),3)
round(sum(hab$pCTREE!=0)/length(hab$pCTREE),3)
round(sum(hab$pGREYS!=0)/length(hab$pGREYS),3)
round(sum(hab$pLANDF!=0)/length(hab$pLANDF),3)
round(sum(hab$pNTREE!=0)/length(hab$pNTREE),3)
round(sum(hab$pRDTRK!=0)/length(hab$pRDTRK),3)
round(sum(hab$pROADS!=0)/length(hab$pROADS),3)
round(sum(hab$pRGRAS!=0)/length(hab$pRGRAS),3)
round(sum(hab$pSCRUB!=0)/length(hab$pSCRUB),3)
round(sum(hab$pWATER!=0)/length(hab$pWATER),3)

### So, the following conclusions:
### 1. pGREYS is pointless to include as it is a 'remaining' classifier and therefore
###     is not infortmative on habitat type.
### 2. pLANDF is also not really informative, and is non-zero in only 6.3% of cases.
### 3. pSCRUB and pWATER also occur in <10% of cases, and are therefore of limited value.
###
### 4.This leaves us with proportional cover variables pBUILD, pCTREE, pNTREE, pRDTRK, pROADS and pRGRAS.
###
### 5. For simplicity, will only consider total extent of tree cover, pTREE = pCTREE+pNTREE

hab$pTREE <- hab$pCTREE+hab$pNTREE

### 6. This leaves us with pBUILD, pTREE, pRDTRK, pROADS, and pGRAS.

### Now what about line length (LLENG)?

# First express this as a measure of 'edge density' by dividing by the total area.
#  (i.e. now in m/m2)

hab$EDGED <- hab$LLENG/hab$AREA_M2
par(mfrow=c(1,1))
# Visualise:
hist(hab$EDGED, breaks=100)


### Finally check and visualise the distance-to variables:
par(mfrow=c(2,2))
hist(hab$D_LIN, breaks=100)
hist(hab$D_BUI, breaks=100)
hist(hab$D_WAT, breaks=100)
hist(hab$D_TRE, breaks=100)


### Prepare data for PCA no. 1:
hab2 <- subset(hab, select=c('pBUILD','pTREE','pRDTRK',
                             'pROADS','pRGRAS','EDGED',
                             'D_LIN','D_BUI','D_WAT','D_TRE'))

### Run PCA:
hab_pca1 <- prcomp(hab2, scale=T, center=T)
summary(hab_pca1)
hab_pca_scores1 <- round(hab_pca1$rotation[,1:4],3)
hab_pca_scores1
# Interpretation of PC1-4:
# PC1:
#  More roads/tracks cover
#  Higher edge density
#  More roadside cover
#  Shorter distance to linear features
#  Shorter distance to buildings
#  == "Edgey habitat with roads, but close to buildings"
#
# PC2 (its inverse!):
#  Higher tree cover
#  Shorter distance from water
#  Shorter distance from trees
#  More rough grassland
# == "Wet woodlands with some rough grassland cover"
#
# PC3:
#  Fewer buildings
#  More roadside cover
#  Longer distance from buildings
#  More road/track cover
# == Isolated roads with few surrounding features, relatively open habitat
#
# PC4 (its inverse):
#  Higher tree cover
#  Less rough grassland
#  Longer distance from water
# == Dry woodlands


### Have a look at the correlation within this subset of habitat data:
library(corrplot)
cm_hab2 <- cor(hab2)
corrplot(cm_hab2, method='circle')


### Prepare data for PCA no. 2:
hab3 <- subset(hab, select=c('EDGED','D_LIN','D_BUI','D_WAT','D_TRE'))
### Have a look at the correlation within this subset of habitat data:
cm_hab3 <- cor(hab3)
corrplot(cm_hab3, method='circle')
hab_pca2 <- prcomp(hab3, scale=T, center=T)
summary(hab_pca2)
hab_pca_scores2 <- round(hab_pca2$rotation[,1:4],3)
hab_pca_scores2
# Interpretation of PC1-4:
# PC1:
#  Higher edge density
#  Shorter distance to linear features
#  Shorter distance to buildings
#  Shorter distance to trees
# == 'Structurally complex habitat close to buildings and trees'
#
# PC2:
#  Longer distance to buildings
#  Shorter distance to water
#  Shorter distance to trees
# == 'Wet habitat near trees'
#
# PC3:
#  Shorter distance to linear features
#  Longer distance to buildings
#  Longer distance to trees
# == 'Habitat with linear features but away from buildings or trees' (i.e. roads)  
#
# PC4 (inverse):
#  Longer distance to buildings
#  Longer distance to water
#  Shorter distance to trees
#  Longer distance to linear features
# == 'Dry habitat near trees'

