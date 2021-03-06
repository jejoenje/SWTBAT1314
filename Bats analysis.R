
####
#### THIS FILE IS NOW DEFUNCT - FILE BUILDING AND ANALYSES PERFORMED IN 01-Build habitat.., 
#### 02-Habitat analysis.., etc

library(xlsx)
library(lattice)
library(plyr)

###
###
### LOAD habitat data
# Load area data:
harea <- read.xlsx('data/Site habitat AREA - OS MasterMap Topo.xlsx',sheetIndex=1,header=T)
# Load line data:
larea <- read.xlsx('data/Site habitat LINE - OS MasterMap Topo.xlsx',sheetIndex=1,header=T)
# Load section data:
sdat <- read.xlsx('data/Transect section data 2013 and 2014.xlsx',sheetIndex=1,header=T)
# Load distance data:
dist <- read.xlsx('data/Site habitat DISTANCE - OS MasterMap Topo.xlsx',sheetIndex=1,header=T)

###
###
### SECTION data:
# Create site-transect-section identifier in section data:
sdat$SiteTBuffer <- paste(sdat$SITE, sdat$TRANSECT, sdat$SECTION, sep='-')
sdat$SiteTBuffer <- factor(sdat$SiteTBuffer)


###
###
### AREA data.
# Summarise area data:
harea$SiteTBuffer <- paste(harea$SITE, harea$Transect, harea$Buffer, sep='-')
harea$SiteTBuffer <- factor(harea$SiteTBuffer)

# Extract first classification from descrTerms and store in new dTerm2:
temp <- strsplit(as.character(as.vector(harea$descrTerms)),'|', fixed=T)
harea$dTerm2 <- unlist(lapply(temp, function(x) return(x[1])))
# Make sure dTerm2 NA's are blanks:
harea$dTerm2[is.na(harea$dTerm2)] <- ''
# Treat "Nonconiferous Trees (Scattered)" as "Nonconiferous Trees":
harea$dTerm2[harea$dTerm2=="Nonconiferous Trees (Scattered)"] <- "Nonconiferous Trees"
# Treat "Coniferous Trees (Scattered)" as "Coniferous Trees":
harea$dTerm2[harea$dTerm2=="Coniferous Trees (Scattered)"] <- "Coniferous Trees"
# Ignore dTerm classifications Archway, Cliff, Heath, Marsh Reeds or Saltmarsh, Pylon and Slope,
#  all of which comprise <1% of total habitat records:
harea$dTerm2[harea$dTerm2=="Archway"] <- ''
harea$dTerm2[harea$dTerm2=="Cliff"] <- ''
harea$dTerm2[harea$dTerm2=="Heath"] <- ''
harea$dTerm2[harea$dTerm2=="Marsh Reeds Or Saltmarsh"] <- ''
harea$dTerm2[harea$dTerm2=="Pylon"] <- ''
harea$dTerm2[harea$dTerm2=="Slope"] <- ''

# As above, extract first classification from descGroups and store in new dGroup2:
temp <- strsplit(as.character(as.vector(harea$descGroups)),'|', fixed=T)
harea$dGroup2 <- unlist(lapply(temp, function(x) return(x[1])))
# Combine dGroup2 and dTerm2 into new classifier.
harea$HABCLASS <- as.character(NA)
harea$HABCLASS[harea$dTerm2==''] <- harea$dGroup2[harea$dTerm2=='']
harea$HABCLASS[harea$dTerm2!=''] <- paste(harea$dGroup2[harea$dTerm2!=''], harea$dTerm2[harea$dTerm2!=''], sep='-')
# Summarise new classifcations:
habclasses <- table(factor(harea$HABCLASS))
# Ignore those classes that occur in <1% of cases:
ignore_harea <-names(habclasses[((habclasses/sum(habclasses))*100)<1])
harea2 <- harea[!(harea$HABCLASS %in% ignore_harea),]
table(factor(harea2$HABCLASS))
# Of remaining categories, treat 'General Surface-Multi Surface' as 'General Surface':
harea2$HABCLASS[harea2$HABCLASS=='General Surface-Multi Surface'] <- 'General Surface'
# Of remaining categories, treat 'Road Or Track-Track' as 'Road Or Track':
harea2$HABCLASS[harea2$HABCLASS=='Road Or Track-Track'] <- 'Road Or Track'
# Rename some categories to shorter names:
harea2$HABCLASS[harea2$HABCLASS=='Inland Water'] <- 'Water'
harea2$HABCLASS[harea2$HABCLASS=='Natural Environment-Coniferous Trees'] <- 'Coniferous Trees'
harea2$HABCLASS[harea2$HABCLASS=='Natural Environment-Nonconiferous Trees'] <- 'Nonconiferous Trees'
harea2$HABCLASS[harea2$HABCLASS=='Natural Environment-Rough Grassland'] <- 'Rough Grassland'
harea2$HABCLASS[harea2$HABCLASS=='Natural Environment-Scrub'] <- 'Scrub'
table(factor(harea2$HABCLASS))
harea2$HABCLASS <- factor(harea2$HABCLASS)
# Rename again to short names:
levels(harea2$HABCLASS) <- c('BUILD','CTREE','GREYS','LANDF','NTREE','RDTRK','ROADS','RGRAS','SCRUB','WATER')
table(factor(harea2$HABCLASS))

# Cross-tabulate Site/Transect/Buffer and Habclass, summing the area of each hab type:
harea3 <- tapply(harea2$TypeArea, list(harea2$SiteTBuffer, harea2$HABCLASS), 
                 function(x) return(sum(x, na.rm=T)))
harea3 <- data.frame(harea3)

# Create columns for sdat file matching above habitat categories:
#temp <- data.frame(matrix(NA, nrow=nrow(sdat), ncol=length(levels(harea2$HABCLASS))))
#names(temp) <- levels(harea2$HABCLASS)
#sdat <- cbind(sdat, temp)

# Check Site/Transect/Buffer names in habitat area data file are all ok with those in sdat:
harea3[!(row.names(harea3) %in% sdat$SiteTBuffer),]

# Match harea3 data columns to sdat:
sdat2 <- cbind(sdat, harea3[match(sdat$SiteTBuffer,row.names(harea3)),])
row.names(sdat2) <- NULL
#write.csv(sdat2, 'data/Temp/Transect Section habitat data 2013 and 2014 WORKING.csv')

### Match habitat AREA data per section to BATS data:
# First create Site/Transect/Section identifier in bats data.
bats$SiteTBuffer <- with(bats, paste(SITE, TRSCT, SECTION, sep='-'))
# These are the SiteTBuffer ID's that occur in the bats data but not in the site/transect/section hab data:
bats[!(bats$SiteTBuffer %in% sdat2$SiteTBuffer),'SiteTBuffer']
# !!! The list of transect sections above are defunct ones that need to be ignored in later analysis!


###
###
### LINE data.
# Summarise line data:
larea$SiteTBuffer <- paste(larea$SITE, larea$Transect, larea$Buffer, sep='-')
larea$SiteTBuffer <- factor(larea$SiteTBuffer)

# Check site/transect/section ID'ers:
levels(larea$SiteTBuffer)[!(levels(larea$SiteTBuffer) %in% sdat$SiteTBuffer)]

# Sum line length per site/transect/buffer:
sumLength <- ddply(larea, .(SiteTBuffer), summarise, sumLength=sum(LineLength))

# Match sumLengths to sdat2:
sdat2$LLENG <- sumLength$sumLength[match(sdat2$SiteTBuffer, sumLength$SiteTBuffer)]
#write.csv(sdat2, 'data/Temp/Transect Section habitat data 2013 and 2014 WORKING.csv')


###
###
### DISTANCE data.
dist$SiteTBuffer <- paste(dist$SITE, dist$Transect, dist$Buffer, sep='-')
dist$SiteTBuffer <- factor(dist$SiteTBuffer)

# Check site/transect/section ID'ers:
levels(dist$SiteTBuffer)[!(levels(dist$SiteTBuffer) %in% sdat$SiteTBuffer)]
dist[!(dist$SITE %in% sdat2$SITE),]

# Remove cols from dist that are not needed for match:
dist$SITE <- NULL
dist$Buffer <- NULL
dist$Transect <- NULL
dist$area <- NULL
row.names(dist) <- dist$SiteTBuffer
dist$SiteTBuffer <- NULL
# Change col names in dist:
names(dist) <- c('D_LIN','D_BUI', 'D_WAT', 'D_TRE')

# Match dist data to sdat2:
sdat3 <- cbind(sdat2, dist[match(sdat2$SiteTBuffer, row.names(dist)),])
#write.csv(sdat3, 'data/Temp/Transect Section habitat data 2013 and 2014 WORKING.csv')

### In the consolidated habitat data file, make sure NA cover == 0:
sdat3[,levels(harea2$HABCLASS)][is.na(sdat3[,levels(harea2$HABCLASS)])] <- 0

### Same applies to LLENG:
sdat3$LLENG[is.na(sdat3$LLENG)] <- 0
#write.csv(sdat3, 'data/Temp/Transect Section habitat data 2013 and 2014 WORKING2.csv')

### Remove unnecessary columns:
sdat3$Mairi_TransectName <- NULL
sdat3$SiteTBuffer <- NULL

###
###
### STORE 