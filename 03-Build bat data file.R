###
### THIS FILE TAKES THE RAW BAT SURVEY DATA, "CLEANS" ALL COLUMNS AND CALCULATES
###  EXTRA DATA FOR ANALYSIS FILE.
### ALSO ADDS IN HABITAT DATA PER TRANSECT SECTION.

###
### METADATA FOR FILE:
### 
### ...
### ...


library(xlsx)
library(lattice)
library(plyr)

###
###
### LOAD bat data:
bats <- read.xlsx('data/Bat data per section 2013 and 2014.xlsx',sheetIndex=1,header=T)

###
### LOAD site data:
site <- read.xlsx('data/Turbine coordinates 2013 and 2014.xlsx',sheetIndex=1,header=T)

###
### LOAD transect section data:
sections <- read.xlsx('data/Transect section data 2013 and 2014.xlsx', sheetIndex=1, header=T)

###
### LOAD habitat data:
hab <- read.csv('data/TRANSECT SECTION HABITAT DATA 2013 and 2014 PROCESSED.csv',header=T)

### 'clean' bat data:

# Wind speed as numeric:
bats$WINDS <- as.numeric(as.vector(bats$WINDS))

# fSECTION is transect section as factor:
bats$fSECTION <- factor(bats$SECTION)

# Make sure SITE and TRSCT are factors:
bats$SITE <- factor(bats$SITE)
bats$TRSCT <- factor(bats$TRSCT)

# Make 'truetime' from s_dtime:
bats$s_dtime <- strptime(bats$s_dtime, format='%Y-%m-%d %H:%M:%S')

# Calculate YEAR-SPECIFIC julian day numbers
bats$DAYNO <- as.numeric(NA)
for (i in 1:nrow(bats)) {
  org <- paste(format(bats$DATE[i],'%Y'),'01-01',sep='-')
  org <- as.Date(org, '%Y-%m-%d')
  bats$DAYNO[i] <- as.numeric(julian(bats$DATE[i],origin=org))
}

# Calculate time minutes to/since midnight w/i given night.
# First extract ISO date/time.
temp <- strptime(bats$s_dtime, format='%Y-%m-%d %H:%M:%S')
# First time difference is difference between current date time and the next day. Leaving out time info
#  in first 'format()' makes sure the returned value is 00:00:00 on the next day.
temp1 <- difftime(temp, format(temp+(3600*24),'%Y-%m-%d'), units='mins')  
# Second time difference is difference between current datetime and 'today', which leaving out time info
#  means the previous midnight.
temp2 <- difftime(temp, format(temp,'%Y-%m-%d'), units='mins')  
# Combine this info into new matrix for checking and calculating:
temp3 <- data.frame(temp, as.vector(temp1), as.vector(temp2))
# Now work out which of the two above differences in minutes is the smallest in absolute terms, 
#  and return that difference.
bats$TTMIDN <- apply(as.matrix(temp3[,2:3]),1,function(x) { x[abs(x)==min(abs(x))] })


# Find and match the number of turbines per site:
site2 <- subset(site, select=c('SITE','T'))
site3 <- site2[site2$T!='Mid',]
site3$SITE <- factor(as.vector(site3$SITE))
site_t <- ddply(site3, .(SITE), nrow)
names(site_t) <- c('SITE','NOTURB')
# Check if site names in site_t correspond to those in bats:
site_t$SITE[!(levels(site_t$SITE) %in% bats$SITE)]
# and vv:
!(levels(bats$SITE) %in% site_t$SITE)
# Match:
bats$NOTURB <- site_t$NOTURB[match(bats$SITE, site_t$SITE)]

# Total number of passes:
bats$PASSES <- bats$SOPPIP+bats$COMPIP+bats$MYOTI+bats$PIP+bats$PAUR+bats$NOCTU

# Total no. of PIPISTRELLUS passes:
bats$ALL_PIPS <- bats$COMPIP+bats$SOPPIP+bats$PIP

# Total no. of 'other' (NON-PIPISTRELLUS: MYOTIS SP., P. AURITUS, N. NOCTULA):
bats$ALL_OTHER <- bats$MYOTI+bats$PAUR+bats$NOCTU

# Remove NW transect for Islabank (too replicative)?
bats <- bats[!(bats$SITE=='Islabank' & bats$TRSCT=='NW'),]

# Remove all 'sections 6-8' (only done at few sites, not comparable).
bats <- bats[bats$SECTION<6,]

# Ignore some specific transects/sections that weren't done consistently/only once:
bats <- bats[!(bats$SITE=='Redlands' & bats$TRSCT=='E'),]
bats <- bats[!(bats$SITE=='Letham Far North T2' & bats$TRSCT=='NW' & bats$SECTION=='5'),]
bats <- bats[!(bats$SITE=='Mid Cambushinnie' & bats$TRSCT=='S' & bats$SECTION=='5'),]
bats <- bats[!(bats$SITE=='Nisbet Hill' & bats$TRSCT=='N' & bats$SECTION=='5'),]
bats <- bats[!(bats$SITE=='Wester Essendy' & bats$TRSCT=='W' & bats$SECTION=='4'),]
bats <- droplevels(bats)

### Match habitat data
bats$id <- paste(bats$SITE, bats$TRSCT, bats$SECTION, sep='-')
bats <- merge(bats, hab, 'id')
bats$AREA_ha <- bats$AREA_M2/10000

### Match coordinates:
sections$id <- paste(sections$SITE, sections$TRANSECT, sections$SECTION, sep='-')
sections2 <- subset(sections, select=c('id','MID_X','MID_Y'))
bats <- merge(bats, sections2, 'id')

### Match altitude data?

### Calculate standardised wind speed per site/visit:
out <- as.data.frame(NULL)
for(i in 1:nlevels(bats$SURVEYID)) {
  survid <- bats[bats$SURVEYID==levels(bats$SURVEYID)[i],]
  out <- rbind(out, data.frame(SURVEYID=survid$SURVEYID,
                               TRSCT=survid$TRSCT,
                               SECTION=survid$SECTION,
                               sWINDS=as.vector(scale(survid$WINDS))))
}
out$id2 <- paste(out$SURVEYID, out$TRSCT, out$SECTION, sep='-')
out <- subset(out, select=c('id2','sWINDS'))
bats$id2 <- paste(bats$SURVEYID, bats$TRSCT, bats$SECTION, sep='-')
bats <- merge(bats, out, 'id2')
bats$id <- NULL
bats$id2 <- NULL

### Which sections are smaller than half the "supposed" size (1ha)?
bats[bats$AREA_M2<5000,]
### Ignore the following "sections", as they are artefacts of the GIS transect plan
###  and don't constitute formal 'whole' sections:
### - Ferrygate-W-5
### - Letham Far North T2-W-5
### - Mosshouses N-4
### - Mosshouses NE-5
### - Mosshouses SW-5
### - Whitehills SE-4
### - Whitehills W-5
bats <- bats[!(bats$AREA_M2<5000),]
bats <- droplevels(bats)

### Which sections are more than twice the "supposed" size (1ha)?
test <- bats[bats$AREA_M2>20000,]
test <- subset(test, select=c('SITE','TRSCT','SECTION','PASSES'))
test <- unique(test)
test
### All the following are retained:
### - East Foturtune N-5; perpendicular to transect dir
### - Mosshouses-SW-3; perpendicular to transect dir
### - Nisbet Hill-N-3; perpendicular to transect dir
### - Park Cottage-E-4; perpendicular to transect dir
### - Redlands-N-3; ; perpendicular to transect dir
### - Sauchie Farm-E-5; perpendicular to transect dir
### - Stewart House-N-3; perpendicular to transect dir
### - West Lodge Balmule-S-2; perpendicular to transect dir
### - Wester Essendy-N-4; perpendicular to transect dir

### Write output for analysis:
write.csv(bats,'data/BAT DATA 2013 and 2014 MASTER.csv',row.names=F)
