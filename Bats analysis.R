library(xlsx)
library(lattice)
library(plyr)
bats <- read.xlsx('data/Bat data per section 2013 and 2014.xlsx',sheetIndex=1,header=T)

# Wind speed as numeric:
bats$WINDS <- as.numeric(as.vector(bats$WINDS))

# fSECTION is transect section as factor:
bats$fSECTION <- factor(bats$SECTION)

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
# Now work out which of the two above differences in minutes is the smallest in absolute terms, and return that
#  difference.
bats$TTMIDN <- apply(as.matrix(temp3[,2:3]),1,function(x) { x[abs(x)==min(abs(x))] })


### Match habitat data
# Load area data:
harea <- read.xlsx('data/Site habitat AREA - OS MasterMap Topo.xlsx',sheetIndex=1,header=T)
# Load line data:
larea <- read.xlsx('data/Site habitat LINE - OS MasterMap Topo.xlsx',sheetIndex=1,header=T)

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
names(habclasses[((habclasses/sum(habclasses))*100)<1])



table(factor(harea$HABCLASS))/sum(table(factor(harea$HABCLASS)))*100



harea[!(harea$SITE %in% bats$SITE),]

temp <- larea[!(larea$SITE %in% bats$SITE),]
levels(factor(as.vector(temp$SITE)))

bats$pips <- bats$SOPPIP+bats$PIP+bats$COMPIP

histogram(~pips | SITE, data=bats)
bwplot(pips ~ fSECTION | SITE, data=bats)

