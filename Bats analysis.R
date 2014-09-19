library(xlsx)
library(lattice)
bats <- read.xlsx('data/Bat data per section 2013 and 2014.xlsx',sheetIndex=1,header=T)

# Wind speed as numeric:
bats$WINDS <- as.numeric(as.vector(bats$WINDS))

# fSECTION is transect section as factor:
bats$fSECTION <- factor(bats$SECTION)

# Make 'truetime' from s_dtime:
bats$s_dtime <- strptime(bats$s_dtime, format='%Y-%m-%d %H:%M:%S')

# For the 2013 surveys, times were only recorded specific to whole transects.
# Make sure that the 2014 data is 'downsampled' to this resolution.
# Construct factor to indicate transect-visit-site:
bats$SiteTrans <- with(bats, paste(SURVEYID,TRSCT,sep='-'))
bats$SiteTrans <- factor(bats$SiteTrans)
# Calculate the earliest time recorded per transect-visit-site:
temp <- ddply(bats, .(SiteTrans), summarise, startt=min(s_dtime))



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
bats$ttMidn <- apply(as.matrix(temp3[,2:3]),1,function(x) { x[abs(x)==min(abs(x))] })




# Julian day number since 1 Jan
bats$jDATE <- round(difftime(temp, as.Date(paste(as.character(format(temp,'%Y')),'-01-01',sep=''),
                                           format='%Y-%m-%d'),units='days'),0)



harea <- read.xlsx('data/Site habitat AREA - OS MasterMap Topo.xlsx',sheetIndex=1,header=T)


bats$pips <- bats$SOPPIP+bats$PIP+bats$COMPIP

histogram(~pips | SITE, data=bats)
bwplot(pips ~ fSECTION | SITE, data=bats)

