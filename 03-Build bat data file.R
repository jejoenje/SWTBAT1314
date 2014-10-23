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

### 'clean' bat data:

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

### Match habitat data

### Match altitude data?

### Write output for analysis:
write.csv(bats,'data/BAT DATA 2013 and 2014 MASTER.csv',row.names=F)
