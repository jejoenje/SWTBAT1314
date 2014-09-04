library(gdata)
library(sp)
library(rgdal)
# Set root folder with data:
rootf <- '../2014/Bats'
# List folders in data folder - should ONLY be individual sites (check!)
sfolders <- list.files(rootf)
# Count no. folders (== no. of sites)
nfolders <- length(sfolders)

# Start with an empty data frame for all site data:
alldat <- as.data.frame(NULL)
alldat_inc_noise <- as.data.frame(NULL)

# Loop through each site i:
for (i in 1:nfolders) {
  # Within site i, list folders - should ONLY be individual dates (check!)
  sdates <- list.files(paste(rootf,sfolders[i],sep='/'))
  # Loop through each date j in site i:
  for (j in 1:length(sdates)) {
    # List all files within date j in site i:
    flist <- list.files(paste(rootf, sfolders[i], sdates[j], sep='/'))
    # Read Anabat header file for date j in site i:
    headerfile <- read.csv(paste(rootf, sfolders[i], sdates[j], 'header.txt', sep='/'),header=T, sep='\t')
    # Create a copy of the headerfile data:
    nwheaderfile <- headerfile
    # 'Clean' individual Anabat file names (rows in header file). 
    # Remove extensions:
    nwheaderfile$Name <- sub('_000.00#','',nwheaderfile$Name)
    # Remove subsequent _, _N_, or _0_ so that file names can be cleanly split:
    nwheaderfile$Name <- sub('___0_', '___', nwheaderfile$Name)
    nwheaderfile$Name <- sub('_N_0_', '_N_', nwheaderfile$Name)
    
    # Find all files with gps.txt extension:
    gpsfilenames <- paste(rootf, sfolders[i], sdates[j],
                          flist[grep('gps.txt', flist)],sep='/')
    # Check if there is more than one:
    if (length(gpsfilenames)>1) {
      # If there is, start with an empty dataframe
      gpsfile <- as.data.frame(NULL)
      # Loop through each gps file n:
      for (n in 1:length(gpsfilenames)) {
        # Append its data to gpsfile dataframe
        gpsfile <- rbind(gpsfile, read.csv(gpsfilenames[n],sep='\t',header=T))
      }
    } else {
      # If not, read gps file:
      gpsfile <- read.csv(gpsfilenames,sep='\t',header=T)  
    }
    # Vectorize and trim LATITUDE and LONGITUDE columsn in gps data:
    gpsfile$LATITUDE <- as.vector(gpsfile$LATITUDE)
    gpsfile$LONGITUDE <- as.vector(gpsfile$LONGITUDE)
    gpsfile$LATITUDE <- trim(gpsfile$LATITUDE)
    gpsfile$LONGITUDE <- trim(gpsfile$LONGITUDE)
    # Take of N suffix from LAT:
    gpsfile$LATITUDE <- sub(' N','',gpsfile$LATITUDE)
    # Check which LON is suffixed 'W':
    minus <- grep('W', gpsfile$LONGITUDE)
    # Remove LON suffixes:
    gpsfile$LONGITUDE <- sub(' W','',gpsfile$LONGITUDE)
    gpsfile$LONGITUDE <- sub(' E','',gpsfile$LONGITUDE)
    # Add '-' to those LONs marked W:
    gpsfile$LONGITUDE[minus] <- paste('-',gpsfile$LONGITUDE[minus],sep='')
    # Project coordinates in WGS84.
    gpsfile$LATITUDE <- as.numeric(gpsfile$LATITUDE)
    gpsfile$LONGITUDE <- as.numeric(gpsfile$LONGITUDE)
    coordinates(gpsfile) <- c('LONGITUDE','LATITUDE')
    proj4string(gpsfile) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")   ### WGS84
    gpsfile_bng <- spTransform(gpsfile, CRS("+init=epsg:27700"))
    gpsfile_bng <- as.data.frame(gpsfile_bng)
    gpsfile$BNG_x <- gpsfile_bng$LONGITUDE
    gpsfile$BNG_y <- gpsfile_bng$LATITUDE
    
    # Vectorise gps file name (rows) and remove extension so to match names with header file names:
    gpsfile$NAME <- as.vector(gpsfile$NAME)
    gpsfile$NAME <- sub('.wav', '',gpsfile$NAME)
    # Match LAT and LON columns from gps file to header file data, on both Name columns:
    nwheaderfile$LAT <- gpsfile$LATITUDE[match(nwheaderfile$Name, gpsfile$NAME)]
    nwheaderfile$LON <- gpsfile$LONGITUDE[match(nwheaderfile$Name, gpsfile$NAME)]
    nwheaderfile$BNG_y <- gpsfile$BNG_y[match(nwheaderfile$Name, gpsfile$NAME)]
    nwheaderfile$BNG_x <- gpsfile$BNG_x[match(nwheaderfile$Name, gpsfile$NAME)]
    
    # Sort species column.
    # First vectorise and trim the column.
    nwheaderfile$Species <- as.vector(nwheaderfile$Species)
    nwheaderfile$Species <- trim(nwheaderfile$Species)
    # Find those values where more than one species was labelled (value split by commas):
    no_mult <- grep(',',nwheaderfile$Species)
    # Find species for those values with more than one:
    spp_mult <- nwheaderfile[no_mult,'Species']
    # Split resulting values by commas:
    spp_mult_split <- strsplit(spp_mult, ',')
    # Count the number of species occurrences in each of these:
    spp_mult_count <- unlist(lapply(spp_mult_split, length))
    # Now repeat original header data X times for each value with X species based on values extracted above:
    headersection <- nwheaderfile[rep(no_mult, spp_mult_count),]
    # Make a list of the species values for each multi-spp occurrence (should be same length as new header
    #  section above), and change species column in new header section to these.
    # End result should be a species column with a single species per occurrence.
    headersection$Species <- unlist(spp_mult_split)
    # Remove the multi-species occurrences from original header file:
    if(nrow(nwheaderfile[-no_mult,])!=0) {
      nwheaderfile <- nwheaderfile[-no_mult,]
      # Add new header section (with multi-ssp occurrences now as repeated rows) to new header data: 
      nwheaderfile <- rbind(nwheaderfile, headersection)
      # Re-order new header data in order of occurrence:
      nwheaderfile <- nwheaderfile[order(as.numeric(row.names(nwheaderfile))),]
    }
    
    # Add column with site/survey id
    nwheaderfile$Loc <- trim(as.vector(nwheaderfile$Loc))
    nwheaderfile$SURVEYID <- paste(nwheaderfile$Loc, j, sep='-')
    
    # Write the new header file data with coords to the site/date folder:
    write.csv(nwheaderfile, paste(rootf, sfolders[i], sdates[j],
                                  'header_coord.csv', sep='/'),
              quote=which(names(nwheaderfile)=='LAT'|names(nwheaderfile)=='LON'))
    # Print site/date name to show progress:
    print(paste(sfolders[i], sdates[j], sep='/'))
  
    # Merge all bat fixes names with all noise fix names, and add coordinates, as reference.
    noisefiles <- list.files(paste(rootf, sfolders[i], sdates[j], 'NOISE',sep='/'))
    noisefiles <- sub('_000.00#', '', noisefiles)
    noisefiles <- sub('___0_', '___', noisefiles)
    allfiles <- c(nwheaderfile$Name, noisefiles)
    allfiles <- data.frame(Name=allfiles)
    allfiles$LAT <- gpsfile$LATITUDE[match(allfiles$Name, gpsfile$NAME)]
    allfiles$LON <- gpsfile$LONGITUDE[match(allfiles$Name, gpsfile$NAME)]
    allfiles$BNG_y <- gpsfile$BNG_y[match(allfiles$Name, gpsfile$NAME)]
    allfiles$BNG_x <- gpsfile$BNG_x[match(allfiles$Name, gpsfile$NAME)]
    allfiles$SURVEYID <- paste(sfolders[i], j, sep='-')
    write.csv(allfiles, paste(rootf, sfolders[i], sdates[j],
                                  'bats and noise fixes.csv', sep='/'),
              quote=which(names(allfiles)=='LAT'|names(allfiles)=='LON'))
    
    # Add current site/date data to 'all data' output:
    alldat <- rbind(alldat, nwheaderfile)
    alldat_inc_noise <- rbind(alldat_inc_noise, allfiles)
  }
  
} # Repeat above for all sites i and dates j.

# Write all site/date file to output folder.
write.csv(alldat, 'data/SWT 2014 all bat fixes with coords.csv', row.names=T)
write.csv(alldat_inc_noise, 'data/SWT 2014 all fixes INC NOISE with coords.csv', row.names=T)
