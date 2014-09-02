library(gdata)
rootf <- '../2014/Bats'
sfolders <- list.files(rootf)
nfolders <- length(sfolders)

for (i in 1:nfolders) {
  sdates <- list.files(paste(rootf,sfolders[i],sep='/'))
  for (j in 1:length(sdates)) {
    flist <- list.files(paste(rootf, sfolders[i], sdates[j], sep='/'))
    headerfile <- read.csv(paste(rootf, sfolders[i], sdates[j], 'header.txt', sep='/'),header=T, sep='\t')
    #     headerfile$Species <- as.vector(headerfile$Species)
    #     headerfile$Species <- trim(headerfile$Species)
    #     no_mult <- grep(',',headerfile$Species)
    #     spp_mult <- headerfile[no_mult,'Species']
    #     spp_mult_split <- strsplit(spp_mult, ',')
    #     spp_mult_count <- unlist(lapply(spp_mult_split, length))
    #     headersection <- headerfile[rep(no_mult, spp_mult_count),]
    #     headersection$Species <- unlist(spp_mult_split)
    
    nwheaderfile$Name <- sub('_000.00#','',nwheaderfile$Name)
    nwheaderfile$Name <- sub('___0_', '___', nwheaderfile$Name)
    nwheaderfile$Name <- sub('_N_0_', '_N_', nwheaderfile$Name)
    
    gpsfilenames <- paste(rootf, sfolders[i], sdates[j],
                          flist[grep('gps.txt', flist)],sep='/')
    if (length(gpsfilenames)>1) {
      gpsfile <- as.data.frame(NULL)
      for (n in 1:length(gpsfilenames)) {
        gpsfile <- rbind(gpsfile, read.csv(gpsfilenames[n],sep='\t',header=T))
      }
    } else {
      gpsfile <- read.csv(gpsfilenames,sep='\t',header=T)  
    }
    gpsfile$LATITUDE <- as.vector(gpsfile$LATITUDE)
    gpsfile$LONGITUDE <- as.vector(gpsfile$LONGITUDE)
    gpsfile$LATITUDE <- trim(gpsfile$LATITUDE)
    gpsfile$LONGITUDE <- trim(gpsfile$LONGITUDE)
    gpsfile$LATITUDE <- sub(' N','',gpsfile$LATITUDE)
    minus <- grep('W', gpsfile$LONGITUDE)
    gpsfile$LONGITUDE <- sub(' W','',gpsfile$LONGITUDE)
    gpsfile$LONGITUDE <- sub(' E','',gpsfile$LONGITUDE)
    gpsfile$LONGITUDE[minus] <- paste('-',gpsfile$LONGITUDE[minus],sep='')
    gpsfile$NAME <- as.vector(gpsfile$NAME)
    gpsfile$NAME <- sub('.wav', '',gpsfile$NAME)
    
    nwheaderfile$LAT <- gpsfile$LATITUDE[match(nwheaderfile$Name, gpsfile$NAME)]
    nwheaderfile$LON <- gpsfile$LONGITUDE[match(nwheaderfile$Name, gpsfile$NAME)]
    
    write.csv(nwheaderfile, paste(rootf, sfolders[i], sdates[j],
                                  'header_coord.csv', sep='/'),
              quote=which(names(nwheaderfile)=='LAT'|names(nwheaderfile)=='LON'), row.names=F)
    print(paste(sfolders[i], sdates[j], sep='/'))
  }
}