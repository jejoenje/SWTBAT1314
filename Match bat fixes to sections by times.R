library(xlsx)
library(gdata)

# Set root folder with data:
rootf <- '../2014/Bats'
# List folders in data folder - should ONLY be individual sites (check!)
sfolders <- list.files(rootf)
# Count no. folders (== no. of sites)
nfolders <- length(sfolders)

sect <- read.xlsx('data/Bat data per section 2014.xlsx', sheetIndex=1)
# Extract and correct date-time:
temp <- paste(sect$DATE, as.character(format(strptime(sect$s_dtime, format="%Y-%m-%d %H:%M:%S"),"%H:%M:%S")))
sect$s_dtime <- strptime(temp, format="%Y-%m-%d %H:%M:%S")

visits <- levels(sect$SURVEYID)
sect_out <- as.data.frame(NULL)
for (i in 1:length(visits)) {
  vst_out <- as.data.frame(NULL)
  vst <- sect[sect$SURVEYID==visits[i], ]
  transects <- levels(factor(vst$TRSCT))
  for (j in 1:length(transects)) {
    vst_t <- vst[vst$TRSCT==transects[j],]
    end <- vst_t$s_dtime[-1]
    end <- c(end, end[length(end)]+max(diff(vst_t$s_dtime)))
    vst_t$e_dtime <- end
    vst_out <- rbind(vst_out, vst_t)
  }
  sect_out <- rbind(sect_out, vst_out)
}
#sect_out$test <- as.numeric(sect_out$e_dtime-sect_out$s_dtime)
#write.csv(sect_out, 'data/Temp/sect_out.csv', row.names=F)

bats <- read.csv('data/SWT 2014 all bat fixes with coords.csv', header=T)
bats_dt <- paste(bats$Y, bats$M, bats$D, sep='-')
bats_t <- paste(bats$H, bats$M.1, bats$S, sep=':')
bats$dtime <- strptime(paste(bats_dt, bats_t), format='%Y-%m-%d %H:%M:%S')
bats$Species <- as.vector(bats$Species)
bats$Species <- sub('jm', '', bats$Species)
bats$Note <- trim(as.vector(bats$Note))
bats$Note[is.na(bats$Note)] <- ""
counts <- c('FLAG','UNK','SOPPIP','COMPIP','MYOTI','PIP','DUNNO','PAUR','NOCTU','NOISE')
sect_out[,counts] <- apply(sect_out[,counts],2,as.numeric)
sect_out[,counts] <- 0


for(i in 1:nrow(sect_out)) {
  s <- sect_out$s_dtime[i]
  e <- sect_out$e_dtime[i]
  bats_sect <- bats[bats$SURVEYID==sect_out$SURVEYID[i],]
  spp <- bats_sect[bats_sect$dtime>s & bats_sect$dtime<e,'Species']
  notes <- bats_sect[bats_sect$dtime>s & bats_sect$dtime<e,'Note']
  notes <- toupper(notes)
  
  for (j in 1:length(notes[notes!=''])) {
    
  }  
  if (notes!='') {
    sptable <- data.frame(spp=gsub("[ 0-9]","",trim(unlist(strsplit(notes[1],',')))),
                          no=as.numeric(gsub("[ A-Z]","",trim(unlist(strsplit(notes[1],','))))))
    sect_out[i,as.vector(sptable$spp)] <- sptable$no
  }

  
}


sect_out2 <- as.data.frame(NULL)
for(i in 1:nrow(bats)) {
  bats_sect <- sect_out[sect_out$SURVEYID==bats$SURVEYID[i],]
  if (sum((bats$dtime[i]>bats_sect$s_dtime)*(bats$dtime[i]<bats_sect$e_dtime))!=0) {
    intvl <- which(as.logical((bats$dtime[i]>bats_sect$s_dtime)*(bats$dtime[i]<bats_sect$e_dtime)))  
    if(bats$Species[i]!='SOCCAL' & bats$Species[i]!='' & bats$Species[i]!='BUZZ') {
      if(bats$Note[i]=='') {
        bats_sect[intvl,bats$Species[i]] <- bats_sect[intvl,bats$Species[i]] + 1
      } else {
        spp <- gsub("[ 0-9]","",unlist(lapply(strsplit(bats$Note[i],','),trim)))
        spp <- toupper(spp)
        no <- as.numeric(gsub("[ a-z]","",unlist(lapply(strsplit(bats$Note[i],','),trim))))
        bats_sect[intvl,spp] <- bats_sect[intvl,spp] + no
      }
    }
  }
  sect_out2 <- rbind(sect_out2, bats_sect[intvl,])
  print(i)
}