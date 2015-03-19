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
# sect_out$test <- as.numeric(sect_out$e_dtime-sect_out$s_dtime)
# write.csv(sect_out, 'data/Temp/sect_out.csv', row.names=F)

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

for (i in 1:nrow(sect_out)) {
  s <- sect_out$s_dtime[i]
  e <- sect_out$e_dtime[i]
  bats_surv <- bats[bats$SURVEYID==sect_out$SURVEYID[i],]
  bats_sect <- bats_surv[bats_surv$dtime>s & bats_surv$dtime<=e,]
  if(nrow(bats_sect)>0) {
    sptab <- data.frame(spp=bats_sect$Species, no=bats_sect$COUNT)
    sptab <- sptab[sptab$spp %in% counts,]
    sptab$spp <- factor(as.vector(sptab$spp))
    sptab <- tapply(sptab$no, sptab$spp, length)
    sect_out[i,names(sptab)] <- sptab
  }
  print(i)
}

write.csv(sect_out, '/data/Bat data per section processed.csv', row.names=F)
