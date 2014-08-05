library(xlsx)
library(lattice)
bats <- read.xlsx('data/Bat data per section.xlsx',sheetIndex=1,header=T)

bats$WINDS <- as.numeric(as.vector(bats$WINDS))
bats$fSECTION <- factor(bats$SECTION)

bats$pips <- bats$SOPPIP+bats$PIP+bats$COMPIP

histogram(~pips | SITE, data=bats)
bwplot(pips ~ fSECTION | SITE, data=bats)

