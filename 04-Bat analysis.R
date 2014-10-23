###
###
### ANALYSE BAT DATA 2013+2014
### Avoid all edits to data file, please do in previous script.

library(lme4)

bats <- read.csv('data/BAT DATA 2013 and 2014 MASTER.csv',header=T)

mod1 <- glmer(PASSES ~ TTMIDN + (1|SITE/TRSCT/SECTION), family='poisson', data=bats)

library(glmmADMB)

bats$SITE <- factor(bats$SITE)
bats$TRSCT <- factor(bats$TRSCT)
bats$SECTION <- factor(bats$SECTION)


mod2 <- glmmadmb(PASSES ~ TTMIDN + (1|SITE/TRSCT/SECTION), family='nbinom1', data=bats)
