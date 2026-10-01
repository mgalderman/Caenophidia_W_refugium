# repeat_proportion_correlations.R

### Goal: test the hypothesis that the density of repeats on the W chromosome is predicted by values on autosomes and the Z chromosome.
### Addressing this is one way of determining whether the W has unique composition at a very broad level, motivating further detailed
### description of W-biased expansions of particular element families.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Format vectors with repeat proportion data across species ----

auto <- c(52.17,49.96,52.65,49.68,53.19,42.76,54.59,54.12,39.58)
chrz <- c(58.35,53.47,60.68,57.51,60.45,50.26,60.76,62.3,51.16)
chrw <- c(83.32,65.92,93.61,81.54,83.56,64.37,87,90.09,84.74)

props <- data.frame(auto,chrz,chrw)

### Perform correlation tests ----
cor.test (props$auto,props$chrz,method='spearman')
cor.test(props$auto,props$chrw,method='spearman')
cor.test(props$chrz,props$chrw,method='spearman')

### Plot ----
par(mfrow=c(1,3))
plot(props$auto,props$chrz,pch=20,xlab='Percent repeats (autosome)',ylab='Percent repeats (Z chromosome)',xlim=c(38,56),ylim=c(48,64))
abline(lm(props$chrz~props$auto))
plot(props$auto,props$chrw,pch=20,xlab='Percent repeats (autosome)',ylab='Percent repeats (W chromosome)',xlim=c(38,56),ylim=c(62,98))
abline(lm(props$chrw~props$auto))
plot(props$chrz,props$chrw,pch=20,xlab='Percent repeats (Z chromosome)',ylab='Percent repeats (W chromosome)',xlim=c(48,64),ylim=c(62,98))
abline(lm(props$chrw~props$chrz))

# Output at 3.25 x 8.75

### Read in repeat element proportion details for comparison between the Z and W ----

zw <- read.table('~/Library/CloudStorage/GoogleDrive-schieldlab@gmail.com/My Drive/projects/caenophidia_W_refugium/analysis/chromosome_repeat_density/chromosome_repeat_proportions_detail.txt',header=T)

### Plot ----
par(mfrow=c(1,5))
boxplot(zw$DNAhaT.Z,zw$DNAhaT.W,names=c('Z','W'),ylab='Percent of chromosome',main='DNA hAT')
boxplot(zw$CR1.Z,zw$CR1.W,names=c('Z','W'),ylab='Percent of chromosome',main='CR1 LINEs')
boxplot(zw$L1.Z,zw$L1.W,names=c('Z','W'),ylab='Percent of chromosome',main='L1 LINEs')
boxplot(zw$mdg4.Z,zw$mdg4.W,names=c('Z','W'),ylab='Percent of chromosome',main='mdg4 LTRs')
boxplot(zw$ERV1.Z,zw$ERV1.W,names=c('Z','W'),ylab='Percent of chromosome',main='ERV1')

# Output at 3.25 x 10

### Perform Mann-Whitney U tests ----

wilcox.test(zw$DNAhaT.Z,zw$DNAhaT.W)
wilcox.test(zw$CR1.Z,zw$CR1.W)
wilcox.test(zw$L1.Z,zw$L1.W)
wilcox.test(zw$mdg4.Z,zw$mdg4.W)
wilcox.test(zw$ERV1.Z,zw$ERV1.W)
