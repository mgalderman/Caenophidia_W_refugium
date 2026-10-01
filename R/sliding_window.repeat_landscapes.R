############################################################################
# Repeat element landscape summaries
############################################################################

### Goal: The purpose of this analysis is to create a sliding window summary of 
### repeat content for the Z and W chromosome of a subset of species to 
### depict local variation of repeat accumulation

### Clear environment-------------------------------------------------------

rm(list = ls())

### Set working directory --------------------------------------------------

setwd("/Users/meganalderman/Library/CloudStorage/GoogleDrive-schieldlab@gmail.com/My Drive/projects/caenophidia_W_refugium/analysis/repeatmasker_landscape/repeat_windows/results")

### Load dependencies-------------------------------------------------------

library(readxl)

### Read in data------------------------------------------------------------

berus_W <- read.table("v_berus_chrW_100kb.txt", header = T)
berus_Z <- read.table("v_berus_chrZ_100kb.txt", header = T)

ursinii_W <- read.table("v_ursinii_chrW_100kb.txt", header = T)
ursinii_Z <- read.table("v_ursinii_chrZ_100kb.txt", header = T)

gasperettii_W <- read.table("c_gasperettii_chrW_100kb.txt", header = T)
gasperettii_Z <- read.table("c_gasperettii_chrZ_100kb.txt", header = T)

adamanteus_W <- read.table("c_adamanteus_chrW_100kb.txt", header = T)
adamanteus_Z <- read.table("c_adamanteus_chrZ_100kb.txt", header = T)

par(mfrow=c(2,1))

plot(berus_Z$start, berus_Z$prop_repeats, type = "l", lwd=1.5, xlim = c(134003728,0), ylim = c(0,1.0))
plot(berus_W$start, berus_W$prop_repeats, type = "l", lwd=1.5, xlim = c(0,134003728), ylim = c(0,1.0))
abline(h = .53, col = "red", lty = 2, lwd = 1) # mean auto repeats
abline(h = .60, col = "blue", lty = 2, lwd = 1) # mean Z repeats

plot(ursinii_Z$start, ursinii_Z$prop_repeats, type = "l", lwd=1.5, xlim = c(0,127138603), ylim = c(0, 1))
plot(ursinii_W$start, ursinii_W$prop_repeats, type = "l", lwd=1.5, xlim = c(0,127138603), ylim = c(0, 1))
abline(h = .50, col = "red", lty = 2, lwd = 1) # mean auto repeats
abline(h = .58, col = "blue", lty = 2, lwd = 1) # mean Z repeats

plot(gasperettii_Z$start, gasperettii_Z$prop_repeats, type = "l", lwd=1.5, xlim = c(135322495,0), ylim = c(0, 1))
plot(gasperettii_W$start, gasperettii_W$prop_repeats, type = "l", lwd=1.5, xlim = c(0,135322495), ylim = c(0, 1))
abline(h = .53, col = "red", lty = 2, lwd = 1) # mean auto repeats
abline(h = .61, col = "blue", lty = 2, lwd = 1) # mean Z repeats

plot(adamanteus_Z$start, adamanteus_Z$prop_repeats, type = "l", lwd=1.5, xlim = c(0,145581062), ylim = c(0, 1))
plot(adamanteus_W$start, adamanteus_W$prop_repeats, type = "l", lwd=1.5, xlim = c(0,145581062), ylim = c(0, 1))
abline(h = .54, col = "red", lty = 2, lwd = 1) # mean auto repeats
abline(h = .62, col = "blue", lty = 2, lwd = 1) # mean Z repeats

# Export 16 x 8
