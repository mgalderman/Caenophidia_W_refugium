# Stastical testing for the refugium hypothesis
# Megan G. Alderman

# Goal: test the hypothesis that sex chromosomes, and the W in particular, have an excess of repeat elements (and ERVs/mdg4s, specifically). 

### Clear environment-------------------------------------------------------

rm(list = ls())

### chi-square tests of uniform distribution of repeats---------------------

# Here, we're testing if the observed bp of repeats in whichever category on autosomes, Z, and W chromosomes matches a uniform distribution, given the lengths of the chromosomes.

# The observed and expected values have been calculated in the Excel table. Plug in real numbers accordingly!

---- #Thamnophis elegans
## set expected and observed values
tot.obs.a <- 782502557
tot.obs.z <- 84690706
tot.obs.w <- 56966401

erv.obs.a <- 1192387
erv.obs.z <- 322586
erv.obs.w <- 718639

gyp.obs.a <- 8803809
gyp.obs.z <- 2116735
gyp.obs.w <- 2630790

l1.obs.a <- 16373401
l1.obs.z <- 3294263
l1.obs.w <- 5363574

# Proportion of genome made up by chrom type
prop.a <- 0.8754
prop.z <- 0.0847
prop.w <- 0.0399

## run chi-square tests

# all repeats
tot.obsfreq <- c(tot.obs.a,tot.obs.z,tot.obs.w)
tot.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(tot.obsfreq,p=tot.nullprobs)

# ERVs
erv.obsfreq <- c(erv.obs.a,erv.obs.z,erv.obs.w)
erv.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(erv.obsfreq,p=erv.nullprobs)

# mdg4 elements
gyp.obsfreq <- c(gyp.obs.a,gyp.obs.z,gyp.obs.w)
gyp.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(gyp.obsfreq,p=gyp.nullprobs)

# L1 elements
l1.obsfreq <- c(l1.obs.a,l1.obs.z,l1.obs.w)
l1.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(l1.obsfreq,p=l1.nullprobs)

---- #Naja naja
## set expected and observed values
tot.obs.a <- 755979110
tot.obs.z <- 82667624
tot.obs.w <- 28211548

erv.obs.a <- 629744
erv.obs.z <- 206029
erv.obs.w <- 240666

gyp.obs.a <- 33683780
gyp.obs.z <- 4030749
gyp.obs.w <- 1521924

l1.obs.a <- 17545570
l1.obs.z <- 4533959
l1.obs.w <- 3299726

# Proportion of genome made up by chrom type
prop.a <- 0.885
prop.z <- 0.090
prop.w <- 0.025

## run chi-square tests

# all repeats
tot.obsfreq <- c(tot.obs.a,tot.obs.z,tot.obs.w)
tot.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(tot.obsfreq,p=tot.nullprobs)

# ERVs
erv.obsfreq <- c(erv.obs.a,erv.obs.z,erv.obs.w)
erv.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(erv.obsfreq,p=erv.nullprobs)

# mdg4 elements
gyp.obsfreq <- c(gyp.obs.a,gyp.obs.z,gyp.obs.w)
gyp.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(gyp.obsfreq,p=gyp.nullprobs)

# L1 elements
l1.obsfreq <- c(l1.obs.a,l1.obs.z,l1.obs.w)
l1.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(l1.obsfreq,p=l1.nullprobs)

---- #Cerastes gasperettii
## set expected and observed values
tot.obs.a <- 754261561
tot.obs.z <- 82114527
tot.obs.w <- 52131096

erv.obs.a <- 1070300
erv.obs.z <- 170416
erv.obs.w <- 433959

gyp.obs.a <- 33683780
gyp.obs.z <- 3733210
gyp.obs.w <- 1872776

l1.obs.a <- 37272790
l1.obs.z <- 5602883
l1.obs.w <- 4099824

# Proportion of genome made up by chrom type
prop.a <- 0.883
prop.z <- 0.083
prop.w <- 0.034

## run chi-square tests

# all repeats
tot.obsfreq <- c(tot.obs.a,tot.obs.z,tot.obs.w)
tot.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(tot.obsfreq,p=tot.nullprobs)

# ERVs
erv.obsfreq <- c(erv.obs.a,erv.obs.z,erv.obs.w)
erv.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(erv.obsfreq,p=erv.nullprobs)

# mdg4 elements
gyp.obsfreq <- c(gyp.obs.a,gyp.obs.z,gyp.obs.w)
gyp.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(gyp.obsfreq,p=gyp.nullprobs)

# L1 elements
l1.obsfreq <- c(l1.obs.a,l1.obs.z,l1.obs.w)
l1.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(l1.obsfreq,p=l1.nullprobs)


---- #Vipera ursinii
## set expected and observed values
tot.obs.a <- 713812980
tot.obs.z <- 73120738
tot.obs.w <- 37478671

erv.obs.a <- 565645
erv.obs.z <- 87383
erv.obs.w <- 447747

gyp.obs.a <- 24284187
gyp.obs.z <- 3738924
gyp.obs.w <- 3260410

l1.obs.a <- 20312967
l1.obs.z <- 3048121
l1.obs.w <- 3391134

# Proportion of genome made up by chrom type
prop.a <- 0.892
prop.z <- 0.079
prop.w <- 0.029

## run chi-square tests

# all repeats
tot.obsfreq <- c(tot.obs.a,tot.obs.z,tot.obs.w)
tot.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(tot.obsfreq,p=tot.nullprobs)

# ERVs
erv.obsfreq <- c(erv.obs.a,erv.obs.z,erv.obs.w)
erv.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(erv.obsfreq,p=erv.nullprobs)

# mdg4 elements
gyp.obsfreq <- c(gyp.obs.a,gyp.obs.z,gyp.obs.w)
gyp.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(gyp.obsfreq,p=gyp.nullprobs)

# L1 elements
l1.obsfreq <- c(l1.obs.a,l1.obs.z,l1.obs.w)
l1.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(l1.obsfreq,p=l1.nullprobs)


---- #Vipera berus
## set expected and observed values
tot.obs.a <- 775505233
tot.obs.z <- 80999203
tot.obs.w <- 64042382

erv.obs.a <- 587564
erv.obs.z <- 115406
erv.obs.w <- 464245

gyp.obs.a <- 25003122
gyp.obs.z <- 3819065
gyp.obs.w <- 3581522

l1.obs.a <- 20772273
l1.obs.z <- 3202902
l1.obs.w <- 3963175

# Proportion of genome made up by chrom type
prop.a <- 0.874
prop.z <- 0.080
prop.w <- 0.046

## run chi-square tests

# all repeats
tot.obsfreq <- c(tot.obs.a,tot.obs.z,tot.obs.w)
tot.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(tot.obsfreq,p=tot.nullprobs)

# ERVs
erv.obsfreq <- c(erv.obs.a,erv.obs.z,erv.obs.w)
erv.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(erv.obsfreq,p=erv.nullprobs)

# mdg4 elements
gyp.obsfreq <- c(gyp.obs.a,gyp.obs.z,gyp.obs.w)
gyp.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(gyp.obsfreq,p=gyp.nullprobs)

# L1 elements
l1.obsfreq <- c(l1.obs.a,l1.obs.z,l1.obs.w)
l1.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(l1.obsfreq,p=l1.nullprobs)

---- #Deinagkistrodon acutus
## set expected and observed values
tot.obs.a <- 576107532
tot.obs.z <- 49363163
tot.obs.w <- 20893494

erv.obs.a <- 592776
erv.obs.z <- 147830
erv.obs.w <- 325855

gyp.obs.a <- 8488869
gyp.obs.z <- 2031139
gyp.obs.w <- 2211497

l1.obs.a <- 10262491
l1.obs.z <- 1502897
l1.obs.w <- 1240761

# Proportion of genome made up by chrom type
prop.a <- 0.912
prop.z <- 0.066
prop.w <- 0.022

## run chi-square tests

# all repeats
tot.obsfreq <- c(tot.obs.a,tot.obs.z,tot.obs.w)
tot.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(tot.obsfreq,p=tot.nullprobs)

# ERVs
erv.obsfreq <- c(erv.obs.a,erv.obs.z,erv.obs.w)
erv.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(erv.obsfreq,p=erv.nullprobs)

# mdg4 elements
gyp.obsfreq <- c(gyp.obs.a,gyp.obs.z,gyp.obs.w)
gyp.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(gyp.obsfreq,p=gyp.nullprobs)

# L1 elements
l1.obsfreq <- c(l1.obs.a,l1.obs.z,l1.obs.w)
l1.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(l1.obsfreq,p=l1.nullprobs)

---- #Crotalus horridus
## set expected and observed values
tot.obs.a <- 740027612
tot.obs.z <- 80277635
tot.obs.w <- 58096140

erv.obs.a <- 1573351
erv.obs.z <- 339165
erv.obs.w <- 1177635

gyp.obs.a <- 20406952
gyp.obs.z <- 3997214
gyp.obs.w <- 4504660

l1.obs.a <- 20594992
l1.obs.z <- 3637106
l1.obs.w <- 4744944

# Proportion of genome made up by chrom type
prop.a <- 0.872
prop.z <- 0.085
prop.w <- 0.043

## run chi-square tests

# all repeats
tot.obsfreq <- c(tot.obs.a,tot.obs.z,tot.obs.w)
tot.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(tot.obsfreq,p=tot.nullprobs)

# ERVs
erv.obsfreq <- c(erv.obs.a,erv.obs.z,erv.obs.w)
erv.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(erv.obsfreq,p=erv.nullprobs)

# mdg4 elements
gyp.obsfreq <- c(gyp.obs.a,gyp.obs.z,gyp.obs.w)
gyp.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(gyp.obsfreq,p=gyp.nullprobs)

# L1 elements
l1.obsfreq <- c(l1.obs.a,l1.obs.z,l1.obs.w)
l1.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(l1.obsfreq,p=l1.nullprobs)

---- #Crotalus adamanteus
## set expected and observed values
tot.obs.a <- 802609365
tot.obs.z <- 90697590
tot.obs.w <- 59241196

erv.obs.a <- 1767615
erv.obs.z <- 362595
erv.obs.w <- 1078079

gyp.obs.a <- 25945723
gyp.obs.z <- 4543369
gyp.obs.w <- 4238810

l1.obs.a <- 27116794
l1.obs.z <- 4493387
l1.obs.w <- 5015008

# Proportion of genome made up by chrom type
prop.a <- 0.875
prop.z <- 0.086
prop.w <- 0.039

## run chi-square tests

# all repeats
tot.obsfreq <- c(tot.obs.a,tot.obs.z,tot.obs.w)
tot.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(tot.obsfreq,p=tot.nullprobs)

# ERVs
erv.obsfreq <- c(erv.obs.a,erv.obs.z,erv.obs.w)
erv.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(erv.obsfreq,p=erv.nullprobs)

# mdg4 elements
gyp.obsfreq <- c(gyp.obs.a,gyp.obs.z,gyp.obs.w)
gyp.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(gyp.obsfreq,p=gyp.nullprobs)

# L1 elements
l1.obsfreq <- c(l1.obs.a,l1.obs.z,l1.obs.w)
l1.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(l1.obsfreq,p=l1.nullprobs)


---- #Crotalus viridis
## set expected and observed values
tot.obs.a <- 485384979
tot.obs.z <- 58313393
tot.obs.w <- 20115701

erv.obs.a <- 722000
erv.obs.z <- 196846
erv.obs.w <- 441115

gyp.obs.a <- 8982568
gyp.obs.z <- 2262576
gyp.obs.w <- 2319344

l1.obs.a <- 9724408
l1.obs.z <- 1975148
l1.obs.w <- 1348144

# Proportion of genome made up by chrom type
prop.a <- 0.899
prop.z <- 0.084
prop.w <- 0.017

## run chi-square tests

# all repeats
tot.obsfreq <- c(tot.obs.a,tot.obs.z,tot.obs.w)
tot.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(tot.obsfreq,p=tot.nullprobs)

# ERVs
erv.obsfreq <- c(erv.obs.a,erv.obs.z,erv.obs.w)
erv.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(erv.obsfreq,p=erv.nullprobs)

# mdg4 elements
gyp.obsfreq <- c(gyp.obs.a,gyp.obs.z,gyp.obs.w)
gyp.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(gyp.obsfreq,p=gyp.nullprobs)

# L1 elements
l1.obsfreq <- c(l1.obs.a,l1.obs.z,l1.obs.w)
l1.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(l1.obsfreq,p=l1.nullprobs)

