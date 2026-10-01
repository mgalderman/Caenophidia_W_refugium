## Dysregulation hypothesis 
# 5/21/26 by Megan Alderman
# Tests of the 'Dysregulation hypothesis' testing framework - 'heterochromatin sink' or 'silencing escape' models

library(tidyverse)
library(readxl)

setwd("~/Library/CloudStorage/GoogleDrive-schieldlab@gmail.com/My Drive/projects/caenophidia_W_refugium/analysis/dysregulation_hypothesis/featureCounts_output")

# Read in  gene-level counts -----

## Thamnophis elegans ---
raw_counts <- read_excel('t_elegans_featureCounts.xlsx') # read in file and ignore comment lines

raw_counts <- raw_counts %>% 
  mutate(length_female = case_when(
    str_detect(Geneid,'_auto') ~ Length * 2, # Make a column of length values to use for female sample (autos length x 2)
    .default= Length # for all non-autosome elements, populate this new column with the existing length value
  )) %>% 
  mutate(length_male = case_when(
    str_detect(Geneid,'_auto') ~ Length * 2, # Make a column of length values to use for male sample (autos and Z length x 2)
    str_detect(Geneid,'_chrZ') ~ Length * 2,
    .default= Length
  )) %>% 
  mutate(female_adj = t_elegans_F1_Aligned.sortedByCoord.out.bam / length_female) %>%  # Normalize female counts by female length
  mutate(male_adj = t_elegans_M1_Aligned.sortedByCoord.out.bam / length_male) # Normalize male counts by male length

raw_counts.mat <- as.matrix(raw_counts[,c('female_adj','male_adj')]) # makes a matrix of all of the count data (selecting just the count columns - 3 and 4)

row.names(raw_counts.mat) <- raw_counts$Geneid # sets the row IDs for the matrix to gene IDs from our initial table

# TPM Calculation
x <- raw_counts.mat # make temporary matrix (x) with raw counts normalized by sex-specific TE length above
tpm.mat <- t( t(x) * 1e6 / colSums(x) ) # scale each column (library) to be out of 1 million, generating a matrix of TPM normalized counts

colSums(tpm.mat) # if we check the total number of normalized counts per sample, we now see that all are scaled to 1e6

tpm.df <- as.data.frame(tpm.mat) # convert TPM matrix to dataframe if desired

write.csv(tpm.df, "t_elegans_TPM_new.csv")

# Separate into chromosome categories

TE_TPM_auto <- tpm.df %>% 
  filter(str_detect(row.names(tpm.df),'_auto')) %>% 
  select(Female=female_adj,Male=male_adj)
TE_TPM_chrZ <- tpm.df %>% 
  filter(str_detect(row.names(tpm.df),'_chrZ')) %>% 
  select(Female=female_adj,Male=male_adj)
TE_TPM_chrW <- tpm.df %>% 
  filter(str_detect(row.names(tpm.df),'_chrW')) %>% 
  select(Female=female_adj,Male=male_adj)

TE_TPM <-tpm.df %>%
select(Female=female_adj,Male=male_adj)

# Plot boxplots of the distributions of reads

TE_log_FA <- log2(TE_TPM_auto$Female) 
TE_log_FA <- TE_log_FA[is.finite((TE_log_FA))]

TE_log_FZ <- log2(TE_TPM_chrZ$Female)
TE_log_FZ <- TE_log_FZ[is.finite((TE_log_FZ))]

TE_log_FW <- log2(TE_TPM_chrW$Female)
TE_log_FW <- TE_log_FW[is.finite((TE_log_FW))]

TE_log_MA <- log2(TE_TPM_auto$Male)
TE_log_MA <- TE_log_MA[is.finite((TE_log_MA))]

TE_log_MZ <- log2(TE_TPM_chrZ$Male)
TE_log_MZ <- TE_log_MZ[is.finite((TE_log_MZ))]

TE_log_F_total <- log2(TE_TPM$Female)
TE_log_F_total <- TE_log_F_total[is.finite((TE_log_F_total))]

TE_log_M_total <- log2(TE_TPM$Male)
TE_log_M_total <- TE_log_M_total[is.finite((TE_log_M_total))]

TE_auto_Z <- rbind(TE_TPM_auto, TE_TPM_chrZ) # Combine auto and Z tables into one

TE_log_FAZ <- log2(TE_auto_Z$Female)
TE_log_FAZ <- TE_log_FAZ[is.finite((TE_log_FAZ))]

TE_log_MAZ <-log2(TE_auto_Z$Male)
TE_log_MAZ <- TE_log_MAZ[is.finite((TE_log_MAZ))]

# Make boxplots of distributions 

par(mfrow=c(1,5))
boxplot(TE_log_FA, TE_log_MA,names=c('Female auto','Male auto'), ylab='Log2 (TPM)',main='Thamnophis elegans')

boxplot(TE_log_FZ, TE_log_FW,names=c('Female Z', 'Female W'),ylab='Log2 (TPM)')

boxplot(TE_log_FZ, TE_log_MZ,names=c('Female Z','Male Z'),ylab='Log2 (TPM)')

boxplot(TE_log_F_total,TE_log_M_total,names=c('Female total','Male total'),ylab='Log2 (TPM)')

boxplot(TE_log_FAZ, TE_log_MAZ,names=c('Female auto+Z','Male auto+Z'), ylab='Log2 (TPM)')


## Cerastes gasperettii ---

raw_counts <- read_excel('c_gasperettii_featureCounts.xlsx') # read in file and ignore comment lines

raw_counts <- raw_counts %>% 
  mutate(length_female = case_when(
    str_detect(Geneid,'_auto') ~ Length * 2, # Make a column of length values to use for female sample (autos length x 2)
    .default= Length # for all non-autosome elements, populate this new column with the existing length value
  )) %>% 
  mutate(length_male = case_when(
    str_detect(Geneid,'_auto') ~ Length * 2, # Make a column of length values to use for male sample (autos and Z length x 2)
    str_detect(Geneid,'_chrZ') ~ Length * 2,
    .default= Length
  )) %>% 
  mutate(female_adj = Female / length_female) %>%  # Normalize female counts by female length
  mutate(male_adj = Male / length_male) # Normalize male counts by male length

raw_counts.mat <- as.matrix(raw_counts[,c('female_adj','male_adj')]) # makes a matrix of all of the count data (selecting just the count columns - 3 and 4)

row.names(raw_counts.mat) <- raw_counts$Geneid # sets the row IDs for the matrix to gene IDs from our initial table

# TPM Calculation
x <- raw_counts.mat # make temporary matrix (x) with raw counts normalized by sex-specific TE length above

tpm.mat <- t( t(x) * 1e6 / colSums(x) ) # scale each column (library) to be out of 1 million, generating a matrix of TPM normalized counts

colSums(tpm.mat) # if we check the total number of normalized counts per sample, we now see that all are scaled to 1e6

tpm.df <- as.data.frame(tpm.mat) # convert TPM matrix to dataframe if desired

write.csv(tpm.df, "c_gasperettii_TPM_new.csv")

CG_TPM_auto <- tpm.df %>% 
  filter(str_detect(row.names(tpm.df),'_auto')) %>% 
  select(Female=female_adj,Male=male_adj)
CG_TPM_chrZ <- tpm.df %>% 
  filter(str_detect(row.names(tpm.df),'_chrZ')) %>% 
  select(Female=female_adj,Male=male_adj)
CG_TPM_chrW <- tpm.df %>% 
  filter(str_detect(row.names(tpm.df),'_chrW')) %>% 
  select(Female=female_adj,Male=male_adj)

CG_TPM <-tpm.df %>%
  select(Female=female_adj,Male=male_adj)

# Plot boxplots of the distributions of reads

CG_log_FA <- log2(CG_TPM_auto$Female) 
CG_log_FA <- CG_log_FA[is.finite((CG_log_FA))]

CG_log_FZ <- log2(CG_TPM_chrZ$Female)
CG_log_FZ <- CG_log_FZ[is.finite((CG_log_FZ))]

CG_log_FW <- log2(CG_TPM_chrW$Female)
CG_log_FW <- CG_log_FW[is.finite((CG_log_FW))]

CG_log_MA <- log2(CG_TPM_auto$Male)
CG_log_MA <- CG_log_MA[is.finite((CG_log_MA))]

CG_log_MZ <- log2(CG_TPM_chrZ$Male)
CG_log_MZ <- CG_log_MZ[is.finite((CG_log_MZ))]

CG_log_F_total <- log2(CG_TPM$Female)
CG_log_F_total <- CG_log_F_total[is.finite((CG_log_F_total))]

CG_log_M_total <- log2(CG_TPM$Male)
CG_log_M_total <- CG_log_M_total[is.finite((CG_log_M_total))]

CG_auto_Z <- rbind(CG_TPM_auto, CG_TPM_chrZ)

CG_FAZ <- sum(CG_auto_Z$Female)
CG_MAZ <- sum(CG_auto_Z$Male)

CG_log_FAZ <- log2(CG_auto_Z$Female)
CG_log_FAZ <- CG_log_FAZ[is.finite((CG_log_FAZ))]

CG_log_MAZ <-log2(CG_auto_Z$Male)
CG_log_MAZ <- CG_log_MAZ[is.finite((CG_log_MAZ))]

par(mfrow=c(1,5),oma = c(0, 0, 3, 0))
boxplot(CG_log_FA, CG_log_MA,names=c('Female auto','Male auto'), ylab='Log2 (TPM)',main='Cerastes gasperettii')

boxplot(CG_log_FZ, CG_log_FW,names=c('Female Z','Female W'),ylab='Log2 (TPM)')

boxplot(CG_log_FZ, CG_log_MZ,names=c('Female Z','Male Z'),ylab='Log2 (TPM)')

boxplot(CG_log_F_total,CG_log_M_total,names=c('Female total','Male total'),ylab='Log2 (TPM)')

boxplot(CG_log_FAZ, CG_log_MAZ,names=c('Female auto+Z','Male auto+Z'), ylab='Log2 (TPM)')

## Crotalus adamanteus ---

raw_counts <- read_excel('c_adamanteus_featureCounts.xlsx') # read in file and ignore comment lines

raw_counts <- raw_counts %>% 
  mutate(length_female = case_when(
    str_detect(Geneid,'_auto') ~ Length * 2, # Make a column of length values to use for female sample (autos length x 2)
    .default= Length # for all non-autosome elements, populate this new column with the existing length value
  )) %>% 
  mutate(length_male = case_when(
    str_detect(Geneid,'_auto') ~ Length * 2, # Make a column of length values to use for male sample (autos and Z length x 2)
    str_detect(Geneid,'_chrZ') ~ Length * 2,
    .default= Length
  )) %>% 
  mutate(female_adj = Female / length_female) %>%  # Normalize female counts by female length
  mutate(male_adj = Male / length_male) # Normalize male counts by male length

raw_counts.mat <- as.matrix(raw_counts[,c('female_adj','male_adj')]) # makes a matrix of all of the count data (selecting just the count columns - 3 and 4)
row.names(raw_counts.mat) <- raw_counts$Geneid # sets the row IDs for the matrix to gene IDs from our initial table

# TPM Calculation
x <- raw_counts.mat # make temporary matrix (x) with raw counts normalized by sex-specific TE length above

tpm.mat <- t( t(x) * 1e6 / colSums(x) ) # scale each column (library) to be out of 1 million, generating a matrix of TPM normalized counts

colSums(tpm.mat) # if we check the total number of normalized counts per sample, we now see that all are scaled to 1e6

tpm.df <- as.data.frame(tpm.mat) # convert TPM matrix to dataframe if desired

# Separate out chromosome categories

CA_TPM_auto <- tpm.df %>% 
  filter(str_detect(row.names(tpm.df),'_auto')) %>% 
  select(Female=female_adj,Male=male_adj)
CA_TPM_chrZ <- tpm.df %>% 
  filter(str_detect(row.names(tpm.df),'_chrZ')) %>% 
  select(Female=female_adj,Male=male_adj)
CA_TPM_chrW <- tpm.df %>% 
  filter(str_detect(row.names(tpm.df),'_chrW')) %>% 
  select(Female=female_adj,Male=male_adj)

CA_TPM <-tpm.df %>%
  select(Female=female_adj,Male=male_adj)

# Plot boxplots of the distributions of reads

CA_log_FA <- log2(CA_TPM_auto$Female) 
CA_log_FA <- CA_log_FA[is.finite((CA_log_FA))]

CA_log_FZ <- log2(CA_TPM_chrZ$Female)
CA_log_FZ <- CA_log_FZ[is.finite((CA_log_FZ))]

CA_log_FW <- log2(CA_TPM_chrW$Female)
CA_log_FW <- CA_log_FW[is.finite((CA_log_FW))]

CA_log_MA <- log2(CA_TPM_auto$Male)
CA_log_MA <- CA_log_MA[is.finite((CA_log_MA))]

CA_log_MZ <- log2(CA_TPM_chrZ$Male)
CA_log_MZ <- CA_log_MZ[is.finite((CA_log_MZ))]

CA_log_F_total <- log2(CA_TPM$Female)
CA_log_F_total <- CA_log_F_total[is.finite((CA_log_F_total))]

CA_log_M_total <- log2(CA_TPM$Male)
CA_log_M_total <- CA_log_M_total[is.finite((CA_log_M_total))]

CA_auto_Z <- rbind(CA_TPM_auto, CA_TPM_chrZ)

CA_FAZ <- sum(CA_auto_Z$Female)
CA_MAZ <- sum(CA_auto_Z$Male)

CA_log_FAZ <- log2(CA_auto_Z$Female)
CA_log_FAZ <- CA_log_FAZ[is.finite((CA_log_FAZ))]

CA_log_MAZ <-log2(CA_auto_Z$Male)
CA_log_MAZ <- CA_log_MAZ[is.finite((CA_log_MAZ))]


par(mfrow=c(1,5),oma = c(0, 0, 3, 0))
boxplot(CA_log_FA, CA_log_MA,names=c('Female auto','Male auto'), ylab='Log2 (TPM)',main='Crotalus adamanteus')

boxplot(CA_log_FZ, CA_log_FW,names=c('Female Z', 'Female W'),ylab='Log2 (TPM)')

boxplot(CA_log_FZ, CA_log_MZ,names=c('Female Z','Male Z'),ylab='Log2 (TPM)')

boxplot(CA_log_F_total,CA_log_M_total,names=c('Female total','Male total'),ylab='Log2 (TPM)')

boxplot(CA_log_FAZ, CA_log_MAZ,names=c('Female auto+Z','Male auto+Z'), ylab='Log2 (TPM)')


# Export 12 x 6

### Perform t-tests ---------------------------------
#install.packages("multcomp")
#install.packages("FSA")
library(multcomp)
library(dplyr)
library(tidyr)
library(FSA)
library(car)

## T. elegans ---

# Male auto and female auto
t.test(TE_log_MA, TE_log_FA) #p-value = 0.127

# Female auto and female W
t.test(TE_log_FA, TE_log_FW) # p-value = 0.01837

# Male Z and female Z
t.test(TE_log_MZ, TE_log_FZ) # p-value = 0.7732

# Female Z and female W
t.test(TE_log_FZ, TE_log_FW) # p-value = 1.635e-11

# Male auto male ZZ
t.test(TE_log_MA, TE_log_MZ) # p-value = 1.806e-05

# Female auto male ZZ
t.test(TE_log_FA, TE_log_MZ) # p-value = 5.67e-08

# Female auto female Z
t.test(TE_log_FA, TE_log_FZ) # p-value = 3.16e-08


## C.gasperettii ---

# Male auto and female auto
t.test(CG_log_MA, CG_log_FA) # p-value = 0.2761

# Female auto and female W
t.test(CG_log_FA, CG_log_FW) # p-value = 0.4955

# Male Z and female Z
t.test(CG_log_MZ, CG_log_FZ) # p-value = 0.3856

# Female Z and female W
t.test(CG_log_FZ, CG_log_FW) # p-value = 0.02283

# Male auto male ZZ
t.test(CG_log_MA, CG_log_MZ) # p-value = 0.8335

# Female auto male ZZ
t.test(CG_log_FA, CG_log_MZ) # p-value = 0.4347

# Female auto female Z
t.test(CG_log_FA, CG_log_FZ) # p-value = 0.08621


## C.adamanteus ---

# Male auto and female auto
t.test(CA_log_MA, CA_log_FA) # p-value = 0.04539

# Female auto and female W
t.test(CA_log_FA, CA_log_FW) # p-value = 0.0142

# Male Z and female Z
t.test(CA_log_MZ, CA_log_FZ) # p-value = 0.6306

# Female Z and female W
t.test(CA_log_FZ, CA_log_FW) # p-value = 0.4427

# Male auto male ZZ
t.test(CA_log_MA, CA_log_MZ) # p-value = 0.08425

# Female auto male ZZ
t.test(CA_log_FA, CA_log_MZ) # p-value = 0.0008519

# Female auto female Z
t.test(CA_log_FA, CA_log_FZ) # p-value = 0.001721

## Plot ------
# Plot distributions of TE expression across chromosome types
# Boxplots of just female, within-sex comparisons ------

## Thamnophis elegans ---
par(mfrow=c(1,3))
thamnophis_dist <- boxplot(TE_log_FA, TE_log_FZ, TE_log_FW,
                           lwd = 0.75, las = 1,
                           col=c('#cccccc','#fcb017','#ffce77'),
                           names=c('Female auto', 'Female Z', 'Female W'),
                           whiskercol = "gray40",
                           ylab='log2 (TPM)', ylim=c(5,20), main='Thamnophis elegans')
#legend("topright", legend = c("p-value < 0.001"), bty = "n", cex=1.3)

t.test(TE_log_FA, TE_log_FZ) #p-value = 3.16e-08
t.test(TE_log_FA, TE_log_FW) #p-value = 0.01837
t.test(TE_log_FZ, TE_log_FW) #p-value = 1.635e-11

# significance between female auto and female Z
segments(x0 = 1, y0 = 16, x1 = 2, y1 = 16)
text(x = 1.5, y = 16.25, labels = "***", cex = 2)

# significance between female auto and female W
segments(x0 = 1, y0 = 17, x1 = 3, y1 = 17)
text(x = 2, y = 17.25, labels = "*", cex = 2)

# significance between female Z and female W
segments(x0 = 2, y0 = 18, x1 = 3, y1 = 18)
text(x = 2.5, y = 18.25, labels = "***", cex = 2)


cerastes_dist <- boxplot(CG_log_FA, CG_log_FZ, CG_log_FW,
                         lwd = 0.75, las = 1,
                         col=c('#cccccc','#fcb017','#ffce77'),
                         names=c('Female auto','Female Z', 'Female W'),
                         whiskercol = "gray40",
                         ylab='log2 (TPM)', ylim=c(5,20), main='Cerastes gasperettii')

t.test(CG_log_FA, CG_log_FZ) #p-value = 0.08621
t.test(CG_log_FA, CG_log_FW) #p-value = 0.4955
t.test(CG_log_FZ, CG_log_FW) #p-value = 0.02283

# significance between female auto and female Z
#segments(x0 = 1, y0 = 16, x1 = 2, y1 = 16)
#text(x = 1.5, y = 16.25, labels = "***", cex = 2)

# significance between female auto and female W
#segments(x0 = 1, y0 = 17, x1 = 3, y1 = 17)
#text(x = 2, y = 17.25, labels = "*", cex = 2)

# significance between female Z and female W
segments(x0 = 2, y0 = 18, x1 = 3, y1 = 18)
text(x = 2.5, y = 18.25, labels = "*", cex = 2)


adamanteus_dist <- boxplot(CA_log_FA, CA_log_FZ, CA_log_FW,
                           lwd = 0.75, las = 1,
                           col=c('#cccccc','#fcb017','#ffce77'),
                           names=c('Female auto','Female Z', 'Female W'),
                           whiskercol = "gray40",
                           ylab='log2 (TPM)', ylim=c(5,20), main='Crotalus adamanteus')

t.test(CA_log_FA, CA_log_FZ) #p-value = 0.001721
t.test(CA_log_FA, CA_log_FW) #p-value = 0.0142
t.test(CA_log_FZ, CA_log_FW) #p-value = 0.4427

# significance between female auto and female Z
segments(x0 = 1, y0 = 16, x1 = 2, y1 = 16)
text(x = 1.5, y = 16.25, labels = "**", cex = 2)

# significance between female auto and female W
segments(x0 = 1, y0 = 17, x1 = 3, y1 = 17)
text(x = 2, y = 17.25, labels = "*", cex = 2)

# Export at 1600 x 600
