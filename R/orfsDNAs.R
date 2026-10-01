# Toxicity hypothesis testing - DNA transposons
# 3/10/26 by Megan Alderman

# Code modified from James Galbraith's script orfChecker.R https://github.com/jamesdgalbraith/OrthologueRegions/blob/master/orfChecker.R
# Code modified from Valentina Peona

# Install packages

install.packages("BSgenome")
install.packages("ORFik")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
BiocManager::install("BSgenome")
BiocManager::install("ORFik")

# import libraries
library(tidyverse)
library(GenomicRanges)
library(BSgenome)
library(ORFik)

# Set working directory
setwd("/Users/meganalderman/Library/CloudStorage/GoogleDrive-schieldlab@gmail.com/My Drive/projects/caenophidia_W_refugium/analysis/toxicity_hypothesis/DNA_transposons")

# Naming: { Genus species chrom} GSC ---- i.e.:
#TEA --> Thamnophis elegans autosome
#TEZ --> Thamnophis elegans chrZ
#TEW --> Thamnophis elegans chrW
--------------------------
  ### Thamnophis elegans ###
  # read in DNAtransp sequences
raw_seq_elegans_auto <- readDNAStringSet("t_elegans_DNAtransp.auto.fasta")
raw_seq_elegans_chrZ <- readDNAStringSet("t_elegans_DNAtransp.chrZ.fasta")
raw_seq_elegans_chrW <- readDNAStringSet("t_elegans_DNAtransp.chrW.fasta")

# create pseudoranges tibble of DNA elements
raw_tbl_elegans_auto <- tibble(seqnames = names(raw_seq_elegans_auto), start = 1, end = width(raw_seq_elegans_auto)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")
raw_tbl_elegans_chrZ <- tibble(seqnames = names(raw_seq_elegans_chrZ), start = 1, end = width(raw_seq_elegans_chrZ)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")
raw_tbl_elegans_chrW <- tibble(seqnames = names(raw_seq_elegans_chrW), start = 1, end = width(raw_seq_elegans_chrW)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")

# create table of DNAtransp seqs with names
seq_group_tbl_elegans_auto <- raw_tbl_elegans_auto %>% dplyr::select(seqnames)
seq_group_tbl_elegans_chrZ <- raw_tbl_elegans_chrZ %>% dplyr::select(seqnames)
seq_group_tbl_elegans_chrW <- raw_tbl_elegans_chrW %>% dplyr::select(seqnames)

# find orfs over 1000bp in sequences
#orfs_1000 <- ORFik::findORFs(raw_seq, startCodon = startDefinition(1), minimumLength = 1000) %>% as_tibble()

orfs_1000_TEA = findORFsFasta("t_elegans_DNAtransp.auto.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_TEA = as_tibble(orfs_1000_TEA)
orfs_1000_TEA$idv_name = paste0(orfs_1000_TEA$seqnames, "#orf", 1:nrow(orfs_1000_TEA))

orfs_1000_TEZ = findORFsFasta("t_elegans_DNAtransp.chrZ.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_TEZ = as_tibble(orfs_1000_TEZ)
orfs_1000_TEZ$idv_name = paste0(orfs_1000_TEZ$seqnames, "#orf", 1:nrow(orfs_1000_TEZ))

orfs_1000_TEW = findORFsFasta("t_elegans_DNAtransp.chrW.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_TEW = as_tibble(orfs_1000_TEW)
orfs_1000_TEW$idv_name = paste0(orfs_1000_TEW$seqnames, "#orf", 1:nrow(orfs_1000_TEW))


# make ranges object of orfs
orfs_1000_TEA_ranges <- GRanges(seqnames = orfs_1000_TEA$seqnames, ranges = IRanges(start = orfs_1000_TEA$start, end = orfs_1000_TEA$end))
orfs_1000_TEZ_ranges <- GRanges(seqnames = orfs_1000_TEZ$seqnames, ranges = IRanges(start = orfs_1000_TEZ$start, end = orfs_1000_TEZ$end))
orfs_1000_TEW_ranges <- GRanges(seqnames = orfs_1000_TEW$seqnames, ranges = IRanges(start = orfs_1000_TEW$start, end = orfs_1000_TEW$end))

# get seq of orfs and name orfs
orfs_seq_TEA <- Biostrings::getSeq(raw_seq_elegans_auto, orfs_1000_TEA_ranges)
names(orfs_seq_TEA) <- orfs_1000_TEA$idv_name

orfs_seq_TEZ <- Biostrings::getSeq(raw_seq_elegans_chrZ, orfs_1000_TEZ_ranges)
names(orfs_seq_TEZ) <- orfs_1000_TEZ$idv_name

orfs_seq_TEW <- Biostrings::getSeq(raw_seq_elegans_chrW, orfs_1000_TEW_ranges)
names(orfs_seq_TEW) <- orfs_1000_TEW$idv_name

# translate orf
orfs_aa_seq_TEA <- translate(orfs_seq_TEA, if.fuzzy.codon = "solve")
names(orfs_aa_seq_TEA) <- orfs_1000_TEA$idv_name

orfs_aa_seq_TEZ <- translate(orfs_seq_TEZ, if.fuzzy.codon = "solve")
names(orfs_aa_seq_TEZ) <- orfs_1000_TEZ$idv_name

orfs_aa_seq_TEW <- translate(orfs_seq_TEW, if.fuzzy.codon = "solve")
names(orfs_aa_seq_TEW) <- orfs_1000_TEW$idv_name

# write aa and nt orfs to file
writeXStringSet(orfs_aa_seq_TEA, paste0("t_elegans_DNAtransp.auto.fasta", "_aa_orfs.auto.fa"))
writeXStringSet(orfs_seq_TEA, paste0("t_elegans_DNAtransp.auto.fasta", "_nt_orfs.auto.fa"))
file.rename(from = "t_elegans_DNAtransp.auto.fasta_aa_orfs.auto.fa", to = "t_elegans_DNAtransp_aa_orfs.auto.fa")
file.rename(from = "t_elegans_DNAtransp.auto.fasta_nt_orfs.auto.fa", to = "t_elegans_DNAtransp_nt_orfs.auto.fa")

writeXStringSet(orfs_aa_seq_TEZ, paste0("t_elegans_DNAtransp.chrZ.fasta", "_aa_orfs.chrZ.fa"))
writeXStringSet(orfs_seq_TEZ, paste0("t_elegans_DNAtransp.chrZ.fasta", "_nt_orfs.chrZ.fa"))
file.rename(from = "t_elegans_DNAtransp.chrZ.fasta_aa_orfs.chrZ.fa", to = "t_elegans_DNAtransp_aa_orfs.chrZ.fa")
file.rename(from = "t_elegans_DNAtransp.chrZ.fasta_nt_orfs.chrZ.fa", to = "t_elegans_DNAtransp_nt_orfs.chrZ.fa")

writeXStringSet(orfs_aa_seq_TEW, paste0("t_elegans_DNAtransp.chrW.fasta", "_aa_orfs.chrW.fa"))
writeXStringSet(orfs_seq_TEW, paste0("t_elegans_DNAtransp.chrW.fasta", "_nt_orfs.chrW.fa"))
file.rename(from = "t_elegans_DNAtransp.chrW.fasta_aa_orfs.chrW.fa", to = "t_elegans_DNAtransp_aa_orfs.chrW.fa")
file.rename(from = "t_elegans_DNAtransp.chrW.fasta_nt_orfs.chrW.fa", to = "t_elegans_DNAtransp_nt_orfs.chrW.fa")

--------------------------------------------------------------------------------
  # scp '_aa_orf_' files back into Xenomorph and run RPSBLAST 
  # Read output back into R and proceed with next steps
  --------------------------------------------------------------------------------
  
  # read in RPSBLAST output
orfs_aa_rps_TEA <- read_tsv("t_elegans_DNAtransp_rpsblast.auto.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))
orfs_aa_rps_TEZ <- read_tsv("t_elegans_DNAtransp_rpsblast.chrZ.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))
orfs_aa_rps_TEW <- read_tsv("t_elegans_DNAtransp_rpsblast.chrW.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))

# split rpsblast title to make usable
orfs_aa_rps_TEA <- orfs_aa_rps_TEA %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")
orfs_aa_rps_TEZ <- orfs_aa_rps_TEZ %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")
orfs_aa_rps_TEW <- orfs_aa_rps_TEW %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")

# determine presence/absence of RT and EN
orfs_aa_rps_rt_TEA <- orfs_aa_rps_TEA %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_TEA <- orfs_aa_rps_TEA %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

orfs_aa_rps_rt_TEZ <- orfs_aa_rps_TEZ %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_TEZ <- orfs_aa_rps_TEZ %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

orfs_aa_rps_rt_TEW <- orfs_aa_rps_TEW %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_TEW <- orfs_aa_rps_TEW %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

# give seqnames their original names without the orf ids
orfs_aa_rps_rt_TEA = orfs_aa_rps_rt_TEA %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_TEA = orfs_aa_rps_en_TEA %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

orfs_aa_rps_rt_TEZ = orfs_aa_rps_rt_TEZ %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_TEZ = orfs_aa_rps_en_TEZ %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

orfs_aa_rps_rt_TEW = orfs_aa_rps_rt_TEW %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_TEW = orfs_aa_rps_en_TEW %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

# select orfs with intact rt and en
intact_orfs_TEA <- orfs_aa_rps_rt_TEA %>% filter(seqnames %in% orfs_aa_rps_en_TEA$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_TEA = intact_orfs_TEA %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_TEA = intact_orfs_bed_TEA[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_TEA), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("t_elegans_auto", "_intact_DNAorfs.bed"))

intact_orfs_TEZ <- orfs_aa_rps_rt_TEZ %>% filter(seqnames %in% orfs_aa_rps_en_TEZ$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_TEZ = intact_orfs_TEZ %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_TEZ = intact_orfs_bed_TEZ[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_TEZ), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("t_elegans_chrZ", "_intact_DNAorfs.bed"))

intact_orfs_TEW <- orfs_aa_rps_rt_TEW %>% filter(seqnames %in% orfs_aa_rps_en_TEW$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_TEW = intact_orfs_TEW %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_TEW = intact_orfs_bed_TEW[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_TEW), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("t_elegans_chrW", "_intact_DNAorfs.bed"))
------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  rm(list = ls())
### Naja naja ###
# read in DNAtransp sequences
raw_seq_naja_auto <- readDNAStringSet("n_naja_DNAtransp.auto.fasta")
raw_seq_naja_chrZ <- readDNAStringSet("n_naja_DNAtransp.chrZ.fasta")
raw_seq_naja_chrW <- readDNAStringSet("n_naja_DNAtransp.chrW.fasta")

# create pseudoranges tibble of CR1s
raw_tbl_naja_auto <- tibble(seqnames = names(raw_seq_naja_auto), start = 1, end = width(raw_seq_naja_auto)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")
raw_tbl_naja_chrZ <- tibble(seqnames = names(raw_seq_naja_chrZ), start = 1, end = width(raw_seq_naja_chrZ)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")
raw_tbl_naja_chrW <- tibble(seqnames = names(raw_seq_naja_chrW), start = 1, end = width(raw_seq_naja_chrW)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")

# create table of CR1 seqs with names
seq_group_tbl_naja_auto <- raw_tbl_naja_auto %>% dplyr::select(seqnames)
seq_group_tbl_naja_chrZ <- raw_tbl_naja_chrZ %>% dplyr::select(seqnames)
seq_group_tbl_naja_chrW <- raw_tbl_naja_chrW %>% dplyr::select(seqnames)

# find orfs over 1000bp in sequences
#orfs_1000 <- ORFik::findORFs(raw_seq, startCodon = startDefinition(1), minimumLength = 1000) %>% as_tibble()

orfs_1000_NNA = findORFsFasta("n_naja_DNAtransp.auto.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_NNA = as_tibble(orfs_1000_NNA)
orfs_1000_NNA$idv_name = paste0(orfs_1000_NNA$seqnames, "#orf", 1:nrow(orfs_1000_NNA))

orfs_1000_NNZ = findORFsFasta("n_naja_DNAtransp.chrZ.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_NNZ = as_tibble(orfs_1000_NNZ)
orfs_1000_NNZ$idv_name = paste0(orfs_1000_NNZ$seqnames, "#orf", 1:nrow(orfs_1000_NNZ))

orfs_1000_NNW = findORFsFasta("n_naja_DNAtransp.chrW.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_NNW = as_tibble(orfs_1000_NNW)
orfs_1000_NNW$idv_name = paste0(orfs_1000_NNW$seqnames, "#orf", 1:nrow(orfs_1000_NNW))

# make ranges object of orfs
orfs_1000_NNA_ranges <- GRanges(seqnames = orfs_1000_NNA$seqnames, ranges = IRanges(start = orfs_1000_NNA$start, end = orfs_1000_NNA$end))
orfs_1000_NNZ_ranges <- GRanges(seqnames = orfs_1000_NNZ$seqnames, ranges = IRanges(start = orfs_1000_NNZ$start, end = orfs_1000_NNZ$end))
orfs_1000_NNW_ranges <- GRanges(seqnames = orfs_1000_NNW$seqnames, ranges = IRanges(start = orfs_1000_NNW$start, end = orfs_1000_NNW$end))

# get seq of orfs and name orfs
orfs_seq_NNA <- Biostrings::getSeq(raw_seq_naja_auto, orfs_1000_NNA_ranges)
names(orfs_seq_NNA) <- orfs_1000_NNA$idv_name

orfs_seq_NNZ <- Biostrings::getSeq(raw_seq_naja_chrZ, orfs_1000_NNZ_ranges)
names(orfs_seq_NNZ) <- orfs_1000_NNZ$idv_name

orfs_seq_NNW <- Biostrings::getSeq(raw_seq_naja_chrW, orfs_1000_NNW_ranges)
names(orfs_seq_NNW) <- orfs_1000_NNW$idv_name

# translate orf
orfs_aa_seq_NNA <- translate(orfs_seq_NNA, if.fuzzy.codon = "solve")
names(orfs_aa_seq_NNA) <- orfs_1000_NNA$idv_name

orfs_aa_seq_NNZ <- translate(orfs_seq_NNZ, if.fuzzy.codon = "solve")
names(orfs_aa_seq_NNZ) <- orfs_1000_NNZ$idv_name

orfs_aa_seq_NNW <- translate(orfs_seq_NNW, if.fuzzy.codon = "solve")
names(orfs_aa_seq_NNW) <- orfs_1000_NNW$idv_name

# write nt and nt orfs to file
writeXStringSet(orfs_aa_seq_NNA, paste0("n_naja_DNAtransp.auto.fasta", "_aa_orfs.auto.fa"))
writeXStringSet(orfs_seq_NNA, paste0("n_naja_DNAtransp.auto.fasta", "_nt_orfs.auto.fa"))
file.rename(from = "n_naja_DNAtransp.auto.fasta_aa_orfs.auto.fa", to = "n_naja_DNAtransp_aa_orfs.auto.fa")
file.rename(from = "n_naja_DNAtransp.auto.fasta_nt_orfs.auto.fa", to = "n_naja_DNAtransp_nt_orfs.auto.fa")

writeXStringSet(orfs_aa_seq_NNZ, paste0("n_naja_DNAtransp.chrZ.fasta", "_aa_orfs.chrZ.fa"))
writeXStringSet(orfs_seq_NNZ, paste0("n_naja_DNAtransp.chrZ.fasta", "_nt_orfs.chrZ.fa"))
file.rename(from = "n_naja_DNAtransp.chrZ.fasta_aa_orfs.chrZ.fa", to = "n_naja_DNAtransp_aa_orfs.chrZ.fa")
file.rename(from = "n_naja_DNAtransp.chrZ.fasta_nt_orfs.chrZ.fa", to = "n_naja_DNAtransp_nt_orfs.chrZ.fa")

writeXStringSet(orfs_aa_seq_NNW, paste0("n_naja_DNAtransp.chrW.fasta", "_aa_orfs.chrW.fa"))
writeXStringSet(orfs_seq_NNW, paste0("n_naja_DNAtransp.chrW.fasta", "_nt_orfs.chrW.fa"))
file.rename(from = "n_naja_DNAtransp.chrW.fasta_aa_orfs.chrW.fa", to = "n_naja_DNAtransp_aa_orfs.chrW.fa")
file.rename(from = "n_naja_DNAtransp.chrW.fasta_nt_orfs.chrW.fa", to = "n_naja_DNAtransp_nt_orfs.chrW.fa")

--------------------------------------------------------------------------------
  # scp '_aa_orf_' files back into Xenomorph and run RPSBLAST 
  # Read output back into R and proceed with next steps
  --------------------------------------------------------------------------------
  
  # read in RPSBLAST output
orfs_aa_rps_NNA <- read_tsv("n_naja_DNAtransp_rpsblast.auto.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))
orfs_aa_rps_NNZ <- read_tsv("n_naja_DNAtransp_rpsblast.chrZ.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))
orfs_aa_rps_NNW <- read_tsv("n_naja_DNAtransp_rpsblast.chrW.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))

# split rpsblast title to make usable
orfs_aa_rps_NNA <- orfs_aa_rps_NNA %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")
orfs_aa_rps_NNZ <- orfs_aa_rps_NNZ %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")
orfs_aa_rps_NNW <- orfs_aa_rps_NNW %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")

# determine presence/absence of RT and EN
orfs_aa_rps_rt_NNA <- orfs_aa_rps_NNA %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_NNA <- orfs_aa_rps_NNA %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

orfs_aa_rps_rt_NNZ <- orfs_aa_rps_NNZ %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_NNZ <- orfs_aa_rps_NNZ %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

orfs_aa_rps_rt_NNW <- orfs_aa_rps_NNW %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_NNW <- orfs_aa_rps_NNW %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)


# give seqnames their original names without the orf ids
orfs_aa_rps_rt_NNA = orfs_aa_rps_rt_NNA %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_NNA = orfs_aa_rps_en_NNA %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

orfs_aa_rps_rt_NNZ = orfs_aa_rps_rt_NNZ %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_NNZ = orfs_aa_rps_en_NNZ %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

orfs_aa_rps_rt_NNW = orfs_aa_rps_rt_NNW %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_NNW = orfs_aa_rps_en_NNW %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

# select orfs with intact rt and en
intact_orfs_NNA <- orfs_aa_rps_rt_NNA %>% filter(seqnames %in% orfs_aa_rps_en_NNA$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_NNA = intact_orfs_NNA %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_NNA = intact_orfs_bed_NNA[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_NNA), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("n_naja_auto", "_intact_DNAorfs.bed"))

intact_orfs_NNZ <- orfs_aa_rps_rt_NNZ %>% filter(seqnames %in% orfs_aa_rps_en_NNZ$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_NNZ = intact_orfs_NNZ %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_NNZ = intact_orfs_bed_NNZ[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_NNZ), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("n_naja_chrZ", "_intact_DNAorfs.bed"))

intact_orfs_NNW <- orfs_aa_rps_rt_NNW %>% filter(seqnames %in% orfs_aa_rps_en_NNW$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_NNW = intact_orfs_NNW %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_NNW = intact_orfs_bed_NNW[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_NNW), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("n_naja_chrW", "_intact_DNAorfs.bed"))
------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  rm(list = ls())

### Cerastes gaperettii ###
# read in DNAtransp sequences
raw_seq_gasperettii_auto <- readDNAStringSet("c_gasperettii_DNAtransp.auto.fasta")
raw_seq_gasperettii_chrZ <- readDNAStringSet("c_gasperettii_DNAtransp.chrZ.fasta")
raw_seq_gasperettii_chrW <- readDNAStringSet("c_gasperettii_DNAtransp.chrW.fasta")

# create pseudoranges tibble of DNAtransp
raw_tbl_gasperettii_auto <- tibble(seqnames = names(raw_seq_gasperettii_auto), start = 1, end = width(raw_seq_gasperettii_auto)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")
raw_tbl_gasperettii_chrZ <- tibble(seqnames = names(raw_seq_gasperettii_chrZ), start = 1, end = width(raw_seq_gasperettii_chrZ)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")
raw_tbl_gasperettii_chrW <- tibble(seqnames = names(raw_seq_gasperettii_chrW), start = 1, end = width(raw_seq_gasperettii_chrW)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")

# create table of DNAtransp seqs with names
seq_group_tbl_gasperettii_auto <- raw_tbl_gasperettii_auto %>% dplyr::select(seqnames)
seq_group_tbl_gasperettii_chrZ <- raw_tbl_gasperettii_chrZ %>% dplyr::select(seqnames)
seq_group_tbl_gasperettii_chrW <- raw_tbl_gasperettii_chrW %>% dplyr::select(seqnames)

# find orfs over 1000bp in sequences
#orfs_1000 <- ORFik::findORFs(raw_seq, startCodon = startDefinition(1), minimumLength = 1000) %>% as_tibble()

orfs_1000_CGA = findORFsFasta("c_gasperettii_DNAtransp.auto.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_CGA = as_tibble(orfs_1000_CGA)
orfs_1000_CGA$idv_name = paste0(orfs_1000_CGA$seqnames, "#orf", 1:nrow(orfs_1000_CGA))

orfs_1000_CGZ = findORFsFasta("c_gasperettii_DNAtransp.chrZ.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_CGZ = as_tibble(orfs_1000_CGZ)
orfs_1000_CGZ$idv_name = paste0(orfs_1000_CGZ$seqnames, "#orf", 1:nrow(orfs_1000_CGZ))

orfs_1000_CGW = findORFsFasta("c_gasperettii_DNAtransp.chrW.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_CGW = as_tibble(orfs_1000_CGW)
orfs_1000_CGW$idv_name = paste0(orfs_1000_CGW$seqnames, "#orf", 1:nrow(orfs_1000_CGW))

# make ranges object of orfs
orfs_1000_CGA_ranges <- GRanges(seqnames = orfs_1000_CGA$seqnames, ranges = IRanges(start = orfs_1000_CGA$start, end = orfs_1000_CGA$end))
orfs_1000_CGZ_ranges <- GRanges(seqnames = orfs_1000_CGZ$seqnames, ranges = IRanges(start = orfs_1000_CGZ$start, end = orfs_1000_CGZ$end))
orfs_1000_CGW_ranges <- GRanges(seqnames = orfs_1000_CGW$seqnames, ranges = IRanges(start = orfs_1000_CGW$start, end = orfs_1000_CGW$end))

# get seq of orfs and name orfs
orfs_seq_CGA <- Biostrings::getSeq(raw_seq_gasperettii_auto, orfs_1000_CGA_ranges)
names(orfs_seq_CGA) <- orfs_1000_CGA$idv_name

orfs_seq_CGZ <- Biostrings::getSeq(raw_seq_gasperettii_chrZ, orfs_1000_CGZ_ranges)
names(orfs_seq_CGZ) <- orfs_1000_CGZ$idv_name

orfs_seq_CGW <- Biostrings::getSeq(raw_seq_gasperettii_chrW, orfs_1000_CGW_ranges)
names(orfs_seq_CGW) <- orfs_1000_CGW$idv_name

# translate orf
orfs_aa_seq_CGA <- translate(orfs_seq_CGA, if.fuzzy.codon = "solve")
names(orfs_aa_seq_CGA) <- orfs_1000_CGA$idv_name

orfs_aa_seq_CGZ <- translate(orfs_seq_CGZ, if.fuzzy.codon = "solve")
names(orfs_aa_seq_CGZ) <- orfs_1000_CGZ$idv_name

orfs_aa_seq_CGW <- translate(orfs_seq_CGW, if.fuzzy.codon = "solve")
names(orfs_aa_seq_CGW) <- orfs_1000_CGW$idv_name

# write nt and nt orfs to file
writeXStringSet(orfs_aa_seq_CGA, paste0("c_gasperettii_DNAtransp.auto.fasta", "_aa_orfs.auto.fa"))
writeXStringSet(orfs_seq_CGA, paste0("c_gasperettii_DNAtransp.auto.fasta", "_nt_orfs.auto.fa"))
file.rename(from = "c_gasperettii_DNAtransp.auto.fasta_aa_orfs.auto.fa", to = "c_gasperettii_DNAtransp_aa_orfs.auto.fa")
file.rename(from = "c_gasperettii_DNAtransp.auto.fasta_nt_orfs.auto.fa", to = "c_gasperettii_DNAtransp_nt_orfs.auto.fa")

writeXStringSet(orfs_aa_seq_CGZ, paste0("c_gasperettii_DNAtransp.chrZ.fasta", "_aa_orfs.chrZ.fa"))
writeXStringSet(orfs_seq_CGZ, paste0("c_gasperettii_DNAtransp.chrZ.fasta", "_nt_orfs.chrZ.fa"))
file.rename(from = "c_gasperettii_DNAtransp.chrZ.fasta_aa_orfs.chrZ.fa", to = "c_gasperettii_DNAtransp_aa_orfs.chrZ.fa")
file.rename(from = "c_gasperettii_DNAtransp.chrZ.fasta_nt_orfs.chrZ.fa", to = "c_gasperettii_DNAtransp_nt_orfs.chrZ.fa")

writeXStringSet(orfs_aa_seq_CGW, paste0("c_gasperettii_DNAtransp.chrW.fasta", "_aa_orfs.chrW.fa"))
writeXStringSet(orfs_seq_CGW, paste0("c_gasperettii_DNAtransp.chrW.fasta", "_nt_orfs.chrW.fa"))
file.rename(from = "c_gasperettii_DNAtransp.chrW.fasta_aa_orfs.chrW.fa", to = "c_gasperettii_DNAtransp_aa_orfs.chrW.fa")
file.rename(from = "c_gasperettii_DNAtransp.chrW.fasta_nt_orfs.chrW.fa", to = "c_gasperettii_DNAtransp_nt_orfs.chrW.fa")

--------------------------------------------------------------------------------
  # scp '_aa_orf_' files back into Xenomorph and run RPSBLAST 
  # Read output back into R and proceed with next steps
  --------------------------------------------------------------------------------
  
  # read in RPSBLAST output
orfs_aa_rps_CGA <- read_tsv("c_gasperettii_DNAtransp_rpsblast.auto.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))
orfs_aa_rps_CGZ <- read_tsv("c_gasperettii_DNAtransp_rpsblast.chrZ.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))
orfs_aa_rps_CGW <- read_tsv("c_gasperettii_DNAtransp_rpsblast.chrW.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))

# split rpsblast title to make usable
orfs_aa_rps_CGA <- orfs_aa_rps_CGA %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")
orfs_aa_rps_CGZ <- orfs_aa_rps_CGZ %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")
orfs_aa_rps_CGW <- orfs_aa_rps_CGW %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")

# determine presence/absence of RT and EN
orfs_aa_rps_rt_CGA <- orfs_aa_rps_CGA %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_CGA <- orfs_aa_rps_CGA %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

orfs_aa_rps_rt_CGZ <- orfs_aa_rps_CGZ %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_CGZ <- orfs_aa_rps_CGZ %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

orfs_aa_rps_rt_CGW <- orfs_aa_rps_CGW %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_CGW <- orfs_aa_rps_CGW %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

# give seqnames their original names without the orf ids
orfs_aa_rps_rt_CGA = orfs_aa_rps_rt_CGA %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_CGA = orfs_aa_rps_en_CGA %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

orfs_aa_rps_rt_CGZ = orfs_aa_rps_rt_CGZ %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_CGZ = orfs_aa_rps_en_CGZ %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

orfs_aa_rps_rt_CGW = orfs_aa_rps_rt_CGW %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_CGW = orfs_aa_rps_en_CGW %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

# select orfs with intact rt and en
intact_orfs_CGA <- orfs_aa_rps_rt_CGA %>% filter(seqnames %in% orfs_aa_rps_en_CGA$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_CGA = intact_orfs_CGA %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_CGA = intact_orfs_bed_CGA[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_CGA), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("c_gasperettii_auto", "_intact_DNAorfs.bed"))

intact_orfs_CGZ <- orfs_aa_rps_rt_CGZ %>% filter(seqnames %in% orfs_aa_rps_en_CGZ$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_CGZ = intact_orfs_CGZ %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_CGZ = intact_orfs_bed_CGZ[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_CGZ), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("c_gasperettii_chrZ", "_intact_DNAorfs.bed"))

intact_orfs_CGW <- orfs_aa_rps_rt_CGW %>% filter(seqnames %in% orfs_aa_rps_en_CGW$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_CGW = intact_orfs_CGW %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_CGW = intact_orfs_bed_CGW[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_CGW), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("c_gasperettii_chrW", "_intact_DNAorfs.bed"))

rm(list = ls())
-----------------------------------------------------------------------------------
  ### Vipera ursinii ###
  # read in DNAtransp sequences
raw_seq_ursinii_auto <- readDNAStringSet("v_ursinii_DNAtransp.auto.fasta")
raw_seq_ursinii_chrZ <- readDNAStringSet("v_ursinii_DNAtransp.chrZ.fasta")
raw_seq_ursinii_chrW <- readDNAStringSet("v_ursinii_DNAtransp.chrW.fasta")

# create pseudoranges tibble of DNAtransp
raw_tbl_ursinii_auto <- tibble(seqnames = names(raw_seq_ursinii_auto), start = 1, end = width(raw_seq_ursinii_auto)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")
raw_tbl_ursinii_chrZ <- tibble(seqnames = names(raw_seq_ursinii_chrZ), start = 1, end = width(raw_seq_ursinii_chrZ)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")
raw_tbl_ursinii_chrW <- tibble(seqnames = names(raw_seq_ursinii_chrW), start = 1, end = width(raw_seq_ursinii_chrW)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")

# create table of DNAtransp seqs with names
seq_group_tbl_ursinii_auto <- raw_tbl_ursinii_auto %>% dplyr::select(seqnames)
seq_group_tbl_ursinii_chrZ <- raw_tbl_ursinii_chrZ %>% dplyr::select(seqnames)
seq_group_tbl_ursinii_chrW <- raw_tbl_ursinii_chrW %>% dplyr::select(seqnames)

# find orfs over 1000bp in sequences
#orfs_1000 <- ORFik::findORFs(raw_seq, startCodon = startDefinition(1), minimumLength = 1000) %>% as_tibble()

orfs_1000_VUA = findORFsFasta("v_ursinii_DNAtransp.auto.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_VUA = as_tibble(orfs_1000_VUA)
orfs_1000_VUA$idv_name = paste0(orfs_1000_VUA$seqnames, "#orf", 1:nrow(orfs_1000_VUA))

orfs_1000_VUZ = findORFsFasta("v_ursinii_DNAtransp.chrZ.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_VUZ = as_tibble(orfs_1000_VUZ)
orfs_1000_VUZ$idv_name = paste0(orfs_1000_VUZ$seqnames, "#orf", 1:nrow(orfs_1000_VUZ))

orfs_1000_VUW = findORFsFasta("v_ursinii_DNAtransp.chrW.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_VUW = as_tibble(orfs_1000_VUW)
orfs_1000_VUW$idv_name = paste0(orfs_1000_VUW$seqnames, "#orf", 1:nrow(orfs_1000_VUW))

# make ranges object of orfs
orfs_1000_VUA_ranges <- GRanges(seqnames = orfs_1000_VUA$seqnames, ranges = IRanges(start = orfs_1000_VUA$start, end = orfs_1000_VUA$end))
orfs_1000_VUZ_ranges <- GRanges(seqnames = orfs_1000_VUZ$seqnames, ranges = IRanges(start = orfs_1000_VUZ$start, end = orfs_1000_VUZ$end))
orfs_1000_VUW_ranges <- GRanges(seqnames = orfs_1000_VUW$seqnames, ranges = IRanges(start = orfs_1000_VUW$start, end = orfs_1000_VUW$end))

# get seq of orfs and name orfs
orfs_seq_VUA <- Biostrings::getSeq(raw_seq_ursinii_auto, orfs_1000_VUA_ranges)
names(orfs_seq_VUA) <- orfs_1000_VUA$idv_name

orfs_seq_VUZ <- Biostrings::getSeq(raw_seq_ursinii_chrZ, orfs_1000_VUZ_ranges)
names(orfs_seq_VUZ) <- orfs_1000_VUZ$idv_name

orfs_seq_VUW <- Biostrings::getSeq(raw_seq_ursinii_chrW, orfs_1000_VUW_ranges)
names(orfs_seq_VUW) <- orfs_1000_VUW$idv_name

# translate orf
orfs_aa_seq_VUA <- translate(orfs_seq_VUA, if.fuzzy.codon = "solve")
names(orfs_aa_seq_VUA) <- orfs_1000_VUA$idv_name

orfs_aa_seq_VUZ <- translate(orfs_seq_VUZ, if.fuzzy.codon = "solve")
names(orfs_aa_seq_VUZ) <- orfs_1000_VUZ$idv_name

orfs_aa_seq_VUW <- translate(orfs_seq_VUW, if.fuzzy.codon = "solve")
names(orfs_aa_seq_VUW) <- orfs_1000_VUW$idv_name

# write nt and nt orfs to file
writeXStringSet(orfs_aa_seq_VUA, paste0("v_ursinii_DNAtransp.auto.fasta", "_aa_orfs.auto.fa"))
writeXStringSet(orfs_seq_VUA, paste0("v_ursinii_DNAtransp.auto.fasta", "_nt_orfs.auto.fa"))
file.rename(from = "v_ursinii_DNAtransp.auto.fasta_aa_orfs.auto.fa", to = "v_ursinii_DNAtransp_aa_orfs.auto.fa")
file.rename(from = "v_ursinii_DNAtransp.auto.fasta_nt_orfs.auto.fa", to = "v_ursinii_DNAtransp_nt_orfs.auto.fa")

writeXStringSet(orfs_aa_seq_VUZ, paste0("v_ursinii_DNAtransp.chrZ.fasta", "_aa_orfs.chrZ.fa"))
writeXStringSet(orfs_seq_VUZ, paste0("v_ursinii_DNAtransp.chrZ.fasta", "_nt_orfs.chrZ.fa"))
file.rename(from = "v_ursinii_DNAtransp.chrZ.fasta_aa_orfs.chrZ.fa", to = "v_ursinii_DNAtransp_aa_orfs.chrZ.fa")
file.rename(from = "v_ursinii_DNAtransp.chrZ.fasta_nt_orfs.chrZ.fa", to = "v_ursinii_DNAtransp_nt_orfs.chrZ.fa")

writeXStringSet(orfs_aa_seq_VUW, paste0("v_ursinii_DNAtransp.chrW.fasta", "_aa_orfs.chrW.fa"))
writeXStringSet(orfs_seq_VUW, paste0("v_ursinii_DNAtransp.chrW.fasta", "_nt_orfs.chrW.fa"))
file.rename(from = "v_ursinii_DNAtransp.chrW.fasta_aa_orfs.chrW.fa", to = "v_ursinii_DNAtransp_aa_orfs.chrW.fa")
file.rename(from = "v_ursinii_DNAtransp.chrW.fasta_nt_orfs.chrW.fa", to = "v_ursinii_DNAtransp_nt_orfs.chrW.fa")

--------------------------------------------------------------------------------
  # scp '_aa_orf_' files back into Xenomorph and run RPSBLAST 
  # Read output back into R and proceed with next steps
  --------------------------------------------------------------------------------
  
  # read in RPSBLAST output
orfs_aa_rps_VUA <- read_tsv("v_ursinii_DNAtransp_rpsblast.auto.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))
orfs_aa_rps_VUZ <- read_tsv("v_ursinii_DNAtransp_rpsblast.chrZ.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))
orfs_aa_rps_VUW <- read_tsv("v_ursinii_DNAtransp_rpsblast.chrW.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))

# split rpsblast title to make usable
orfs_aa_rps_VUA <- orfs_aa_rps_VUA %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")
orfs_aa_rps_VUZ <- orfs_aa_rps_VUZ %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")
orfs_aa_rps_VUW <- orfs_aa_rps_VUW %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")

# determine presence/absence of RT and EN
orfs_aa_rps_rt_VUA <- orfs_aa_rps_VUA %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_VUA <- orfs_aa_rps_VUA %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

orfs_aa_rps_rt_VUZ <- orfs_aa_rps_VUZ %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_VUZ <- orfs_aa_rps_VUZ %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

orfs_aa_rps_rt_VUW <- orfs_aa_rps_VUW %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_VUW <- orfs_aa_rps_VUW %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

# give seqnames their original names without the orf ids
orfs_aa_rps_rt_VUA = orfs_aa_rps_rt_VUA %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_VUA = orfs_aa_rps_en_VUA %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

orfs_aa_rps_rt_VUZ = orfs_aa_rps_rt_VUZ %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_VUZ = orfs_aa_rps_en_VUZ %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

orfs_aa_rps_rt_VUW = orfs_aa_rps_rt_VUW %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_VUW = orfs_aa_rps_en_VUW %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

# select orfs with intact rt and en
intact_orfs_VUA <- orfs_aa_rps_rt_VUA %>% filter(seqnames %in% orfs_aa_rps_en_VUA$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_VUA = intact_orfs_VUA %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_VUA = intact_orfs_bed_VUA[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_VUA), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("v_ursinii_auto", "_intact_DNAorfs.bed"))

intact_orfs_VUZ <- orfs_aa_rps_rt_VUZ %>% filter(seqnames %in% orfs_aa_rps_en_VUZ$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_VUZ = intact_orfs_VUZ %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_VUZ = intact_orfs_bed_VUZ[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_VUZ), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("v_ursinii_chrZ", "_intact_DNAorfs.bed"))

intact_orfs_VUW <- orfs_aa_rps_rt_VUW %>% filter(seqnames %in% orfs_aa_rps_en_VUW$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_VUW = intact_orfs_VUW %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_VUW = intact_orfs_bed_VUW[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_VUW), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("v_ursinii_chrW", "_intact_DNArfs.bed"))

rm(list = ls())
-----------------------------------------------------------------------------------
  ### Vipera berus ###
  # read in DNAtransp sequences
raw_seq_berus_auto <- readDNAStringSet("v_berus_DNAtransp.auto.fasta")
raw_seq_berus_chrZ <- readDNAStringSet("v_berus_DNAtransp.chrZ.fasta")
raw_seq_berus_chrW <- readDNAStringSet("v_berus_DNAtransp.chrW.fasta")

# create pseudoranges tibble of DNAtransp
raw_tbl_berus_auto <- tibble(seqnames = names(raw_seq_berus_auto), start = 1, end = width(raw_seq_berus_auto)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")
raw_tbl_berus_chrZ <- tibble(seqnames = names(raw_seq_berus_chrZ), start = 1, end = width(raw_seq_berus_chrZ)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")
raw_tbl_berus_chrW <- tibble(seqnames = names(raw_seq_berus_chrW), start = 1, end = width(raw_seq_berus_chrW)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")

# create table of DNAtransp seqs with names
seq_group_tbl_berus_auto <- raw_tbl_berus_auto %>% dplyr::select(seqnames)
seq_group_tbl_berus_chrZ <- raw_tbl_berus_chrZ %>% dplyr::select(seqnames)
seq_group_tbl_berus_chrW <- raw_tbl_berus_chrW %>% dplyr::select(seqnames)

# find orfs over 1000bp in sequences
#orfs_1000 <- ORFik::findORFs(raw_seq, startCodon = startDefinition(1), minimumLength = 1000) %>% as_tibble()

orfs_1000_VBA = findORFsFasta("v_berus_DNAtransp.auto.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_VBA = as_tibble(orfs_1000_VBA)
orfs_1000_VBA$idv_name = paste0(orfs_1000_VBA$seqnames, "#orf", 1:nrow(orfs_1000_VBA))

orfs_1000_VBZ = findORFsFasta("v_berus_DNAtransp.chrZ.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_VBZ = as_tibble(orfs_1000_VBZ)
orfs_1000_VBZ$idv_name = paste0(orfs_1000_VBZ$seqnames, "#orf", 1:nrow(orfs_1000_VBZ))

orfs_1000_VBW = findORFsFasta("v_berus_DNAtransp.chrW.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_VBW = as_tibble(orfs_1000_VBW)
orfs_1000_VBW$idv_name = paste0(orfs_1000_VBW$seqnames, "#orf", 1:nrow(orfs_1000_VBW))

# make ranges object of orfs
orfs_1000_VBA_ranges <- GRanges(seqnames = orfs_1000_VBA$seqnames, ranges = IRanges(start = orfs_1000_VBA$start, end = orfs_1000_VBA$end))
orfs_1000_VBZ_ranges <- GRanges(seqnames = orfs_1000_VBZ$seqnames, ranges = IRanges(start = orfs_1000_VBZ$start, end = orfs_1000_VBZ$end))
orfs_1000_VBW_ranges <- GRanges(seqnames = orfs_1000_VBW$seqnames, ranges = IRanges(start = orfs_1000_VBW$start, end = orfs_1000_VBW$end))

# get seq of orfs and name orfs
orfs_seq_VBA <- Biostrings::getSeq(raw_seq_berus_auto, orfs_1000_VBA_ranges)
names(orfs_seq_VBA) <- orfs_1000_VBA$idv_name

orfs_seq_VBZ <- Biostrings::getSeq(raw_seq_berus_chrZ, orfs_1000_VBZ_ranges)
names(orfs_seq_VBZ) <- orfs_1000_VBZ$idv_name

orfs_seq_VBW <- Biostrings::getSeq(raw_seq_berus_chrW, orfs_1000_VBW_ranges)
names(orfs_seq_VBW) <- orfs_1000_VBW$idv_name

# translate orf
orfs_aa_seq_VBA <- translate(orfs_seq_VBA, if.fuzzy.codon = "solve")
names(orfs_aa_seq_VBA) <- orfs_1000_VBA$idv_name

orfs_aa_seq_VBZ <- translate(orfs_seq_VBZ, if.fuzzy.codon = "solve")
names(orfs_aa_seq_VBZ) <- orfs_1000_VBZ$idv_name

orfs_aa_seq_VBW <- translate(orfs_seq_VBW, if.fuzzy.codon = "solve")
names(orfs_aa_seq_VBW) <- orfs_1000_VBW$idv_name

# write nt and nt orfs to file
writeXStringSet(orfs_aa_seq_VBA, paste0("v_berus_DNAtransp.auto.fasta", "_aa_orfs.auto.fa"))
writeXStringSet(orfs_seq_VBA, paste0("v_berus_DNAtransp.auto.fasta", "_nt_orfs.auto.fa"))
file.rename(from = "v_berus_DNAtransp.auto.fasta_aa_orfs.auto.fa", to = "v_berus_DNAtransp_aa_orfs.auto.fa")
file.rename(from = "v_berus_DNAtransp.auto.fasta_nt_orfs.auto.fa", to = "v_berus_DNAtransp_nt_orfs.auto.fa")

writeXStringSet(orfs_aa_seq_VBZ, paste0("v_berus_DNAtransp.chrZ.fasta", "_aa_orfs.chrZ.fa"))
writeXStringSet(orfs_seq_VBZ, paste0("v_berus_DNAtransp.chrZ.fasta", "_nt_orfs.chrZ.fa"))
file.rename(from = "v_berus_DNAtransp.chrZ.fasta_aa_orfs.chrZ.fa", to = "v_berus_DNAtransp_aa_orfs.chrZ.fa")
file.rename(from = "v_berus_DNAtransp.chrZ.fasta_nt_orfs.chrZ.fa", to = "v_berus_DNAtransp_nt_orfs.chrZ.fa")

writeXStringSet(orfs_aa_seq_VBW, paste0("v_berus_DNAtransp.chrW.fasta", "_aa_orfs.chrW.fa"))
writeXStringSet(orfs_seq_VBW, paste0("v_berus_DNAtransp.chrW.fasta", "_nt_orfs.chrW.fa"))
file.rename(from = "v_berus_DNAtransp.chrW.fasta_aa_orfs.chrW.fa", to = "v_berus_DNAtransp_aa_orfs.chrW.fa")
file.rename(from = "v_berus_DNAtransp.chrW.fasta_nt_orfs.chrW.fa", to = "v_berus_DNAtransp_nt_orfs.chrW.fa")

--------------------------------------------------------------------------------
  # scp '_aa_orf_' files back into Xenomorph and run RPSBLAST 
  # Read output back into R and proceed with next steps
  --------------------------------------------------------------------------------
  
  # read in RPSBLAST output
orfs_aa_rps_VBA <- read_tsv("v_berus_DNAtransp_rpsblast.auto.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))
orfs_aa_rps_VBZ <- read_tsv("v_berus_DNAtransp_rpsblast.chrZ.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))
orfs_aa_rps_VBW <- read_tsv("v_berus_DNAtransp_rpsblast.chrW.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))

# split rpsblast title to make usable
orfs_aa_rps_VBA <- orfs_aa_rps_VBA %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")
orfs_aa_rps_VBZ <- orfs_aa_rps_VBZ %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")
orfs_aa_rps_VBW <- orfs_aa_rps_VBW %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")

# determine presence/absence of RT and EN
orfs_aa_rps_rt_VBA <- orfs_aa_rps_VBA %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_VBA <- orfs_aa_rps_VBA %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

orfs_aa_rps_rt_VBZ <- orfs_aa_rps_VBZ %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_VBZ <- orfs_aa_rps_VBZ %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

orfs_aa_rps_rt_VBW <- orfs_aa_rps_VBW %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_VBW <- orfs_aa_rps_VBW %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

# give seqnames their original names without the orf ids
orfs_aa_rps_rt_VBA = orfs_aa_rps_rt_VBA %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_VBA = orfs_aa_rps_en_VBA %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

orfs_aa_rps_rt_VBZ = orfs_aa_rps_rt_VBZ %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_VBZ = orfs_aa_rps_en_VBZ %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

orfs_aa_rps_rt_VBW = orfs_aa_rps_rt_VBW %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_VBW = orfs_aa_rps_en_VBW %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

# select orfs with intact rt and en
intact_orfs_VBA <- orfs_aa_rps_rt_VBA %>% filter(seqnames %in% orfs_aa_rps_en_VBA$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_VBA = intact_orfs_VBA %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_VBA = intact_orfs_bed_VBA[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_VBA), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("v_berus_auto", "_intact_DNAorfs.bed"))

intact_orfs_VBZ <- orfs_aa_rps_rt_VBZ %>% filter(seqnames %in% orfs_aa_rps_en_VBZ$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_VBZ = intact_orfs_VBZ %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_VBZ = intact_orfs_bed_VBZ[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_VBZ), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("v_berus_chrZ", "_intact_DNAorfs.bed"))

intact_orfs_VBW <- orfs_aa_rps_rt_VBW %>% filter(seqnames %in% orfs_aa_rps_en_VBW$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_VBW = intact_orfs_VBW %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_VBW = intact_orfs_bed_VBW[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_VBW), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("v_berus_chrW", "_intact_DNAorfs.bed"))

rm(list = ls())
-----------------------------------------------------------------------------------
  ### Deinagkistrodon acutus ###
  # read in DNAtransp sequences
raw_seq_acutus_auto <- readDNAStringSet("d_acutus_DNAtransp.auto.fasta")
raw_seq_acutus_chrZ <- readDNAStringSet("d_acutus_DNAtransp.chrZ.fasta")
raw_seq_acutus_chrW <- readDNAStringSet("d_acutus_DNAtransp.chrW.fasta")

# create pseudoranges tibble of DNAtransp
raw_tbl_acutus_auto <- tibble(seqnames = names(raw_seq_acutus_auto), start = 1, end = width(raw_seq_acutus_auto)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")
raw_tbl_acutus_chrZ <- tibble(seqnames = names(raw_seq_acutus_chrZ), start = 1, end = width(raw_seq_acutus_chrZ)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")
raw_tbl_acutus_chrW <- tibble(seqnames = names(raw_seq_acutus_chrW), start = 1, end = width(raw_seq_acutus_chrW)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")

# create table of DNAtransp seqs with names
seq_group_tbl_acutus_auto <- raw_tbl_acutus_auto %>% dplyr::select(seqnames)
seq_group_tbl_acutus_chrZ <- raw_tbl_acutus_chrZ %>% dplyr::select(seqnames)
seq_group_tbl_acutus_chrW <- raw_tbl_acutus_chrW %>% dplyr::select(seqnames)

# find orfs over 1000bp in sequences
#orfs_1000 <- ORFik::findORFs(raw_seq, startCodon = startDefinition(1), minimumLength = 1000) %>% as_tibble()

orfs_1000_DAA = findORFsFasta("d_acutus_DNAtransp.auto.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_DAA = as_tibble(orfs_1000_DAA)
orfs_1000_DAA$idv_name = paste0(orfs_1000_DAA$seqnames, "#orf", 1:nrow(orfs_1000_DAA))
# There were no orfs over 1000bp in autosome sequences. When I changed to 200bp, there were many. At 500bp, there were 3.

orfs_1000_DAZ = findORFsFasta("d_acutus_DNAtransp.chrZ.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_DAZ = as_tibble(orfs_1000_DAZ)
orfs_1000_DAZ$idv_name = paste0(orfs_1000_DAZ$seqnames, "#orf", 1:nrow(orfs_1000_DAZ))
# no orfs over 1000bp in chrZ

orfs_1000_DAW = findORFsFasta("d_acutus_DNAtransp.chrW.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_DAW = as_tibble(orfs_1000_DAW)
orfs_1000_DAW$idv_name = paste0(orfs_1000_DAW$seqnames, "#orf", 1:nrow(orfs_1000_DAW))

# make ranges object of orfs
orfs_1000_DAA_ranges <- GRanges(seqnames = orfs_1000_DAA$seqnames, ranges = IRanges(start = orfs_1000_DAA$start, end = orfs_1000_DAA$end))
orfs_1000_DAZ_ranges <- GRanges(seqnames = orfs_1000_DAZ$seqnames, ranges = IRanges(start = orfs_1000_DAZ$start, end = orfs_1000_DAZ$end))
orfs_1000_DAW_ranges <- GRanges(seqnames = orfs_1000_DAW$seqnames, ranges = IRanges(start = orfs_1000_DAW$start, end = orfs_1000_DAW$end))

# get seq of orfs and name orfs
orfs_seq_DAA <- Biostrings::getSeq(raw_seq_acutus_auto, orfs_1000_DAA_ranges)
names(orfs_seq_DAA) <- orfs_1000_DAA$idv_name

orfs_seq_DAZ <- Biostrings::getSeq(raw_seq_acutus_chrZ, orfs_1000_DAZ_ranges)
names(orfs_seq_DAZ) <- orfs_1000_DAZ$idv_name

orfs_seq_DAW <- Biostrings::getSeq(raw_seq_acutus_chrW, orfs_1000_DAW_ranges)
names(orfs_seq_DAW) <- orfs_1000_DAW$idv_name

# translate orf
orfs_aa_seq_DAA <- translate(orfs_seq_DAA, if.fuzzy.codon = "solve")
names(orfs_aa_seq_DAA) <- orfs_1000_DAA$idv_name

orfs_aa_seq_DAZ <- translate(orfs_seq_DAZ, if.fuzzy.codon = "solve")
names(orfs_aa_seq_DAZ) <- orfs_1000_DAZ$idv_name

orfs_aa_seq_DAW <- translate(orfs_seq_DAW, if.fuzzy.codon = "solve")
names(orfs_aa_seq_DAW) <- orfs_1000_DAW$idv_name

# write nt and nt orfs to file
writeXStringSet(orfs_aa_seq_DAA, paste0("d_acutus_DNAtransp.auto.fasta", "_aa_orfs.auto.fa"))
writeXStringSet(orfs_seq_DAA, paste0("d_acutus_DNAtransp.auto.fasta", "_nt_orfs.auto.fa"))
file.rename(from = "d_acutus_DNAtransp.auto.fasta_aa_orfs.auto.fa", to = "d_acutus_DNAtransp_aa_orfs.auto.fa")
file.rename(from = "d_acutus_DNAtransp.auto.fasta_nt_orfs.auto.fa", to = "d_acutus_DNAtransp_nt_orfs.auto.fa")

writeXStringSet(orfs_aa_seq_DAZ, paste0("d_acutus_DNAtransp.chrZ.fasta", "_aa_orfs.chrZ.fa"))
writeXStringSet(orfs_seq_DAZ, paste0("d_acutus_DNAtransp.chrZ.fasta", "_nt_orfs.chrZ.fa"))
file.rename(from = "d_acutus_DNAtransp.chrZ.fasta_aa_orfs.chrZ.fa", to = "d_acutus_DNAtransp_aa_orfs.chrZ.fa")
file.rename(from = "d_acutus_DNAtransp.chrZ.fasta_nt_orfs.chrZ.fa", to = "d_acutus_DNAtransp_nt_orfs.chrZ.fa")

writeXStringSet(orfs_aa_seq_DAW, paste0("d_acutus_DNAtransp.chrW.fasta", "_aa_orfs.chrW.fa"))
writeXStringSet(orfs_seq_DAW, paste0("d_acutus_DNAtransp.chrW.fasta", "_nt_orfs.chrW.fa"))
file.rename(from = "d_acutus_DNAtransp.chrW.fasta_aa_orfs.chrW.fa", to = "d_acutus_DNAtransp_aa_orfs.chrW.fa")
file.rename(from = "d_acutus_DNAtransp.chrW.fasta_nt_orfs.chrW.fa", to = "d_acutus_DNAtransp_nt_orfs.chrW.fa")

--------------------------------------------------------------------------------
  # scp '_aa_orf_' files back into Xenomorph and run RPSBLAST 
  # Read output back into R and proceed with next steps
  --------------------------------------------------------------------------------
  
  # read in RPSBLAST output
orfs_aa_rps_DAA <- read_tsv("d_acutus_DNAtransp_rpsblast.auto.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))
orfs_aa_rps_DAZ <- read_tsv("d_acutus_DNAtransp_rpsblast.chrZ.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))
orfs_aa_rps_DAW <- read_tsv("d_acutus_DNAtransp_rpsblast.chrW.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))

# split rpsblast title to make usable
orfs_aa_rps_DAA <- orfs_aa_rps_DAA %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")
orfs_aa_rps_DAZ <- orfs_aa_rps_DAZ %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")
orfs_aa_rps_DAW <- orfs_aa_rps_DAW %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")

# determine presence/absence of RT and EN
orfs_aa_rps_rt_DAA <- orfs_aa_rps_DAA %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_DAA <- orfs_aa_rps_DAA %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

orfs_aa_rps_rt_DAZ <- orfs_aa_rps_DAZ %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_DAZ <- orfs_aa_rps_DAZ %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

orfs_aa_rps_rt_DAW <- orfs_aa_rps_DAW %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_DAW <- orfs_aa_rps_DAW %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

# give seqnames their original names without the orf ids
orfs_aa_rps_rt_DAA = orfs_aa_rps_rt_DAA %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_DAA = orfs_aa_rps_en_DAA %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

orfs_aa_rps_rt_DAZ = orfs_aa_rps_rt_DAZ %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_DAZ = orfs_aa_rps_en_DAZ %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

orfs_aa_rps_rt_DAW = orfs_aa_rps_rt_DAW %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_DAW = orfs_aa_rps_en_DAW %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

# select orfs with intact rt and en
intact_orfs_DAA <- orfs_aa_rps_rt_DAA %>% filter(seqnames %in% orfs_aa_rps_en_DAA$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_DAA = intact_orfs_DAA %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_DAA = intact_orfs_bed_DAA[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_DAA), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("d_acutus_auto", "_intact_DNAorfs.bed"))

intact_orfs_DAZ <- orfs_aa_rps_rt_DAZ %>% filter(seqnames %in% orfs_aa_rps_en_DAZ$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_DAZ = intact_orfs_DAZ %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_DAZ = intact_orfs_bed_DAZ[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_DAZ), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("d_acutus_chrZ", "_intact_DNAorfs.bed"))

intact_orfs_DAW <- orfs_aa_rps_rt_DAW %>% filter(seqnames %in% orfs_aa_rps_en_DAW$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_DAW = intact_orfs_DAW %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_DAW = intact_orfs_bed_DAW[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_DAW), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("d_acutus_chrW", "_intact_DNAorfs.bed"))


rm(list = ls())
-----------------------------------------------------------------------------------
  ### Crotalus adamanteus ###
  # read in DNAtransp sequences
raw_seq_adamanteus_auto <- readDNAStringSet("c_adamanteus_DNAtransp.auto.fasta")
raw_seq_adamanteus_chrZ <- readDNAStringSet("c_adamanteus_DNAtransp.chrZ.fasta")
raw_seq_adamanteus_chrW <- readDNAStringSet("c_adamanteus_DNAtransp.chrW.fasta")

# create pseudoranges tibble of DNAtransp
raw_tbl_adamanteus_auto <- tibble(seqnames = names(raw_seq_adamanteus_auto), start = 1, end = width(raw_seq_adamanteus_auto)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")
raw_tbl_adamanteus_chrZ <- tibble(seqnames = names(raw_seq_adamanteus_chrZ), start = 1, end = width(raw_seq_adamanteus_chrZ)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")
raw_tbl_adamanteus_chrW <- tibble(seqnames = names(raw_seq_adamanteus_chrW), start = 1, end = width(raw_seq_adamanteus_chrW)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")

# create table of DNAtransp seqs with names
seq_group_tbl_adamanteus_auto <- raw_tbl_adamanteus_auto %>% dplyr::select(seqnames)
seq_group_tbl_adamanteus_chrZ <- raw_tbl_adamanteus_chrZ %>% dplyr::select(seqnames)
seq_group_tbl_adamanteus_chrW <- raw_tbl_adamanteus_chrW %>% dplyr::select(seqnames)

# find orfs over 1000bp in sequences
#orfs_1000 <- ORFik::findORFs(raw_seq, startCodon = startDefinition(1), minimumLength = 1000) %>% as_tibble()

orfs_1000_CAA = findORFsFasta("c_adamanteus_DNAtransp.auto.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_CAA = as_tibble(orfs_1000_CAA)
orfs_1000_CAA$idv_name = paste0(orfs_1000_CAA$seqnames, "#orf", 1:nrow(orfs_1000_CAA))

orfs_1000_CAZ = findORFsFasta("c_adamanteus_DNAtransp.chrZ.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_CAZ = as_tibble(orfs_1000_CAZ)
orfs_1000_CAZ$idv_name = paste0(orfs_1000_CAZ$seqnames, "#orf", 1:nrow(orfs_1000_CAZ))

orfs_1000_CAW = findORFsFasta("c_adamanteus_DNAtransp.chrW.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_CAW = as_tibble(orfs_1000_CAW)
orfs_1000_CAW$idv_name = paste0(orfs_1000_CAW$seqnames, "#orf", 1:nrow(orfs_1000_CAW))

# make ranges object of orfs
orfs_1000_CAA_ranges <- GRanges(seqnames = orfs_1000_CAA$seqnames, ranges = IRanges(start = orfs_1000_CAA$start, end = orfs_1000_CAA$end))
orfs_1000_CAZ_ranges <- GRanges(seqnames = orfs_1000_CAZ$seqnames, ranges = IRanges(start = orfs_1000_CAZ$start, end = orfs_1000_CAZ$end))
orfs_1000_CAW_ranges <- GRanges(seqnames = orfs_1000_CAW$seqnames, ranges = IRanges(start = orfs_1000_CAW$start, end = orfs_1000_CAW$end))

# get seq of orfs and name orfs
orfs_seq_CAA <- Biostrings::getSeq(raw_seq_adamanteus_auto, orfs_1000_CAA_ranges)
names(orfs_seq_CAA) <- orfs_1000_CAA$idv_name

orfs_seq_CAZ <- Biostrings::getSeq(raw_seq_adamanteus_chrZ, orfs_1000_CAZ_ranges)
names(orfs_seq_CAZ) <- orfs_1000_CAZ$idv_name

orfs_seq_CAW <- Biostrings::getSeq(raw_seq_adamanteus_chrW, orfs_1000_CAW_ranges)
names(orfs_seq_CAW) <- orfs_1000_CAW$idv_name

# translate orf
orfs_aa_seq_CAA <- translate(orfs_seq_CAA, if.fuzzy.codon = "solve")
names(orfs_aa_seq_CAA) <- orfs_1000_CAA$idv_name

orfs_aa_seq_CAZ <- translate(orfs_seq_CAZ, if.fuzzy.codon = "solve")
names(orfs_aa_seq_CAZ) <- orfs_1000_CAZ$idv_name

orfs_aa_seq_CAW <- translate(orfs_seq_CAW, if.fuzzy.codon = "solve")
names(orfs_aa_seq_CAW) <- orfs_1000_CAW$idv_name

# write nt and nt orfs to file
writeXStringSet(orfs_aa_seq_CAA, paste0("c_adamanteus_DNAtransp.auto.fasta", "_aa_orfs.auto.fa"))
writeXStringSet(orfs_seq_CAA, paste0("c_adamanteus_DNAtransp.auto.fasta", "_nt_orfs.auto.fa"))
file.rename(from = "c_adamanteus_DNAtransp.auto.fasta_aa_orfs.auto.fa", to = "c_adamanteus_DNAtransp_aa_orfs.auto.fa")
file.rename(from = "c_adamanteus_DNAtransp.auto.fasta_nt_orfs.auto.fa", to = "c_adamanteus_DNAtransp_nt_orfs.auto.fa")

writeXStringSet(orfs_aa_seq_CAZ, paste0("c_adamanteus_DNAtransp.chrZ.fasta", "_aa_orfs.chrZ.fa"))
writeXStringSet(orfs_seq_CAZ, paste0("c_adamanteus_DNAtransp.chrZ.fasta", "_nt_orfs.chrZ.fa"))
file.rename(from = "c_adamanteus_DNAtransp.chrZ.fasta_aa_orfs.chrZ.fa", to = "c_adamanteus_DNAtransp_aa_orfs.chrZ.fa")
file.rename(from = "c_adamanteus_DNAtransp.chrZ.fasta_nt_orfs.chrZ.fa", to = "c_adamanteus_DNAtransp_nt_orfs.chrZ.fa")

writeXStringSet(orfs_aa_seq_CAW, paste0("c_adamanteus_DNAtransp.chrW.fasta", "_aa_orfs.chrW.fa"))
writeXStringSet(orfs_seq_CAW, paste0("c_adamanteus_DNAtransp.chrW.fasta", "_nt_orfs.chrW.fa"))
file.rename(from = "c_adamanteus_DNAtransp.chrW.fasta_aa_orfs.chrW.fa", to = "c_adamanteus_DNAtransp_aa_orfs.chrW.fa")
file.rename(from = "c_adamanteus_DNAtransp.chrW.fasta_nt_orfs.chrW.fa", to = "c_adamanteus_DNAtransp_nt_orfs.chrW.fa")

--------------------------------------------------------------------------------
  # scp '_aa_orf_' files back into Xenomorph and run RPSBLAST 
  # Read output back into R and proceed with next steps
  --------------------------------------------------------------------------------
  
  # read in RPSBLAST output
orfs_aa_rps_CAA <- read_tsv("c_adamanteus_DNAtransp_rpsblast.auto.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))
orfs_aa_rps_CAZ <- read_tsv("c_adamanteus_DNAtransp_rpsblast.chrZ.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))
orfs_aa_rps_CAW <- read_tsv("c_adamanteus_DNAtransp_rpsblast.chrW.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))

# split rpsblast title to make usable
orfs_aa_rps_CAA <- orfs_aa_rps_CAA %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")
orfs_aa_rps_CAZ <- orfs_aa_rps_CAZ %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")
orfs_aa_rps_CAW <- orfs_aa_rps_CAW %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")

# determine presence/absence of RT and EN
orfs_aa_rps_rt_CAA <- orfs_aa_rps_CAA %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_CAA <- orfs_aa_rps_CAA %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

orfs_aa_rps_rt_CAZ <- orfs_aa_rps_CAZ %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_CAZ <- orfs_aa_rps_CAZ %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

orfs_aa_rps_rt_CAW <- orfs_aa_rps_CAW %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_CAW <- orfs_aa_rps_CAW %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

# give seqnames their original names without the orf ids
orfs_aa_rps_rt_CAA = orfs_aa_rps_rt_CAA %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_CAA = orfs_aa_rps_en_CAA %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

orfs_aa_rps_rt_CAZ = orfs_aa_rps_rt_CAZ %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_CAZ = orfs_aa_rps_en_CAZ %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

orfs_aa_rps_rt_CAW = orfs_aa_rps_rt_CAW %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_CAW = orfs_aa_rps_en_CAW %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

# select orfs with intact rt and en
intact_orfs_CAA <- orfs_aa_rps_rt_CAA %>% filter(seqnames %in% orfs_aa_rps_en_CAA$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_CAA = intact_orfs_CAA %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_CAA = intact_orfs_bed_CAA[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_CAA), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("c_adamanteus_auto", "_intact_DNAorfs.bed"))

intact_orfs_CAZ <- orfs_aa_rps_rt_CAZ %>% filter(seqnames %in% orfs_aa_rps_en_CAZ$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_CAZ = intact_orfs_CAZ %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_CAZ = intact_orfs_bed_CAZ[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_CAZ), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("c_adamanteus_chrZ", "_intact_DNAorfs.bed"))

intact_orfs_CAW <- orfs_aa_rps_rt_CAW %>% filter(seqnames %in% orfs_aa_rps_en_CAW$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_CAW = intact_orfs_CAW %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_CAW = intact_orfs_bed_CAW[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_CAW), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("c_adamanteus_chrW", "_intact_DNAorfs.bed"))

rm(list = ls())
-----------------------------------------------------------------------------------
  ### Crotalus viridis ###
  # read in DNAtransp sequences
raw_seq_viridis_auto <- readDNAStringSet("c_viridis_DNAtransp.auto.fasta")
raw_seq_viridis_chrZ <- readDNAStringSet("c_viridis_DNAtransp.chrZ.fasta")
raw_seq_viridis_chrW <- readDNAStringSet("c_viridis_DNAtransp.chrW.fasta")

# create pseudoranges tibble of DNAtransp
raw_tbl_viridis_auto <- tibble(seqnames = names(raw_seq_viridis_auto), start = 1, end = width(raw_seq_viridis_auto)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")
raw_tbl_viridis_chrZ <- tibble(seqnames = names(raw_seq_viridis_chrZ), start = 1, end = width(raw_seq_viridis_chrZ)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")
raw_tbl_viridis_chrW <- tibble(seqnames = names(raw_seq_viridis_chrW), start = 1, end = width(raw_seq_viridis_chrW)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")

# create table of DNAtransp seqs with names
seq_group_tbl_viridis_auto <- raw_tbl_viridis_auto %>% dplyr::select(seqnames)
seq_group_tbl_viridis_chrZ <- raw_tbl_viridis_chrZ %>% dplyr::select(seqnames)
seq_group_tbl_viridis_chrW <- raw_tbl_viridis_chrW %>% dplyr::select(seqnames)

# find orfs over 1000bp in sequences
#orfs_1000 <- ORFik::findORFs(raw_seq, startCodon = startDefinition(1), minimumLength = 1000) %>% as_tibble()

orfs_1000_CVA = findORFsFasta("c_viridis_DNAtransp.auto.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_CVA = as_tibble(orfs_1000_CVA)
orfs_1000_CVA$idv_name = paste0(orfs_1000_CVA$seqnames, "#orf", 1:nrow(orfs_1000_CVA))

orfs_1000_CVZ = findORFsFasta("c_viridis_DNAtransp.chrZ.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_CVZ = as_tibble(orfs_1000_CVZ)
orfs_1000_CVZ$idv_name = paste0(orfs_1000_CVZ$seqnames, "#orf", 1:nrow(orfs_1000_CVZ))

orfs_1000_CVW = findORFsFasta("c_viridis_DNAtransp.chrW.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_CVW = as_tibble(orfs_1000_CVW)
orfs_1000_CVW$idv_name = paste0(orfs_1000_CVW$seqnames, "#orf", 1:nrow(orfs_1000_CVW))

# make ranges object of orfs
orfs_1000_CVA_ranges <- GRanges(seqnames = orfs_1000_CVA$seqnames, ranges = IRanges(start = orfs_1000_CVA$start, end = orfs_1000_CVA$end))
orfs_1000_CVZ_ranges <- GRanges(seqnames = orfs_1000_CVZ$seqnames, ranges = IRanges(start = orfs_1000_CVZ$start, end = orfs_1000_CVZ$end))
orfs_1000_CVW_ranges <- GRanges(seqnames = orfs_1000_CVW$seqnames, ranges = IRanges(start = orfs_1000_CVW$start, end = orfs_1000_CVW$end))

# get seq of orfs and name orfs
orfs_seq_CVA <- Biostrings::getSeq(raw_seq_viridis_auto, orfs_1000_CVA_ranges)
names(orfs_seq_CVA) <- orfs_1000_CVA$idv_name

orfs_seq_CVZ <- Biostrings::getSeq(raw_seq_viridis_chrZ, orfs_1000_CVZ_ranges)
names(orfs_seq_CVZ) <- orfs_1000_CVZ$idv_name

orfs_seq_CVW <- Biostrings::getSeq(raw_seq_viridis_chrW, orfs_1000_CVW_ranges)
names(orfs_seq_CVW) <- orfs_1000_CVW$idv_name

# translate orf
orfs_aa_seq_CVA <- translate(orfs_seq_CVA, if.fuzzy.codon = "solve")
names(orfs_aa_seq_CVA) <- orfs_1000_CVA$idv_name

orfs_aa_seq_CVZ <- translate(orfs_seq_CVZ, if.fuzzy.codon = "solve")
names(orfs_aa_seq_CVZ) <- orfs_1000_CVZ$idv_name

orfs_aa_seq_CVW <- translate(orfs_seq_CVW, if.fuzzy.codon = "solve")
names(orfs_aa_seq_CVW) <- orfs_1000_CVW$idv_name

# write nt and nt orfs to file
writeXStringSet(orfs_aa_seq_CVA, paste0("c_viridis_DNAtransp.auto.fasta", "_aa_orfs.auto.fa"))
writeXStringSet(orfs_seq_CVA, paste0("c_viridis_DNAtransp.auto.fasta", "_nt_orfs.auto.fa"))
file.rename(from = "c_viridis_DNAtransp.auto.fasta_aa_orfs.auto.fa", to = "c_viridis_DNAtransp_aa_orfs.auto.fa")
file.rename(from = "c_viridis_DNAtransp.auto.fasta_nt_orfs.auto.fa", to = "c_viridis_DNAtransp_nt_orfs.auto.fa")

writeXStringSet(orfs_aa_seq_CVZ, paste0("c_viridis_DNAtransp.chrZ.fasta", "_aa_orfs.chrZ.fa"))
writeXStringSet(orfs_seq_CVZ, paste0("c_viridis_DNAtransp.chrZ.fasta", "_nt_orfs.chrZ.fa"))
file.rename(from = "c_viridis_DNAtransp.chrZ.fasta_aa_orfs.chrZ.fa", to = "c_viridis_DNAtransp_aa_orfs.chrZ.fa")
file.rename(from = "c_viridis_DNAtransp.chrZ.fasta_nt_orfs.chrZ.fa", to = "c_viridis_DNAtransp_nt_orfs.chrZ.fa")

writeXStringSet(orfs_aa_seq_CVW, paste0("c_viridis_DNAtransp.chrW.fasta", "_aa_orfs.chrW.fa"))
writeXStringSet(orfs_seq_CVW, paste0("c_viridis_DNAtransp.chrW.fasta", "_nt_orfs.chrW.fa"))
file.rename(from = "c_viridis_DNAtransp.chrW.fasta_aa_orfs.chrW.fa", to = "c_viridis_DNAtransp_aa_orfs.chrW.fa")
file.rename(from = "c_viridis_DNAtransp.chrW.fasta_nt_orfs.chrW.fa", to = "c_viridis_DNAtransp_nt_orfs.chrW.fa")

--------------------------------------------------------------------------------
  # scp '_aa_orf_' files back into Xenomorph and run RPSBLAST 
  # Read output back into R and proceed with next steps
  --------------------------------------------------------------------------------
  
  # read in RPSBLAST output
orfs_aa_rps_CVA <- read_tsv("c_viridis_DNAtransp_rpsblast.auto.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))
# no hits in auto
orfs_aa_rps_CVZ <- read_tsv("c_viridis_DNAtransp_rpsblast.chrZ.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))
# no hits in Z
orfs_aa_rps_CVW <- read_tsv("c_viridis_DNAtransp_rpsblast.chrW.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))

# split rpsblast title to make usable
orfs_aa_rps_CVA <- orfs_aa_rps_CVA %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")
orfs_aa_rps_CVZ <- orfs_aa_rps_CVZ %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")
orfs_aa_rps_CVW <- orfs_aa_rps_CVW %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")

# determine presence/absence of RT and EN
orfs_aa_rps_rt_CVA <- orfs_aa_rps_CVA %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_CVA <- orfs_aa_rps_CVA %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

orfs_aa_rps_rt_CVZ <- orfs_aa_rps_CVZ %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_CVZ <- orfs_aa_rps_CVZ %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

orfs_aa_rps_rt_CVW <- orfs_aa_rps_CVW %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_CVW <- orfs_aa_rps_CVW %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

# give seqnames their original names without the orf ids
orfs_aa_rps_rt_CVA = orfs_aa_rps_rt_CVA %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_CVA = orfs_aa_rps_en_CVA %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

orfs_aa_rps_rt_CVZ = orfs_aa_rps_rt_CVZ %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_CVZ = orfs_aa_rps_en_CVZ %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

orfs_aa_rps_rt_CVW = orfs_aa_rps_rt_CVW %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_CVW = orfs_aa_rps_en_CVW %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

# select orfs with intact rt and en
intact_orfs_CVA <- orfs_aa_rps_rt_CVA %>% filter(seqnames %in% orfs_aa_rps_en_CVA$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_CVA = intact_orfs_CVA %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_CVA = intact_orfs_bed_CVA[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_CVA), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("c_viridis_auto", "_intact_DNAorfs.bed"))

intact_orfs_CVZ <- orfs_aa_rps_rt_CVZ %>% filter(seqnames %in% orfs_aa_rps_en_CVZ$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_CVZ = intact_orfs_CVZ %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_CVZ = intact_orfs_bed_CVZ[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_CVZ), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("c_viridis_chrZ", "_intact_DNAorfs.bed"))

intact_orfs_CVW <- orfs_aa_rps_rt_CVW %>% filter(seqnames %in% orfs_aa_rps_en_CVW$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_CVW = intact_orfs_CVW %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_CVW = intact_orfs_bed_CVW[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_CVW), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("c_viridis_chrW", "_intact_DNAorfs.bed"))

rm(list = ls())
-----------------------------------------------------------------------------------
  ### Crotalus horridus ###
  # read in DNAtransp sequences
raw_seq_horridus_auto <- readDNAStringSet("c_horridus_DNAtransp.auto.fasta")
raw_seq_horridus_chrZ <- readDNAStringSet("c_horridus_DNAtransp.chrZ.fasta")
raw_seq_horridus_chrW <- readDNAStringSet("c_horridus_DNAtransp.chrW.fasta")

# create pseudoranges tibble of DNAtransp
raw_tbl_horridus_auto <- tibble(seqnames = names(raw_seq_horridus_auto), start = 1, end = width(raw_seq_horridus_auto)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")
raw_tbl_horridus_chrZ <- tibble(seqnames = names(raw_seq_horridus_chrZ), start = 1, end = width(raw_seq_horridus_chrZ)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")
raw_tbl_horridus_chrW <- tibble(seqnames = names(raw_seq_horridus_chrW), start = 1, end = width(raw_seq_horridus_chrW)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")

# create table of DNAtransp seqs with names
seq_group_tbl_horridus_auto <- raw_tbl_horridus_auto %>% dplyr::select(seqnames)
seq_group_tbl_horridus_chrZ <- raw_tbl_horridus_chrZ %>% dplyr::select(seqnames)
seq_group_tbl_horridus_chrW <- raw_tbl_horridus_chrW %>% dplyr::select(seqnames)

# find orfs over 1000bp in sequences
#orfs_1000 <- ORFik::findORFs(raw_seq, startCodon = startDefinition(1), minimumLength = 1000) %>% as_tibble()

orfs_1000_CHA = findORFsFasta("c_horridus_DNAtransp.auto.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_CHA = as_tibble(orfs_1000_CHA)
orfs_1000_CHA$idv_name = paste0(orfs_1000_CHA$seqnames, "#orf", 1:nrow(orfs_1000_CHA))

orfs_1000_CHZ = findORFsFasta("c_horridus_DNAtransp.chrZ.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_CHZ = as_tibble(orfs_1000_CHZ)
orfs_1000_CHZ$idv_name = paste0(orfs_1000_CHZ$seqnames, "#orf", 1:nrow(orfs_1000_CHZ))

orfs_1000_CHW = findORFsFasta("c_horridus_DNAtransp.chrW.fasta", startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)
orfs_1000_CHW = as_tibble(orfs_1000_CHW)
orfs_1000_CHW$idv_name = paste0(orfs_1000_CHW$seqnames, "#orf", 1:nrow(orfs_1000_CHW))

# make ranges object of orfs
orfs_1000_CHA_ranges <- GRanges(seqnames = orfs_1000_CHA$seqnames, ranges = IRanges(start = orfs_1000_CHA$start, end = orfs_1000_CHA$end))
orfs_1000_CHZ_ranges <- GRanges(seqnames = orfs_1000_CHZ$seqnames, ranges = IRanges(start = orfs_1000_CHZ$start, end = orfs_1000_CHZ$end))
orfs_1000_CHW_ranges <- GRanges(seqnames = orfs_1000_CHW$seqnames, ranges = IRanges(start = orfs_1000_CHW$start, end = orfs_1000_CHW$end))

# get seq of orfs and name orfs
orfs_seq_CHA <- Biostrings::getSeq(raw_seq_horridus_auto, orfs_1000_CHA_ranges)
names(orfs_seq_CHA) <- orfs_1000_CHA$idv_name

orfs_seq_CHZ <- Biostrings::getSeq(raw_seq_horridus_chrZ, orfs_1000_CHZ_ranges)
names(orfs_seq_CHZ) <- orfs_1000_CHZ$idv_name

orfs_seq_CHW <- Biostrings::getSeq(raw_seq_horridus_chrW, orfs_1000_CHW_ranges)
names(orfs_seq_CHW) <- orfs_1000_CHW$idv_name

# translate orf
orfs_aa_seq_CHA <- translate(orfs_seq_CHA, if.fuzzy.codon = "solve")
names(orfs_aa_seq_CHA) <- orfs_1000_CHA$idv_name

orfs_aa_seq_CHZ <- translate(orfs_seq_CHZ, if.fuzzy.codon = "solve")
names(orfs_aa_seq_CHZ) <- orfs_1000_CHZ$idv_name

orfs_aa_seq_CHW <- translate(orfs_seq_CHW, if.fuzzy.codon = "solve")
names(orfs_aa_seq_CHW) <- orfs_1000_CHW$idv_name

# write nt and nt orfs to file
writeXStringSet(orfs_aa_seq_CHA, paste0("c_horridus_DNAtransp.auto.fasta", "_aa_orfs.auto.fa"))
writeXStringSet(orfs_seq_CHA, paste0("c_horridus_DNAtransp.auto.fasta", "_nt_orfs.auto.fa"))
file.rename(from = "c_horridus_DNAtransp.auto.fasta_aa_orfs.auto.fa", to = "c_horridus_DNAtransp_aa_orfs.auto.fa")
file.rename(from = "c_horridus_DNAtransp.auto.fasta_nt_orfs.auto.fa", to = "c_horridus_DNAtransp_nt_orfs.auto.fa")

writeXStringSet(orfs_aa_seq_CHZ, paste0("c_horridus_DNAtransp.chrZ.fasta", "_aa_orfs.chrZ.fa"))
writeXStringSet(orfs_seq_CHZ, paste0("c_horridus_DNAtransp.chrZ.fasta", "_nt_orfs.chrZ.fa"))
file.rename(from = "c_horridus_DNAtransp.chrZ.fasta_aa_orfs.chrZ.fa", to = "c_horridus_DNAtransp_aa_orfs.chrZ.fa")
file.rename(from = "c_horridus_DNAtransp.chrZ.fasta_nt_orfs.chrZ.fa", to = "c_horridus_DNAtransp_nt_orfs.chrZ.fa")

writeXStringSet(orfs_aa_seq_CHW, paste0("c_horridus_DNAtransp.chrW.fasta", "_aa_orfs.chrW.fa"))
writeXStringSet(orfs_seq_CHW, paste0("c_horridus_DNAtransp.chrW.fasta", "_nt_orfs.chrW.fa"))
file.rename(from = "c_horridus_DNAtransp.chrW.fasta_aa_orfs.chrW.fa", to = "c_horridus_DNAtransp_aa_orfs.chrW.fa")
file.rename(from = "c_horridus_DNAtransp.chrW.fasta_nt_orfs.chrW.fa", to = "c_horridus_DNAtransp_nt_orfs.chrW.fa")

--------------------------------------------------------------------------------
  # scp '_aa_orf_' files back into Xenomorph and run RPSBLAST 
  # Read output back into R and proceed with next steps
  --------------------------------------------------------------------------------
  
  # read in RPSBLAST output
orfs_aa_rps_CHA <- read_tsv("c_horridus_DNAtransp_rpsblast.auto.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))
orfs_aa_rps_CHZ <- read_tsv("c_horridus_DNAtransp_rpsblast.chrZ.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))
orfs_aa_rps_CHW <- read_tsv("c_horridus_DNAtransp_rpsblast.chrW.tsv", col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))

# split rpsblast title to make usable
orfs_aa_rps_CHA <- orfs_aa_rps_CHA %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")
orfs_aa_rps_CHZ <- orfs_aa_rps_CHZ %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")
orfs_aa_rps_CHW <- orfs_aa_rps_CHW %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")

# determine presence/absence of RT and EN
orfs_aa_rps_rt_CHA <- orfs_aa_rps_CHA %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_CHA <- orfs_aa_rps_CHA %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

orfs_aa_rps_rt_CHZ <- orfs_aa_rps_CHZ %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_CHZ <- orfs_aa_rps_CHZ %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

orfs_aa_rps_rt_CHW <- orfs_aa_rps_CHW %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en_CHW <- orfs_aa_rps_CHW %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

# give seqnames their original names without the orf ids
orfs_aa_rps_rt_CHA = orfs_aa_rps_rt_CHA %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_CHA = orfs_aa_rps_en_CHA %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

orfs_aa_rps_rt_CHZ = orfs_aa_rps_rt_CHZ %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_CHZ = orfs_aa_rps_en_CHZ %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

orfs_aa_rps_rt_CHW = orfs_aa_rps_rt_CHW %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en_CHW = orfs_aa_rps_en_CHW %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

# select orfs with intact rt and en
intact_orfs_CHA <- orfs_aa_rps_rt_CHA %>% filter(seqnames %in% orfs_aa_rps_en_CHA$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_CHA = intact_orfs_CHA %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_CHA = intact_orfs_bed_CHA[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_CHA), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("c_horridus_auto", "_intact_DNAorfs.bed"))

intact_orfs_CHZ <- orfs_aa_rps_rt_CHZ %>% filter(seqnames %in% orfs_aa_rps_en_CHZ$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_CHZ = intact_orfs_CHZ %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_CHZ = intact_orfs_bed_CHZ[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_CHZ), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("c_horridus_chrZ", "_intact_DNAorfs.bed"))

intact_orfs_CHW <- orfs_aa_rps_rt_CHW %>% filter(seqnames %in% orfs_aa_rps_en_CHW$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

intact_orfs_bed_CHW = intact_orfs_CHW %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
intact_orfs_bed_CHW = intact_orfs_bed_CHW[,c(2:4,1)]
write.table(x = as.data.frame(intact_orfs_bed_CHW), sep = "\t", quote = F, col.names = T, row.names = F, file = paste0("c_horridus_chrW", "_intact_DNAorfs.bed"))
