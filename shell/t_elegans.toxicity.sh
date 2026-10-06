
#!/bin/bash

## Toxicity analysis of the T. elegans genome for fl-LTRs 

cp /media/queen/extradrive2/caenophidia_W_refugium/genomes/t_elegans/t_elegans_auto.LTR.fasta /media/queen/extradrive2/caenophidia_W_refugium/analysis/toxicity_analysis/LTRharvest
cp /media/queen/extradrive2/caenophidia_W_refugium/genomes/t_elegans/t_elegans_chrZ.LTR.fasta /media/queen/extradrive2/caenophidia_W_refugium/analysis/toxicity_analysis/LTRharvest
cp /media/queen/extradrive2/caenophidia_W_refugium/genomes/t_elegans/t_elegans_chrW.LTR.fasta /media/queen/extradrive2/caenophidia_W_refugium/analysis/toxicity_analysis/LTRharvest

echo "genomes acquired, running suffixerator" 

gt suffixerator -db t_elegans_auto.LTR.fasta -indexname t_elegans_auto.LTR.fasta -tis -suf -lcp -des -ssp -sds -dna
gt suffixerator -db t_elegans_chrZ.LTR.fasta -indexname t_elegans_chrZ.LTR.fasta -tis -suf -lcp -des -ssp -sds -dna
gt suffixerator -db t_elegans_chrW.LTR.fasta -indexname t_elegans_chrW.LTR.fasta -tis -suf -lcp -des -ssp -sds -dna

echo "suffixerator complete" 

gt ltrharvest -index t_elegans_auto.LTR.fasta -gff3 t_elegans_auto.LTR.fasta.gff -out t_elegans_auto.LTR.fasta.ltr.fa
gt ltrharvest -index t_elegans_auto.LTR.fasta -seqids yes -tabout no > t_elegans_auto.LTR.ltrharvest.out

gt ltrharvest -index t_elegans_chrZ.LTR.fasta -gff3 t_elegans_chrZ.LTR.fasta.gff -out t_elegans_chrZ.LTR.fasta.ltr.fa
gt ltrharvest -index t_elegans_chrZ.LTR.fasta -seqids yes -tabout no > t_elegans_chrZ.LTR.ltrharvest.out

gt ltrharvest -index t_elegans_chrW.LTR.fasta -gff3 t_elegans_chrW.LTR.fasta.gff -out t_elegans_chrW.LTR.fasta.ltr.fa
gt ltrharvest -index t_elegans_chrW.LTR.fasta -seqids yes -tabout no > t_elegans_chrW.LTR.ltrharvest.out

echo "LTRharvest complete!"

# Run LTRdigest to extract full-length LTR elements.

gt gff3 -sortlines yes -retainids yes -tidy yes -fixregionboundaries yes -checkids yes t_elegans_auto.LTR.fasta.gff > t_elegans_auto.LTR.fasta.sorted.gff
gt gff3 -sortlines yes -retainids yes -tidy yes -fixregionboundaries yes -checkids yes t_elegans_chrZ.LTR.fasta.gff > t_elegans_chrZ.LTR.fasta.sorted.gff
gt gff3 -sortlines yes -retainids yes -tidy yes -fixregionboundaries yes -checkids yes t_elegans_chrW.LTR.fasta.gff > t_elegans_chrW.LTR.fasta.sorted.gff

mv t_elegans_auto.LTR.fasta.sorted.gff t_elegans_auto.LTR.fasta.gff
mv t_elegans_chrZ.LTR.fasta.sorted.gff t_elegans_chrZ.LTR.fasta.gff
mv t_elegans_chrW.LTR.fasta.sorted.gff t_elegans_chrW.LTR.fasta.gff

gt ltrdigest -hmms ./hmm/*.hmm -aaout yes -outfileprefix t_elegans_auto.LTR.fasta_ltrdigest t_elegans_auto.LTR.fasta.gff t_elegans_auto.LTR.fasta > t_elegans_auto.LTR.fasta_ltrdigest_output_gff
gt ltrdigest -hmms ./hmm/*.hmm -aaout yes -outfileprefix t_elegans_chrZ.LTR.fasta_ltrdigest t_elegans_chrZ.LTR.fasta.gff t_elegans_chrZ.LTR.fasta > t_elegans_chrZ.LTR.fasta_ltrdigest_output_gff
gt ltrdigest -hmms ./hmm/*.hmm -aaout yes -outfileprefix t_elegans_chrW.LTR.fasta_ltrdigest t_elegans_chrW.LTR.fasta.gff t_elegans_chrW.LTR.fasta > t_elegans_chrW.LTR.fasta_ltrdigest_output_gff

echo "LTRdigest complete!"

# Move individual digest sequence files to the fasta subdirectory.

mv *.fas ./fasta
echo "fastas moved"

# Extract GFF entries for full-length LTRs.

echo "almost done!"
grep -v '#' t_elegans_auto.LTR.fasta_ltrdigest_output_gff | awk '$3 == "LTR_retrotransposon" {print}' > fl-LTR.t_elegans.auto.fasta.gff
grep -v '#' t_elegans_chrZ.LTR.fasta_ltrdigest_output_gff | awk '$3 == "LTR_retrotransposon" {print}' > fl-LTR.t_elegans.chrZ.fasta.gff
grep -v '#' t_elegans_chrW.LTR.fasta_ltrdigest_output_gff | awk '$3 == "LTR_retrotransposon" {print}' > fl-LTR.t_elegans.chrW.fasta.gff

rm t_elegans_auto.LTR.fasta
rm t_elegans_chrZ.LTR.fasta
rm t_elegans_chrW.LTR.fasta

mkdir t_elegans
mv t_elegans_* ./t_elegans
mv fl-LTR* ./t_elegans/

