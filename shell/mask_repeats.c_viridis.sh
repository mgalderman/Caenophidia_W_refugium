
#!/bin/bash

cd ./c_viridis

mkdir -p 1_bovb_cr1_mask 2_repbase_tetrapoda_mask 3_snake_known_mask 4_snake_unknown_mask 5_full_mask

# run BovB/CR1 masking
~/tmp/repeat-annotation/RepeatMasker/RepeatMasker -pa 8 -engine ncbi -lib ../libraries/CR1_BovB_Squamates_TElib.fasta -a -dir 1_bovb_cr1_mask /media/queen/extradrive2/caenophidia_W_refugium/genomes/c_viridis/c_viridis_final.fasta 2>&1 | tee 1_bovb_cr1_mask/1_bovb_cr1_mask.log

rename 's/fasta/bovb_mask/g' 1_bovb_cr1_mask/*
rename 's/.masked$/.masked.fasta/g' 1_bovb_cr1_mask/*

# run RepBase Tetrapoda masking based on repease 20181026
~/tmp/repeat-annotation/RepeatMasker/RepeatMasker -pa 8 -engine ncbi -species tetrapoda -a -dir 2_repbase_tetrapoda_mask 1_bovb_cr1_mask/c_viridis_final.bovb_mask.masked.fasta 2>&1 | tee 2_repbase_tetrapoda_mask/2_repbase_tetrapoda_mask.log

rename 's/bovb_mask.masked.fasta/tetrapoda_mask/g' 2_repbase_tetrapoda_mask/*
rename 's/masked$/masked.fasta/g' 2_repbase_tetrapoda_mask/*

# run "allsnake2025" known masking
~/tmp/repeat-annotation/RepeatMasker/RepeatMasker -pa 8 -engine ncbi -lib ../libraries/allsnakes2025.known.clust.fasta -a -dir 3_snake_known_mask 2_repbase_tetrapoda_mask/c_viridis_final.tetrapoda_mask.masked.fasta 2>&1 | tee 3_snake_known_mask/3_snake_known_mask.log

rename 's/tetrapoda_mask.masked.fasta/snake_known_mask/g' 3_snake_known_mask/*
rename 's/masked$/masked.fasta/g' 3_snake_known_mask/*

# run "allsnake2025" unknown masking
~/tmp/repeat-annotation/RepeatMasker/RepeatMasker -pa 8 -engine ncbi -lib ../libraries/allsnakes2025.unknown.clust.fasta -a -dir 4_snake_unknown_mask 3_snake_known_mask/c_viridis_final.snake_known_mask.masked.fasta 2>&1 | tee 4_snake_unknown_mask/4_snake_unknown_mask.log

rename 's/snake_known_mask.masked.fasta/snake_unknown_mask/g' 4_snake_unknown_mask/*
rename 's/masked$/masked.fasta/g' 4_snake_unknown_mask/*

# summarize/combine full output
# .fasta
cp 4_snake_unknown_mask/c_viridis_final.snake_unknown_mask.masked.fasta 5_full_mask/c_viridis_final.full_mask.masked.fasta

# .out
cat <(cat 1_bovb_cr1_mask/c_viridis_final.bovb_mask.out) <(cat 2_repbase_tetrapoda_mask/c_viridis_final.tetrapoda_mask.out | tail -n +4) <(cat 3_snake_known_mask/c_viridis_final.snake_known_mask.out | tail -n +4) <(cat 4_snake_unknown_mask/c_viridis_final.snake_unknown_mask.out | tail -n +4) > 5_full_mask/c_viridis_final.full_mask.out

# .cat
cat 1_bovb_cr1_mask/c_viridis_final.bovb_mask.cat.gz 2_repbase_tetrapoda_mask/c_viridis_final.tetrapoda_mask.cat.gz 3_snake_known_mask/c_viridis_final.snake_known_mask.cat.gz 4_snake_unknown_mask/c_viridis_final.snake_unknown_mask.cat.gz > 5_full_mask/c_viridis_final.full_mask.cat.gz

# .align
cat 1_bovb_cr1_mask/c_viridis_final.bovb_mask.align 2_repbase_tetrapoda_mask/c_viridis_final.tetrapoda_mask.align 3_snake_known_mask/c_viridis_final.snake_known_mask.align 4_snake_unknown_mask/c_viridis_final.snake_unknown_mask.align > 5_full_mask/c_viridis_final.full_mask.align

# .tbl
~/tmp/repeat-annotation/RepeatMasker/ProcessRepeats -species tetrapoda 5_full_mask/c_viridis_final.full_mask.cat.gz

# create repeat landscape files
~/tmp/repeat-annotation/faToTwoBit /media/queen/extradrive2/caenophidia_W_refugium/genomes/c_viridis/c_viridis_final.fasta ./c_viridis_final.2bit
~/tmp/repeat-annotation/RepeatMasker/util/calcDivergenceFromAlign.pl -s 5_full_mask/c_viridis_final.full_mask.landscape 5_full_mask/c_viridis_final.full_mask.align
~/tmp/repeat-annotation/RepeatMasker/util/createRepeatLandscape.pl -div 5_full_mask/c_viridis_final.full_mask.landscape -twoBit ./c_viridis_final.2bit > 5_full_mask/c_viridis_final.full_mask.landscape.html

# create GFFs
~/tmp/repeat-annotation/RepeatMasker/util/rmOutToGFF3.pl 5_full_mask/c_viridis_final.full_mask.out > 5_full_mask/c_viridis_final.full_mask.gff3
cat 5_full_mask/c_viridis_final.full_mask.gff3 | perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' > 5_full_mask/c_viridis_final.full_mask.reformat.gff3

# filter out simple repeats
grep -v -e "Satellite" -e ")n" -e "-rich" 5_full_mask/c_viridis_final.full_mask.reformat.gff3 > 5_full_mask/c_viridis_final.full_mask.reformat.complex.gff3

# compress outputs
gzip */*.out
gzip */*.fasta
gzip */*.align

cd ..
