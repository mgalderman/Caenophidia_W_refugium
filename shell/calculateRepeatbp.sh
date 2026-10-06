gff=$1
awk 'BEGIN{OFS="\t"}{print $1,$4-1,$5}' $gff | bedtools sort | bedtools merge > tmp.bed.merge
awk '{print $3-$2}' tmp.bed.merge > tmp.bed.lengths
total=$(awk '{sum+=$1;}END{print sum;}' tmp.bed.lengths)
echo Total length of annotated repeats in $gff is $total bp
rm tmp.bed.merge
rm tmp.bed.lengths