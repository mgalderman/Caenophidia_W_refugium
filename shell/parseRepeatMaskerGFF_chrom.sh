list=$1
gff=$2
for chrom in `cat $list`; do
	grep -w $chrom $gff
done