# setup the abricate database and run abricate
cp reference.fa /home/jason/miniconda3/envs/python3.6/db/test/sequences
abricate --setupdb
abricate --db test --threads 18 contigs.fasta > results.tsv

# target genes extraction
cut -f 1-5 results.tsv > test.tsv

while read A B C D E
do if  
[[ "$E" == "+" ]]
then seqkit grep -r -p $B $A | seqkit subseq -r $C:$D > ../fna/arnF/"${A##*/}"
else seqkit grep -r -p $B $A | seqkit subseq -r $C:$D | seqkit seq -r -p > ../fna/arnF/"${A##*/}"
fi
done < test.tsv