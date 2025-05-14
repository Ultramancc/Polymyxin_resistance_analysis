# Polymyxin_resistance (PR) analysis #
This repository outlines the analysis pipeline for PR mechanism and genomic diversity in _Klebsiella pneumoniae_.  

## PR related gene screening ##
Set up the abricate database and run abricate (https://github.com/tseemann/abricate)  
```
cp PR-related_nucl.fna /home/jason/miniconda3/envs/python3.6/db/PR_screen/sequences
abricate --setupdb
abricate --db PR_screen --threads 10 genome_fna/*.fna > PR_gene.tsv
```
Gene extraction
```
cut -f 1-5 PR_gene.tsv > list.tsv

while read A B C D E
do if  
[[ "$E" == "+" ]]
then seqkit grep -r -p $B $A | seqkit subseq -r $C:$D > output/"${A##*/}"
else seqkit grep -r -p $B $A | seqkit subseq -r $C:$D | seqkit seq -r -p > output/"${A##*/}"
fi
done < list.tsv
```
