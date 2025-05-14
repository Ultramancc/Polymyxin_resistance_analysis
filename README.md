# Polymyxin_resistance (PR) analysis #
This repository outlines the analysis pipeline for PR mechanism and genomic diversity in _Klebsiella pneumoniae_.  

## PR related gene screening ##
Set up the abricate database and run abricate (https://github.com/tseemann/abricate). The PR-related genes was extracted from polymyxin susceptible genome _K. pneumoniae_ MGH 78578 (GenBank accession number: CP000647).  
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

## Mutation site output ##
All extracted PR genes were aligned using muscle (https://github.com/rcedgar/muscle) and output the site mutation for each genome.  
```
python PR_mutation.py -h
usage: PR_mutation.py [-h] -i INPUT -o OUTPUT -p PROTEIN -t TABLE

Translate nucleotide sequences to proteins, align them, and find mutations.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input nucleotide FASTA file (e.g., sequence.all.fna)
  -o OUTPUT, --output OUTPUT
                        Output alignment file (e.g., protein.aln.fna)
  -p PROTEIN, --protein PROTEIN
                        Output translated protein file (e.g., protein.all.fna)
  -t TABLE, --table TABLE
                        Output mutation table file (e.g., mutations.tsv)
```
Example usage
```
python PR_mutation.py -i PR_gene.all.fna -o PR_gene.pro.aln -p PR_gene.pro.fna -t PR_gene.mutation.tsv
```

## IS disruption detection ##
