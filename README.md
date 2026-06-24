# Metabarcoding sequencing from two Marine Rocky Shorelines.
Bioinformatic pipeline for metabarcoding sequencing.
# 1) Trimming and quality control
## 1.CUTADAPT
Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10–12. DOI: 10.14806/ej.17.1.200. 
```
##########################################################################
  1. Use Cutadapt to select the sequences that contains the target primers
##########################################################################


#!/bin/bash

# Forward and reverse primers
############for 16S###########
PRIMERF="AGAGTTTGATCMTGGCTCAG"
PRIMERR="RGYTACCTTGTTACGACTT"
PRIMERFC="CTGAGCCAKGATCAAACTCT"
PRIMERRC="AAGTCGTAACAAGGTARCY"
############for 16S###########

############for 18S###########
PRIMERF="GAAACTGCGAATGGCTC"
PRIMERR="CYGCAGGTTCACCTAC"
PRIMERFC="GAGCCATTCGCAGTTTC"
PRIMERRC="GTAGGTGAACCTGCRG"
############for 18S###########

# Loop over all FASTQ.gz files in the current directory
for file in *.fastq.gz; do
if [[ -f "$file" ]]; then
base=$(basename "$file" .fastq.gz)
echo "Trimming: $file"

cutadapt \
-g "$PRIMERF...$PRIMERRC" \
-g "$PRIMERR...$PRIMERFC" \
-o "${base}_trimmed.fastq.gz" \
"$file" \
--error-rate 0.2 \
--cores=40 \
--discard-untrimmed
fi
done
```

## 2.NANOFILT
De Coster, W., D’Hert, S., Schultz, D. T., Cruts, M., Van Broeckhoven, C. (2018) NanoPack: visualizing and processing long-read sequencing data, Bioinformatics, Volume 34, Issue 15, 2666–2669, https://doi.org/10.1093/bioinformatics/bty149.
```
########################################################
  2.Use NanoFilt to quality trimm the selected sequences
########################################################

#!/bin/bash

for file in *.fastq.gz; do
if [[ -f "$file" ]]; then
base=$(basename "$file" _trimmed.fastq.gz)
echo "Trimming: $file"
gunzip -c "$file" | NanoFilt -q 15 -l xxx --maxlength xxxx --tailcrop 50 --headcrop 50 > "${base}_filtered_crop_q15.fastq"
fi
done

#################################
-l 800 --maxlength 1600 for 16S
-l 800 --maxlength 1900 for 18S
#################################
```
# 2) Select representative sequences, taxonomic assignemnet and OTU table construction
##  1.VSEARCH
Rognes T, Flouri T, Nichols B, Quince C, Mahé F. VSEARCH: a versatile open-source tool for metagenomics. PeerJ. 2016 Oct 18;4:e2584.
```
#############################################################
1. Use VSEARCH to orient reads and to remove the chimeras
############################################################

vsearch --orient filtered_q15.fa --db ref_db.fasta --fastaout oriented_q15.fa
##############################
ref_db=silva_16S.fasta for 16S
ref_db=silva_18S.fasta for 18S
##############################

```
```
vsearch --uchime_ref input.fasta \
        --nonchimeras output.good.fasta \
        --chimeras output.chimeras.fasta \
        --db ref_db.fasta
##############################
ref_db=silva_16S.fasta for 16S
ref_db=silva_18S.fasta for 18S
##############################

```
## 2.Amplicon Sorter
Vierstraete, A. R., Braeckman, B. P. (2022). Amplicon_sorter: A tool for reference-free amplicon sorting based on sequence similarity and for building consensus sequences. Ecology and Evolution, 12, e8603. https://doi.org/10.1002/ece3.8603.
```
###########################################################################################################
2.Use Amplicon Sorter to cluster filtered reads into groups of closely related species and generate consensus sequences.
###########################################################################################################

./amplicon_sorter.py   -i all_oriented.fasta -o amplicon_clustered_oriented -np 64   -min xxx   -max xxxx -maxr=xxxxx --allreads
########################################
-min 800 --max 1600 for 16S -maxr=482656
-min 800 --max 1900 for 18S -maxr=285075
########################################
```
## 3.VSEARCH taxonomy 
```
##################################################################################################
3. Use VSERACH to assign taxonomy to the consensus sequences against a reference database (SILVA_nr99)
###################################################################################################

vsearch -usearch_global otus.fa -db ref_db\
  -id 0.90 \
  --output_no_hits \ #### keep the nohit in the output #####
  -blast6out taxonomy_assignment_silva.txt
##############################
ref_db=silva_16S.fasta for 16S
ref_db=silva_18S.fasta for 18S
##############################

```
## 4.Create an OTU table  
```
#####################################
3. Use VSERACH to create an OTU table 
#####################################
1.Rename the fasta header of each fasta file >barcode01_1;Sample=barcode01

#!/bin/bash
# Step 1: rename sequences and append sample=barcode to each header

for f in *_combined_filtered_crop_q15.fasta; do
    sample=$(basename "$f" "_combined_filtered_crop_q15.fasta")
    renamed="${sample}_renamed.fasta"

    # 1. Rename: >SAMPLE_1, SAMPLE_2, ...
    awk -v prefix="$sample" '/^>/ {i++; print ">"prefix"_"i} !/^>/' "$f" > "$renamed"

    # 2. Append ;sample=barcode to each header
    awk -v barcode="$sample" '/^>/{sub(/$/, ";sample=" barcode)} 1' "$renamed" > "${sample}_renamed_append.fasta"

    echo "Processed: $f -> ${sample}_renamed_append.fasta"
done
```
2.Map reads to reference FASTA
```
vsearch -usearch_global all_reads.fasta -db consensus_ref_seqs.fasta -strand both -id 0.90 -otutabout otu.tsv
```

The taxonomy and OTU table files can now be imported into R or converted into a BIOM file for use in other pipelines.

# Publication
Microbial community and biodiversity meta-barcoding sequencing data from rocky shorelines on the northeast and southwest coasts of England.  
