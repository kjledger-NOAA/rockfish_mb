# rockfish_mb

step 1: transfer all raw sequencing reads into /genetics/edna/workdir/sebastes/tissues/

step 2: remove primer sequences using cutadapt 
- activated cutadapt env (conda activate cutadaptenv)
- set DATA to be directory with files (DATA=/genetics/edna/workdir/sebastes/tissues)
 - create array of file names:  NAMELIST=$(ls ${DATA} | sed 's/e*_L001.*//' | uniq)
 - double check: echo "${NAMELIST}"
- make folder for trimmed reads
- do: for i in ${NAMELIST}; do cutadapt --discard-untrimmed -g ATNACCATATCTAGGNTTNAACC -G TGRRCTTGTTGGTCGGYT -o trimmed/${i}_R1.fastq.gz -p trimmed/${i}_R2.fastq.gz "$DATA/${i}_L001_R1_001.fastq.gz" "$DATA/${i}_L001_R2_001.fastq.gz"; done
- unzip trimmed reads: pigz -d trimmed/*.gz

step 3: use DADA2 to filter and merge sequencing reads 
- make folder for filtered reads (mkdir ../filtered)
- run sequence_filtering.Rmd 

step 4: use blastn of custom rockfish db (v.519) for taxonomic assignment 
blastn -query /genetics/edna/workdir/sebastes/tissues/trimmed/filtered/outputs/myasvs.fasta -db /home/kimberly.ledger/rockfish_mb/custom_db/rockfish_db_519 -out blastn_custom519_20241106.txt -perc_identity 96 -qcov_hsp_perc 98 -num_threads 10 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames staxids'

** still need to make final adjustments to custom reference db - there are missing species 

step 5: use "rockfishdb_taxonomic_assignment.Rmd" to figure out ASV assignments 

step 6: compare taxonomic ID to known rockfish tissue ID using "ablg_tissue_test.Rmd"