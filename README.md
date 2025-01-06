# rockfish_mb

step 1: transfer all raw sequencing reads into /genetics/edna/workdir/sebastes/tissues/
12/31/24 update - added tissue sequences from 12/19/24 sequencing run 

step 2: FIND primer sequences using cutadapt (not removing as of 12/31/24 updates)
- activated cutadapt env (conda activate cutadaptenv)
- set DATA to be directory with files (DATA=/genetics/edna/workdir/sebastes/tissues)
 - create array of file names:  NAMELIST=$(ls ${DATA} | sed 's/e*_L001.*//' | uniq)
 - double check: echo "${NAMELIST}"
- make folder for untrimmed reads with primer regions 
- do: for i in ${NAMELIST}; do cutadapt --discard-untrimmed --no-trim -g ATNACCATATCTAGGNTTNAACC -G TGRRCTTGTTGGTCGGYT -o untrimmed_20241231/${i}_R1.fastq.gz -p untrimmed_20241231/${i}_R2.fastq.gz "$DATA/${i}_L001_R1_001.fastq.gz" "$DATA/${i}_L001_R2_001.fastq.gz"; done
- unzip untrimmed reads: pigz -d untrimmed_20241231/*.gz

step 3: use DADA2 to filter and merge sequencing reads 
- make folder for filtered reads (mkdir ../filtered)
- run sequence_filtering.Rmd 

step 4: use blastn of custom rockfish db (v.536) for taxonomic assignment (updated after 12/2/2024 miseq run) 

v.536 has the addition of dark and dusky ref seqs and the subtraction of a few ref seqs that could not map primers or tree placement suggested mis-id's 

(base) [kimberly.ledger@akc0ss-vu-134 data]$ blastn -query /genetics/edna/workdir/sebastes/tissues/untrimmed_20241231/filtered/outputs/myasvs.fasta -db /home/kimberly.ledger/rockfish_mb/custom_db/rockfish_db_536 -out blastn_custom536_20250103.txt -perc_identity 92 -qcov_hsp_perc 98 -num_threads 10 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames staxids'

## also running it on the version of the sequencing data where the primers were removed. 
(base) [kimberly.ledger@akc0ss-vu-134 data]$ blastn -query /genetics/edna/workdir/sebastes/tissues/trimmed_20241204/filtered/outputs/myasvs.fasta -db /home/kimberly.ledger/rockfish_mb/custom_db/rockfish_db_536 -out blastn_custom536_20250103_trimmed.txt -perc_identity 92 -qcov_hsp_perc 98 -num_threads 10 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames staxids'


note: if the custom rockfish db needs to be updated.
- export fasta and metadata csv (Name and Organism) from Geneious folder of rockfish
-example: makeblastdb -in rockfish_reference_db_529.fasta -dbtype nucl -out rockfish_db_529

step 5: use "rockfishdb_taxonomic_assignment.Rmd" to figure out ASV assignments 

step 6: compare taxonomic ID to known rockfish tissue ID using "ablg_tissue_test.Rmd"