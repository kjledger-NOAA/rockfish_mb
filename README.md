# rockfish_mb

step 1: combined all raw sequencing reads (tissue, mock communities, and eDNA samples) into a single folder and created sample sheet. uploaded folder to sedna. 

step 2: processed samples using dadasnake (config.rkfish.yaml) and saved output in 'rkfish_tissue_mock_bt22'

step 3: use blastn of custom rockfish db (v.536) for taxonomic assignment

v.536 has the addition of dark and dusky ref seqs and the subtraction of a few ref seqs that could not map primers or tree placement suggested mis-id's 
blastn -query /home/kimberly.ledger/rockfish_mb/data/dadasnake_output/filtered.seqs.fasta -db /home/kimberly.ledger/rockfish_mb/custom_db/rockfish_db_536 -out blastn_tissue_mock_bt22_db536.txt -perc_identity 92 -qcov_hsp_perc 98 -num_threads 10 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames staxids'

note: if the custom rockfish db needs to be updated.
- export fasta and metadata csv (Name and Organism) from Geneious folder of rockfish
- example: makeblastdb -in rockfish_reference_db_529.fasta -dbtype nucl -out rockfish_db_529

step 4: use "1_rockfishdb_taxonomic_assignment.Rmd" to figure out ASV assignments 

step 5: run "2_decontamination.Rmd" 

step 6: compare taxonomic ID to known rockfish tissue ID using "ablg_tissue_test.Rmd"





to try and resolve the dusky, dark, northern group, i made a haplotype network and then choose representative sequences for each haplotype for another custom db. this was called cpv_db 

blastn -query /home/kimberly.ledger/rockfish_mb/data/dadasnake_output/filtered.seqs.fasta -db /home/kimberly.ledger/rockfish_mb/custom_db/cpv_db -out blastn_tissue_mock_bt22_cpv.txt -perc_identity 92 -qcov_hsp_perc 98 -num_threads 10 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames staxids'

after doing this i removed the few intial ciliatus and variabilis seq i added to the ref db and replaced them with the haplotypes 

the new rfdb is v 534_20250117

blastn -query /home/kimberly.ledger/rockfish_mb/data/dadasnake_output/filtered.seqs.fasta -db /home/kimberly.ledger/rockfish_mb/custom_db/rockfish_db_534_20250117 -out blastn_tissue_mock_bt22_534_20250117.txt -perc_identity 92 -qcov_hsp_perc 98 -num_threads 10 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames staxids'