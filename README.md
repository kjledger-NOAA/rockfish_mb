# Development and validation of an eDNA primer set for rockfishes in the Alaska Current provides a new tool for fisheries management

Kimberly J. Ledger<sup>1*</sup>, Meredith Everett<sup>2</sup>, Krista M. Nichols<sup>2</sup>, Wes Larson<sup>1</sup>, and Diana S. Baetscher<sup>1</sup>

<sup>1</sup> NOAA, National Marine Fisheries Service, Alaska Fisheries Science Center, Auke Bay Laboratories, Juneau, Alaska, USA <br />
<sup>2</sup> NOAA, National Marine Fisheries Service, Northwest Fisheries Science Center, Seattle, Washington, USA <br />
* Corresponding author: kimberly.ledger@noaa.gov <br />

## Overview 
This repo is dedicated to hosting data and code for the manuscript 

### 1. pre-processing 
step 1: combined all raw sequencing reads (tissue, mock communities, and eDNA samples) into a single folder and created sample sheet. uploaded folder to HPCC   

step 2: processed samples using dadasnake (config.rkfish.yaml) and saved output in 'rkfish_tissue_mock_bt22'

step 3: use blastn of custom rockfish db for taxonomic assignment

blastn -query /home/kimberly.ledger/rockfish_mb/data/dadasnake_output/filtered.seqs.fasta -db /home/kimberly.ledger/rockfish_mb/custom_db/rockfish_db_534_20250117 -out blastn_tissue_mock_bt22_534_20250117.txt -perc_identity 92 -qcov_hsp_perc 98 -num_threads 10 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames staxids' 

note: to make custom rockfish db,
- export fasta and metadata csv (Name and Organism) from Geneious
- example: makeblastdb -in rockfish_reference_db_534_20250117.fasta -dbtype nucl -out rockfish_db_534_20250117

### 2. R code 

*1_rockfish_taxonomic_assignment.Rmd:* This code takes the blastn output and determines the taxonomic assignment of each ASV   
*2_decontamination.Rmd:* This code removes ASVs without a rockfish taxonomic assignment (all samples) and discards low read depth replicates based on ASV accumulation curve (for eDNA samples only)   
*3A_format_mocks.Rmd:* This code formats expected mock community species compositions with the observed species proportions following Sebastes D-loop metabarcoding for quantitative metabarcoding model input   
*3B_qm_mock.Rm:* This code runs the quantitative metabarcoding model from Shelton et al. 2023 to generate estimates of amplification efficiency   
*3C_qm_mock_figures:* This code plots expected vs observed read proportions and amplification efficiency estimates from the mock communities (Figure 3)   
*4_tissue_test.Rmd:* This code evaluates the taxonomic assignment of known rockfish tissue samples using the Sebastes D-loop metabarcoding primers and plots Figure 2   
*5_aquaria_bt22.Rmd:* This code evaluates the eDNA samples from the ABL aquarium tank and bottom trawl survey (includes Figures 4 and 5)  
*6_tree_and_extra_stats.Rmd:* This code plots the phylogenetic tree (Figure 1A) and calculates some data summaries  
*7_haplotypes.Rmd:* This code computes and plots haplotypes for dusky, dark, and northern rockfishes (Figure 1B)  

### 3. Data 

data/blastn_tissue_mock_bt22_534_20250117.txt - blastn output  
data/rkfish_ref_dbs/rockfish_reference_db_534_20250117.csv - IDs of reference sequences   
data/metadata_tissue_mock_bt22.csv - sample metadata   
data/dadasnake_output/filtered.seqTab.RDS - ASV by sample table (dadasnake output)  
data/rockfish_mock.csv - expected species composition of mock communities    
data/ABLG_BURKE_rockfishDNA.csv - species IDs of rockfish tissue samples  
data/miniDloop_tissue.newick - tree used to order species in Figure 2   
data/miniDloop_AKspp_25supportthres.newick - tree of AK rockfish species used in Figure 1A   
data/cpv_fasta/cpv_minidloop_aligned.fasta - alignment of S. ciliatus, S. polyspinis, and S. variablis SebDLoop amplicon sequences   
data/cpv_fasta/cpv_meta.csv - metadata for cpv_minidloop_aligned.fasta  

note: due to the file size, rockfish_reference_db_534_20250117.fasta is included as a supplementary file
-----