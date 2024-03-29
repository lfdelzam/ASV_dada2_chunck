# ASV_dada2_chunck
Rscript for Dada2 taxonomy affiliation. It reduces the memory requirements and is faster.
Primers and adapters must be removed before using this script with the following snakemake pipeline: https://github.com/biodiversitydata-se/amplicon-multi-cutadapt

# Installation 
 Create conda environment:
 
    conda create -n dada2_env -c bioconda -c conda-forge r-base=4.1.3 bioconductor-dada2=1.22.0 python=3.11.0
   
    conda activate dada2_env
   
    conda install -c conda-forge r-argparser=0.7.1

Clone repository:

    git clone https://github.com/lfdelzam/ASV_dada2_chunck

# Usage

Go to ASV_dada2_chunck directory and activate conda environment:

     cd ASV_dada2_chunck
     conda activate dada2_env

Run the R script:

     Rscript ASV_mergable_dada2_with_taxa_chunk.R [options]
   use -h flag to know all the options 

## When Reads using dada2 are not merged but concatenated, assign species to the input sequences by exact matching against a reference fasta using blastn

Create blast environment:

     conda create -n blast_env -c bioconda -c conda-forge blast=2.13.0 python=3.10.6 -y
     conda activate blast_env

Convert the dadas2 ASV_seqs.tsv output file into a fasta file with split fwd and rev sequences:

     python split_concatenated_reads.py [options]
   use -h flag to know all the options
   
For example:

      python split_concatenated_reads.py -i vibrio_from_gtdb_asv_sequences.tsv -o ASV_split_file

Create blastn database with your reference database, for instances, the custom reference detabase 16S_rRNA_from_BS_and_RefSeq_complete_genomes_db.fasta 

       mkdir -p blast_db
       makeblastdb -in 16S_rRNA_from_BS_and_RefSeq_complete_genomes_db.fasta -out blast_db/custom.Ntdb -dbtype nucl
   
Now, run blastn with 100% identity filter:

       blastn -query ASV_split_file -db blast_db/custom.Ntdb -outfmt 7 -out Custom_hits_16s -num_threads 8 -evalue 0.01 -strand 'both' -task 'blastn' -word_size 11 -max_target_seqs 1000 -perc_identity 100
   
Finally, obtain the hits with exact match of both fwd and rev reads to an unique specie:

      python parse_blast_assign_Specie.py [options]
   use -h flag to know all the options.

Here is an example:

       python parse_blast_assign_Specie.py -i Custom_hits_16s -a ASV_split_file -o custom_vibrio_from_gtdb -g 16S_rRNA_from_BS_and_RefSeq_complete_genomes_db.fasta -b custom

## Optional 

If you want to combine dada2 taxonomy (up to genus level) with Blastn specie affiliation, and/or create an ASV count-taxonomy table grouped by taxonomy affiliation.  

      python one_table.py [options] 
   use -h flag to display all the options 
