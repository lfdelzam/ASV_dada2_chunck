# ASV_dada2_chunck
Rscript for Dada2 taxonomy affiliation. It reduces the memory requirements and it is faster.
Primers and adapters must be removed before using this script (please use the snakemake pipeline: https://github.com/biodiversitydata-se/amplicon-multi-cutadapt)

# Installation 

    conda create -n dada2_env -c bioconda -c conda-forge r-base=4.1.3 bioconductor-dada2=1.22.0 python=3.11.0
   
    conda activate dada2_env
   
    conda install -c conda-forge r-argparser=0.7.1

# Usage

     Rscript ASV_mergable_dada2_with_taxa_chunk.R [options]
   use -h flag to know all the options 

     python one_table.py [options]
   use -h flag to know all the options 