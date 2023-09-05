rm(list = ls())
options(warn=-1)
suppressPackageStartupMessages(library(argparser))

# arguments
p <- arg_parser("ASV_mergable_dada2_with_taxa_chunk.R - use it after removing adapters and primers using cutadapt")
p <- add_argument(p, "-d", help="directory containing the fastq files, *.gz. This folder should contain trimmed reads. It may contain none or several folders (from each illumina lane used that will be merged into one final table)", default="Cutadapt_results")
p <- add_argument(p, "-f", help="forward reads extention", default="_L001_R1_001.fastq.fq.gz")
p <- add_argument(p, "-r", help="reverse reads extention", default="_L001_R2_001.fastq.fq.gz")
p <- add_argument(p, "-l", help="trunlenleft, according to FastQC quality figures - 3' end of trimmed R1 reads", default=265) 
p <- add_argument(p, "-t", help="trunlenright, according to FastQC quality figures - 3' end of trimmed R2 reads", default=240) 
p <- add_argument(p, "-s", help="mergeSequenceTables -repeats, when merging Illumina runs", default="sum")
p <- add_argument(p, "-c", help="assingTaxonomy - dada2 parameter: tryRC", default=TRUE)
p <- add_argument(p, "-m", help="allowMultiple - dada2 parameter", default=TRUE)
p <- add_argument(p, "-g", help="path to taxonomy database(s), containing compressed files , *.gz",default="DBs_reference_amplicon/gtdb_16s_db") #DBs_reference_amplicon/SILVA138_db,
p <- add_argument(p, "-b", help="reference database(s), options: gtdb, silva, pr2. Pr2 version >= 5.0", default="gtdb")
p <- add_argument(p, "-u", help="Chunck size when assingTaxonomy", default=4000)
p <- add_argument(p, "-z", help="JustConcatenate -dada2 mergePairs, TRUE or FALSE. When using Novaseq, set it to TRUE", default=FALSE)
p <- add_argument(p, "-O", help="Outdir", default="Dada2_Results")

argv <- parse_args(p)

species_add_step=TRUE
#set species_add_step to FALSE if reads will be concatenated (argv$z=TRUE)
if (argv$z) species_add_step=FALSE 

set.seed(123)
dada2_ASV_table<-function(path_run){
	
  FR_outfolder=paste(argv$O, paste("Run",basename(path_run), sep="_"), sep="/")
  fnFs <- sort(list.files(path=path_run, pattern=argv$f, full.names = TRUE))
  fnRs <- sort(list.files(path=path_run, pattern=argv$r, full.names = TRUE))
  # Extract sample names
  sample.names <- gsub(argv$f, "", basename(fnFs))
  #plotQualityProfile(fnFs[1:2])
  # Place filtered files in filtered/ subdirectory

  filtFs <- file.path(FR_outfolder,  "filtered_reads", paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(FR_outfolder,  "filtered_reads", paste0(sample.names, "_R_filt.fastq.gz"))
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names

  out = filterAndTrim(
    fnFs, filtFs, fnRs, filtRs,truncLen=c(argv$l,argv$t),
    trimLeft=c(0,0),
    maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
    compress=TRUE, multithread=TRUE
  )

  ###
  not.lost <- file.exists(filtFs)
  filtFs <- filtFs[not.lost]
  filtRs <- filtRs[not.lost]
  if (length(sample.names[!not.lost]) > 0) { cat("The following samples did not pass the filter:\n",sample.names[!not.lost],"\n" ) }
  sample.names <- sample.names[not.lost]
  out <- out[not.lost,]


  ###

  errF <- learnErrors(filtFs, multithread=TRUE)
  pdf(paste(FR_outfolder,paste(argv$l,"Dada2_errF.pdf",sep="_"), sep="/"))
  plotErrors(errF, nominalQ=TRUE)
  dev.off()

  errR <- learnErrors(filtRs, multithread=TRUE)
  pdf(paste(FR_outfolder,paste(argv$t,"Dada2_errR.pdf",sep="_"), sep="/"))
  plotErrors(errR, nominalQ=TRUE)
  dev.off()

  dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
  dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
  mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, justConcatenate=argv$z)
  seqtab <- makeSequenceTable(mergers)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

  getN <- function(x) sum(getUniques(x))
  track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track) <- sample.names

  write.table(track, file =paste(FR_outfolder,"Dada2_track.tsv", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = TRUE,
              col.names = TRUE, qmethod = c("escape", "double"),
              fileEncoding = "")

  return(seqtab.nochim)
}

#libraries
list_of_packages <- c("dada2")
for (s in list_of_packages) { suppressPackageStartupMessages(library(s, character.only = TRUE))}

dir.create(argv$O, showWarnings = FALSE)


if (!file.exists(paste(argv$O,"mergetab_R",sep="/"))) {

  runs_list=list.dirs(argv$d, recursive = FALSE)

  if (length(runs_list) > 1){
	input_tables=list()
	counter=1
	for (folder in runs_list) {
 		cat("Analysing run", folder, "\n")
  		input_tables[[counter]]= dada2_ASV_table(folder)
  		counter = counter + 1
	        }

  #Merging tables

  mergetab <- mergeSequenceTables(tables=input_tables, repeats = argv$s)
  #repeats (Optional). Default "error". Specifies how merging should proceed in the presence of repeated sample names. Valid values: "error", "sum". If "sum", then samples with the same name are summed together in the merged table.

  } else {
	cat("Analysing data in ", argv$d, "\n")
        mergetab=dada2_ASV_table(argv$d)
  }

saveRDS(mergetab, file=paste(argv$O,"mergetab_R",sep="/"))

} else {

  mergetab=readRDS(file=paste(argv$O,"mergetab_R",sep="/"))
}

#COnverting seqs to ASV
  if (file.exists(paste(argv$O,"ASV_seqs.tsv", sep="/"))) {
  ASV_df=read.delim2(paste(argv$O,"ASV_seqs.tsv", sep="/") , header = TRUE, sep = "\t", quote = "",
              dec = ",", fill = TRUE, comment.char = "")
   } else {
    seqtab_copy = data.frame(ASV=paste("ASV", 1:ncol(mergetab), sep = "") , t(mergetab))
    ASV_df=data.frame(ASV=seqtab_copy$ASV,Sequence=colnames(mergetab))

    cat("Printing out count tables\n")
    write.table(mergetab, file = paste(argv$O,"ASVseqs_count_table.tsv", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

    write.table(seqtab_copy, file =paste(argv$O,"ASV_counts.tsv", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

    cat("Printing out ASV sequences\n")
    write.table(ASV_df, file =paste(argv$O,"ASV_seqs.tsv", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
    }


#taxonomy
DBs=strsplit(argv$b, ",")[[1]]
Path_dbs=strsplit(argv$g, ",")[[1]]

for (i in 1:length(DBs)) {

  DB=DBs[i]
  cat("assignTaxonomy with database", DB, "\n")
  Path_db=Path_dbs[i]

  if (DB == "gtdb") {
    file1=paste(Path_db, "gtdb-sbdi-sativa.r06rs202.assignTaxonomy.fna.gz", sep="/")
    file2=paste(Path_db, "gtdb-sbdi-sativa.r06rs202.addSpecies.fna.gz", sep="/")
    taxLEVELS= c("Kingdom","Domain", "Phylum", "Class", "Order", "Family", "Genus")
  }

  if (DB == "silva") {

    file1=paste(Path_db, "silva_nr99_v138.1_train_set.fa.gz", sep="/")
    file2=paste(Path_db, "silva_species_assignment_v138.1.fa.gz", sep="/")
    taxLEVELS= c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
  }

  if (DB == "pr2") {

    file1=paste(Path_db, "pr2_version_5.0.0_SSU_dada2.fasta.gz", sep="/")
    taxLEVELS= c("Domain","Supergroup","Division","Subdivision","Class","Order","Family","Genus","Species") # new tax leves with pr2v5
    species_add_step=FALSE
  }

if (!file.exists(paste(argv$O,"tax_R",sep="/"))) {

      cat("assignTaxonomy in chunks\n")
      Seqs=colnames(mergetab)
      eND=length(Seqs)
      sTEP=argv$u

      if (species_add_step) {
        COLNAMES=c(taxLEVELS, "Species", paste("Boot",taxLEVELS, sep="_"))
      } else {
      COLNAMES=c(taxLEVELS, paste("Boot",taxLEVELS, sep="_"))
      }

      taxadf=data.frame(matrix(ncol=length(COLNAMES), nrow=0))

      for (i in seq(1, eND, sTEP ) ) {

          chunk=i+sTEP-1

          if (chunk >= eND) { chunk=eND }
          cat("   assigning taxonomy seqs ", i, " to ", chunk, "(",round(chunk*100/eND,2), "% )\n")

          taxasub=assignTaxonomy(Seqs[i:chunk], file1, multithread=TRUE, taxLevels = taxLEVELS, tryRC = argv$c, outputBootstraps=TRUE)
          if (species_add_step) {
            cat("   addSpecies\n")
            taxasub$tax <- addSpecies(taxasub$tax, file2, allowMultiple = argv$m, tryRC = argv$c, n=argv$u)
            #allowMultiple (Optional). Default FALSE. Defines the behavior when exact matches to multiple (different) species are found. By default, only unambiguous identifications are returned. If set to TRUE, a concatenated string of all exactly matched species is returned.
          }

        	tdftax=as.data.frame(taxasub$tax)
        	tdfboot=as.data.frame(taxasub$boot)
        	tdf= merge(tdftax, tdfboot,by = 'row.names', all = TRUE)
        	row.names(tdf)=tdf$Row.names
        	tdf=tdf[,-1]
        	taxadf=rbind(taxadf,tdf)
      }


      names(taxadf)=COLNAMES
      taxadf$Sequence=row.names(taxadf)

  taxa=merge(ASV_df, taxadf, by="Sequence" )
  taxa=taxa[order(as.integer(gsub("ASV","",taxa$ASV))),]

  saveRDS(taxa, file=paste(argv$O,"tax_R",sep="/"))
} else {

taxa=readRDS(file=paste(argv$O,"tax_R",sep="/"))

}

  cat("printing out taxonomy tables\n")
  write.table(taxa, file = paste(argv$O,paste(DB,"NEW_ASVseqs_taxa_table.tsv", sep="_"), sep="/"), append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = c("Sequence","ASV",COLNAMES), qmethod = c("escape", "double"),
              fileEncoding = "")


  write.table(taxa[!names(taxa) %in% "Sequence" ], file =paste(argv$O,paste(DB,"NEW_ASV_taxonomy.tsv",sep="_"), sep="/"), append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = c("ASV", COLNAMES), qmethod = c("escape", "double"),
              fileEncoding = "")


}
