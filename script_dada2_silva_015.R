###23May2020
###I resequence my samples using macrogem service
###I modified this from script_dada2_silva_013.R because it is taking forever
###This script turn on multithread to expedite the process
##the initial script #refer to Daniel command script "script_dada2_033_full.R" to run my 16S R1R2 reads with dada2 using SILVA database version 132

# this will run with raw sequences
# revise the parameter as per Fede's suggestion - using default dada2 setting

# 24Jun - adding cutadapt command in this to remove primer and N from the sequences
# 25Jun - lost a lot after remove primer, adjust the parameter
# 30Jun - this is the finalized version with right directory to raw sequences and primer and N are properly removed from the raw sequences before merging


# cd /home/wenshu/reDPM3a2PT/R/
# R CMD BATCH script_dada2_silva_015.R script_dada2_silva_015.out

# =================================
# Swtiches to adjust the processing
# =================================

testing <- FALSE
do_cutadapt     <- TRUE
do_summary      <- TRUE
do_plot_quality <- TRUE
do_filtering <- TRUE  # If TRUE and primers are present must also do_cutadapt = TRUE
do_dada2    <- TRUE
do_taxo <- TRUE

#=================================
#   Parameters
#=================================

dataset_code <- "Tsunami_Wenshu_FedePrimer"
dataset_id <- 015
#dataset_path <- "/Users/wenshu/Dropbox/India-DPM3a2/dada2/16S_R1R2_002/"
dataset_path <- "/data/WENSHU/reDPM3aPT/raw/015/"
#these dada2-formatted training fasta files were derived from Silva project version 132 release.
#https://zenodo.org/record/1172783#.XbPofC-p3dc
silva_file <- "/home/wenshu/database/silva/silva_nr_v132_train_set.fa"
species_file <- "/home/wenshu/database/silva/silva_species_assignment_v132.fa"


# -- File structure

file_identifier = ".fastq"  # String to identify the files to be processed.

paired_reads = TRUE

R1_identifier = "_1"
R2_identifier = "_2"
file_name_separator = "_"
# This the first character of the file name to consider to extract the sample name (usually = 1)
sample.names_first_character = 1 


# --  Other parameters
# target
gene = "16S rRNA"
gene_region = "V4"
organelle = "nucleus"

# primer sets (926F-1392R)
FWD = "AAACTYAAAKGAATTGRCGG"
REV = "GGCGGTGTGTRC"
anchor = "^"  # Put in front of primer to anchor primer at start of sequence

# parameters for filterAndTrim
sequencer = "Illumina"

# parameters for filterAndTrim
truncLen = c(280, 240) # This influences the number of ASVs and the percent of asv recovered (choose based on the quality plot)
minLen = c(280,240)
truncQ = 2    # Usually use 2     
maxEE = c(2, 5) # Usually use 2-5  

if (paired_reads == FALSE) {
  truncLen = 250
  minLen = 250
  maxEE = maxEE[1]
}

# This is for 454 to remove long and bad reads
maxLen = 400  

# parameters for removeBimeraDenovo
method_chimera = "pooled"



#=================================
#   Load all required libraries
#=================================

library(dada2)
library(Biostrings)
library(ShortRead)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)



#=================================
#   Set up file directories
#=================================

path_dataset <- function(file_name) str_c(dataset_path, file_name)
file_dataset <- function(file_end)  str_c(dataset_path, "dada2/", file_end)

dir_fastq <-    path_dataset("fastq") # for raw sequences from illumina
dir_fastqN <-   path_dataset("fastqN/") # for filtered files removing reads with N
dir_cutadapt <- path_dataset("cutadapt/") # files after cutadapt
dir_filtered <- path_dataset("fastq_filtered/")
dir_qual <-     path_dataset("qual_pdf/")
dir_dada2 <-    path_dataset("dada2/")
dir_blast <-    path_dataset("blast/")



#=============================
#  Start running the script
#=============================

if(!dir.exists(dir_filtered)) dir.create(dir_filtered)
if(!dir.exists(dir_fastqN)) dir.create(dir_fastqN)
if(!dir.exists(dir_cutadapt)) dir.create(dir_cutadapt)
if(!dir.exists(dir_qual)) dir.create(dir_qual)
if(!dir.exists(dir_dada2)) dir.create(dir_dada2)
if(!dir.exists(dir_blast)) dir.create(dir_blast)

primer_length_fwd <- str_length(FWD)
primer_length_rev <- str_length(REV)



# =================================
#    Get the file names 
# =================================


fns <- sort(list.files(dir_fastq, full.names = TRUE))
fns <- fns[str_detect( basename(fns),file_identifier)]
if (testing) fns <- fns[1:12]

# fastq
if (paired_reads){
  
  fns_R1 <- fns[str_detect( basename(fns),R1_identifier)]
  fns_R2 <- fns[str_detect( basename(fns),R2_identifier)]
  
  fns_R1.fastq <- fns_R1
  fns_R2.fastq <- fns_R2
  
  # filters with reads with N removed
   fns_R1.filtN <- str_c(dir_fastqN, basename(fns_R1))  # Put N-filterd files in filtN/ subdirectory
   fns_R2.filtN <- str_c(dir_fastqN, basename(fns_R2))
  
  # after cutadapt
   fns_R1.cut <- str_c(dir_cutadapt, basename(fns_R1))  # Put files in /cutadapt subdirectory
   fns_R2.cut <- str_c(dir_cutadapt, basename(fns_R2))
  
  # after FilterAndTrim
  fns_R1.filt <- str_c(dir_filtered, basename(fns_R1))
  fns_R2.filt <- str_c(dir_filtered, basename(fns_R2))
  
  sample.names <- str_split(basename(fns_R1), pattern = file_name_separator, simplify = TRUE)
  
}  else {
  if (sequencer == "Illumina") {
    fns_R1 <- fns[str_detect( basename(fns),R1_identifier)]
    fns <- fns_R1
  }
  
  fns.fastq <- fns
  fns.filtN <- str_c(dir_fastqN, basename(fns))
  fns.cut <- str_c(dir_cutadapt, basename(fns))
  fns.filt <- str_c(dir_filtered, basename(fns))
  sample.names <- str_split(basename(fns), pattern = file_name_separator, simplify = TRUE)
}

sample.names <- sample.names[,1]
sample.names <- str_sub(sample.names, start=sample.names_first_character)
sample.names



# =================================
#   Get the number of reads in each file (just R1 if paired_reads)
# =================================

if (do_summary){
  
  print("=== Number of reads in each file ===")  
  
  summary <- data.frame()
  
  if (paired_reads){ 
    fns.summary <- fns_R1
  } else {
    fns.summary <- fns
  }
  
  for(i in 1:length(fns.summary)) {
    geom <- fastq.geometry(fns.summary[i])
    summary_one_row <- data.frame (n_seq=geom[1], file_name=basename(fns.summary[i]))
    summary <- bind_rows(summary, summary_one_row)
    print(paste("Finished with file", fns.summary[i], ". ", round(i/length(fns.summary)*100, 2), "%", sep=""))
  }
  
  write_tsv(summary, file_dataset( "_summary_after_cutadapt.tsv"))
}


if (do_cutadapt){

# =================================
#   Remove primers using cutadapt 
# =================================

# following: https://benjjneb.github.io/dada2/ITS_workflow.html

# ==========
# Step # 1   - Create all orientations of the input primers
# ==========

allOrients <- function(primer) {
  
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

# ==========
# Step # 2 - Remove any reads that contain N - can also be  done by cutadapt
# ==========

#if (paired_reads){
#  out_N <- filterAndTrim(fns_R1.fastq, fns_R1.filtN, fns_R2.fastq, fns_R2.filtN,
#                        maxN = 0, minQ = -10, multithread = TRUE)
#} else {
#  out_N <- filterAndTrim(fns.fastq, fns.filtN,
#                        maxN = 0, minQ = -10, multithread = TRUE)
#}
#out_N

# ==========
# Step # 3 - Check primers in one sample
# ==========

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

if (paired_reads){  
  rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fns_R1.fastq[[1]]), 
        FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fns_R2.fastq[[1]]), 
        REV.ForwardReads = sapply(REV.orients, primerHits, fn = fns_R1.fastq[[1]]), 
        REV.ReverseReads = sapply(REV.orients, primerHits, fn = fns_R2.fastq[[1]]))
} else {
  rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fns.fastq[[1]]),  
        REV.ForwardReads = sapply(REV.orients, primerHits, fn = fns.fastq[[1]]))
}

# =========
# Step #4 - Run cutadapt - to remove primer and N - don't work in R studio
# =========
 
# install cutadapt with python
#python3.6 -m pip install --user --upgrade cutadapt

        # cutadapt <- "/usr/bin/cutadapt"   # the cutadapt path on fede server
        # system2(cutadapt, args = "--version")   # Run sheel commands from R to see if cutadapt is working

 cutadapt <- "/usr/bin/cutadapt"
 
# CRITICAL - primers need to be in the right orientation!!

 FWD.RC <- dada2:::rc(FWD)
 REV.RC <- dada2:::rc(REV)

 if (paired_reads){
     # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
     R1.flags <- paste("-g", FWD, "-a", REV.RC)
     # R1.flags <- str_c("-g", anchor, FWD) # just FWD

     # Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
     R2.flags <- paste("-G", REV, "-A", FWD.RC)
     # R2.flags <- str_c("-G", anchor, REV)

     # Run cutadapt
    for (i in seq_along (fns_R1)) {
      system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,              # -n 2 required to remove FWD and REV from reads
                                 "--max-n=0",                              # remove N
                                 "-o", fns_R1.cut[i], "-p", fns_R2.cut[i], # output files
                                 fns_R1.fastq[i], fns_R2.fastq[i],         # input files
                                 "--discard-untrimmed",                    # remove all reads where primer not found
                                 "--minimum-length 50"))                   # remove reads that are too short (will cause error in Plot Quality)
 
# ================================
# Using now the cutadapt files as the starting files
# ===============================

  fns_R1 <- fns_R1.cut
  fns_R2 <- fns_R2.cut
  fns <- c(fns_R1, fns_R2)

  }
} else {
  # assume all reads to be in the good orientation
  # trim FWD and the reverse-complement of REV off of reads. The reverse cannot be anchored
  flags_FWD <- str_c("-g", anchor, FWD)
  flags_REV <- str_c("-a", REV.RC)

  ## Must run cutadapt in 2 steps to make sure both forward and reverse are removed and that all files
  # temp_file <- str_c(dir_cutadapt, "temp.fastq.gz")

  # Run cutadapt
  for (i in seq_along(fns)) {
    system2(cutadapt, args = c(flags_FWD,
                               flags_REV,
                               "-o", fns.cut[i],        #output files
                               fns.fastq[i],            #input files
                               "--max-n=0",             #remove any N
                               "--discard-untrimmed"))  #remove all reads where primer not found

    system2(cutadapt, args = c(flags_REV,
                               "-o", fns.cut[i],        #output files
                               temp_file))              #keep reads even if reverse primer not found

   }


# ================================
# Using cutadapt files as the starting files
# ================================
   fns <- fns.cut
 }
}


# =================================
#   Get the number of reads in each file (just R1)
# =================================

if (do_summary){
  
  print("=== Number of reads in each file ===")  
  
  summary <- data.frame()
  
  if (paired_reads){ 
    fns.summary <- fns_R1
  } else {
    fns.summary <- fns
  }
  
  for(i in 1:length(fns.summary)) {
    geom <- fastq.geometry(fns.summary[i])
    summary_one_row <- data.frame (n_seq=geom[1], file_name=basename(fns.summary[i]))
    summary <- bind_rows(summary, summary_one_row)
    print(paste("Finished with file", fns.summary[i], ". ", round(i/length(fns.summary)*100, 2), "%", sep=""))
  }
  
  write_tsv(summary, file_dataset( "_summary_after_cutadapt.tsv"))
}

# =================================
#   Plot quality after cutadapt
# =================================
if (do_plot_quality) {
  
  print("=== Plot quality ===")   
  
  for(i in 1:length(fns)) {
    print(str_c("i = ", i))
    p1 <- plotQualityProfile(fns[i])
    # if (i <= 2) {print(p1)}
    p1_file <- paste0(dir_qual, basename(fns[i]),".qual.pdf")
    ggsave( plot=p1, filename= p1_file, device = "pdf", width = 15, height = 15, scale=1, units="cm")
    
    read_length <- data.frame(length = width(ShortRead::readFastq(fns[i]))) # Read the fastQ file and get the length of all the reads...
    print(str_c("File before filtration", fns[i], "- Read length min=", min(read_length$length),"max=", max(read_length$length), "mean=", mean(read_length$length, na.rm=TRUE),  sep=" "))
    
    print(paste("Finished with file", fns[i], ". ", sep=""))
  }
  
}

# =================================
#   Filtering
# =================================
if (do_filtering) {
  
  print("=== Filtering ===") 
  
  if (paired_reads){
    out <- filterAndTrim(fns_R1, fns_R1.filt,fns_R2, fns_R2.filt, 
                         maxN=0, rm.phix=TRUE,
                         truncLen=truncLen, minLen=minLen, maxEE=maxEE,truncQ=truncQ,
                         compress=TRUE, multithread = TRUE)
    
    fns.filt <- c(fns_R1.filt, fns_R2.filt)
    
  } else { 
    out <- filterAndTrim(fns, fns.filt, 
                         maxN=0, rm.phix=TRUE,
                         truncLen=truncLen[1], maxLen = maxLen, minLen=minLen[1], maxEE=maxEE[1],truncQ=truncQ,
                         compress=TRUE, multithread = FALSE)
  }  
  
  print(out)
  write_tsv(data.frame(out), file_dataset( "_summary_filtered_files.tsv")) 
}

# =================================
#   Dada2 block
# =================================  

if (do_dada2) {  
  
  print("=== Start Dada2 ===") 
  
  # =================================
  #   Error
  # =================================
  
  print("=== Error computing ===") 
  
  if (paired_reads){ 
    err_R1 <- learnErrors(fns_R1.filt, multithread = TRUE)
    p <- plotErrors(err_R1, nominalQ=TRUE)
    p_file <- file_dataset("_LearnErrors_R1.pdf")
    ggsave( plot=p, filename= p_file, device = "pdf", 
            width = 15, height = 15, scale=1, units="cm")
    
    
    
    err_R2 <- learnErrors(fns_R2.filt, multithread = TRUE)
    p <- plotErrors(err_R2, nominalQ=TRUE)
    p_file <- file_dataset("_LearnErrors_R2.pdf")
    ggsave( plot=p, filename= p_file, device = "pdf", 
            width = 15, height = 15, scale=1, units="cm")
  } else {
    err_R1 <- learnErrors(fns.filt, multithread = FALSE)
    p <- plotErrors(err_R1, nominalQ=TRUE)
    p_file <- file_dataset("_LearnErrors_R1.pdf")
    ggsave( plot=p, filename= p_file, device = "pdf", 
            width = 15, height = 15, scale=1, units="cm")
  }
  
  
  
  # =================================
  #   Dereplicate
  # =================================
  
  print("=== dereplicate ===") 
  
  if (paired_reads){ 
    derep_R1 <- derepFastq(fns_R1.filt, verbose=T)
    derep_R2 <- derepFastq(fns_R2.filt, verbose=T)
    
    names(derep_R1) <- sample.names
    names(derep_R2) <- sample.names
  } else {
    derep_R1 <- derepFastq(fns.filt, verbose=T)
    
    names(derep_R1) <- sample.names
  }
  
  
  # =================================
  #   dada
  # ================================= 
  
  print("=== dada ===")  
  
  
  if (paired_reads){   
    dada_R1 <- dada(derep_R1, err=err_R1, multithread = TRUE)
    dada_R2 <- dada(derep_R2, err=err_R2, multithread = TRUE)
    
    # dada_R1[[1]]
    # dada_R2[[1]]
  }  else {
    # following https://benjjneb.github.io/dada2/faq.html#can-i-use-dada2-with-my-454-or-ion-torrent-data
    if (sequencer == "454") {
      dada_R1 <- dada(derep_R1, err=err_R1, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32, multithread = FALSE)
    } else{
      dada_R1 <- dada(derep_R1, err=err_R1, multithread = FALSE)
    }
    seqtab <- makeSequenceTable(dada_R1)
  }
  # =================================
  #  Merge pairs
  # =================================  
  
    print("=== Merging pairs ===")  
  
    if (paired_reads){ 
      mergers <- mergePairs(dada_R1, derep_R1, dada_R2, derep_R2, verbose=TRUE)
  
      seqtab <- makeSequenceTable(mergers)
    } 
  
    t_seqtab <- t(seqtab)
  
    print(table(nchar(getSequences(seqtab))))
  
    print(sprintf("Mean asv length : %.2f", mean(nchar(getSequences(seqtab)))))
  
  
  # =================================
  #  Remove chimera
  # =================================  
  
  print("=== Remove Chimera ===")  
  
  seqtab.nochim <- removeBimeraDenovo(seqtab, method=method_chimera, multithread=TRUE, verbose=TRUE)
  
  p <- ggplot(data.frame(seq_length=nchar(getSequences(seqtab.nochim)))) +
    geom_histogram(aes(x=seq_length)) +
    ggtitle(str_c("Number of asv: ", ncol(seqtab.nochim)))
  p_file <- file_dataset("_asv_length_hist.pdf")
  ggsave( plot=p, filename= p_file, device = "pdf", 
          width = 15, height = 15, scale=1, units="cm")
  
  # =================================
  #  saveRDS
  # ================================= 
  
  saveRDS(seqtab.nochim, file_dataset("_seqtab.nochim.rds")) 
  
  # =================================
  #  Compile number of reads at each step
  # ================================= 
  
  print("=== Compile number of reads at each step ===")  
  
  getN <- function(x) sum(getUniques(x))
  
  if(paired_reads){
    track <- cbind(sapply(dada_R1, getN), 
                   sapply(mergers, getN), 
                   rowSums(seqtab), 
                   rowSums(seqtab.nochim))
    colnames(track) <- c("denoised", "merged", "tabled", "nonchim")
  } else {
    track <- cbind(sapply(dada_R1, getN), 
                   rowSums(seqtab), 
                   rowSums(seqtab.nochim))    
    colnames(track) <- c("denoised", "tabled", "nonchim")
  }
  
  
  track <- data.frame(track) %>% 
    mutate(file_code = sample.names)
  
  write_tsv(track, file_dataset("_summary_dada2.txt"))
  
  # =================================
  #    Write fasta file without taxo
  # =================================   
  
  seqtab.nochim_trans <- as.data.frame(t(seqtab.nochim)) %>% 
    rownames_to_column(var = "sequence") %>%
    rowid_to_column(var = "asv_number") %>%
    mutate(asv_code = sprintf("asv_%03d_%05d", dataset_id, asv_number)) %>% 
    mutate(sequence = str_replace_all(sequence, "(-|\\.)",""))
  
  seq_out <- Biostrings::DNAStringSet(seqtab.nochim_trans$sequence)
  names(seq_out) <- seqtab.nochim_trans$asv_code
  Biostrings::writeXStringSet(seq_out, file_dataset("_no_taxo.fasta"), 
                              compress=FALSE, width = 20000)
}



if (do_taxo) {
  
  print("=== Assigning Taxonomy ===")
  
  # =================================
  #  Reload the files
  # =================================
  
  seqtab.nochim <- readRDS(file_dataset("_seqtab.nochim.rds"))
  
  # =================================
  #   Assign taxonomy
  # =================================  
  
  taxa <- assignTaxonomy(seqs=seqtab.nochim,
                         refFasta=silva_file,
                         minBoot = 0, outputBootstraps = TRUE,
                         multithread = TRUE, verbose = TRUE)
 
  #taxa <- addSpecies(taxa, refFasta = species_file)  
  saveRDS(taxa, file_dataset("_taxa.rds")) 
  
  
  # =================================
  #  For debugging reload the files
  # =================================
  
  #  taxa <- readRDS(file_dataset("_taxa.rds"))
  
  # =================================
  #    Create the ASV table
  # =================================  
  
  boot <- data.frame(taxa$boot) %>% 
    rename_all(dplyr::funs(str_c(.,"_boot"))) %>% 
    rownames_to_column(var = "sequence")
  taxa.df <- data.frame(taxa$tax)  %>% 
    rownames_to_column(var = "sequence")
  
  seqtab.nochim_trans <- as.data.frame(t(seqtab.nochim)) %>% 
    rownames_to_column(var = "sequence") %>%
    rowid_to_column(var = "asv_number") %>%
    mutate(asv_code = sprintf("asv_%03d_%05d", dataset_id, asv_number)) %>% 
    mutate(sequence = str_replace_all(sequence, "(-|\\.)","")) %>% 
    left_join(taxa.df) %>% 
    left_join(boot) 
  
  write_tsv(seqtab.nochim_trans, file_dataset("_dada2.tsv"), na="")
  
  # =================================
  #    Create tables for import into database
  # =================================  
  
  metapr2_asv <- seqtab.nochim_trans %>% 
    select(asv_code,sequence, asv_code:dataset_id)
  
  metapr2_asv$sequence_hash = purrr::map_chr(metapr2_asv$sequence,digest::sha1)
  
  write_tsv(metapr2_asv, file_dataset("_metapr2_asv.txt"), na="")
  
  metapr2_asv_abundance <- seqtab.nochim_trans %>% 
    select(-asv_number, -sequence, -(Kingdom:Genus_boot)) %>% 
    gather("file_code", "n_reads", -contains("asv_code")) %>% 
    filter(n_reads > 0 )
  
  write_tsv(metapr2_asv_abundance, file_dataset("_metapr2_asv_abundance.txt"), na="")
  
  # =================================
  #    Write fasta file with taxo
  # =================================   
  
  seq_out <- Biostrings::DNAStringSet(seqtab.nochim_trans$sequence)
  names(seq_out) <- str_c(seqtab.nochim_trans$asv_code,seqtab.nochim_trans$species, sep="|")
  Biostrings::writeXStringSet(seq_out, file_dataset("_taxo.fasta"), 
                              compress=FALSE, width = 20000)
}

# =================================
#    Clean up
# =================================  

rm(list = ls())
