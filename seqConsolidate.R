#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
suppressMessages(library("argparse"))
suppressMessages(library("pander"))
suppressMessages(library("stringr"))

code_dir <- dirname(
  sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=T)))

#' Set up and gather command line arguments
parser <- ArgumentParser(
  description = "R-based nucleotide sequence consolidater. Consolidate nucleotide sequences down to unique sequences and produce a key to revert back.")
parser$add_argument(
  "seqFile", nargs = 1, type = "character", default = NULL,
  help = "Sequence file to trim, either fasta or fastq format.")
parser$add_argument(
  "-o", "--output", nargs = 1, type = "character", default = NULL,
  help = "Output fasta file name. Ex. sample.consolidated.fasta")
parser$add_argument(
  "-k", "--keyFile", nargs = 1, type = "character", default = NULL,
  help = "Key file output name. Ex. sample.r1.csv")
parser$add_argument(
  "-l", "--seqName", nargs = 1, type = "character", default = NULL,
  help = "Name to append to unique sequences. Ex. sample.r1")
parser$add_argument(
  "--compress", action = "store_true", help = "Output fasta file is gzipped.") 

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

# Check I/O file types
seqType <- str_extract(args$seqFile, "fa[\\w]*")
if(!seqType %in% c("fa", "fasta", "fastq")){
  stop("Unrecognized sequence file type, please convert to '*.fasta' or '*.fastq'. Gzip compression is acceptable as well.")
}
seqType <- ifelse(seqType %in% c("fa", "fasta"), "fasta", "fastq")

if(!is.null(args$output)){
  outType <- str_extract(args$output, "fa[\\w]*")
  args$output <- unlist(strsplit(args$output, outType))[1]
  if(!outType %in% c("fa", "fasta", "fastq")){
    stop("Unrecognized output file type, please choose '*.fasta' or '*.fastq'.")
  }
  outType <- ifelse(outType %in% c("fa", "fasta"), "fasta", "fastq")
  if(outType == "fastq"){
    message("Since consolidation of sequences is only based on sequences, quality scores will be dropped. Output file will be in fasta format.")
    outType <- "fasta"
  }
  args$output <- paste0(args$output, outType)
}

if(!is.null(args$keyFile)){
  keyType <- str_extract(args$keyFile, "[\\w]+$")
  if(!keyType %in% c("csv", "tsv", "rds", "RData")){
    stop("Output key file type not supported. Please use csv, tsv, rds, or RData.")
  }
}
  
# Check sequence name lead
if(is.null(args$seqName)){
  args$seqName <- unlist(strsplit(args$seqFile, "fa[\\w]*"))[1]
}

# Print inputs to table
input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply(1:length(args), function(i){
    paste(args[[i]], collapse = ", ")}))
input_table <- input_table[
  match(c("seqFile :", "output :", "keyFile :", "seqName :"),
        input_table$Variables),]
pandoc.title("seqConsolidateR Inputs")
pandoc.table(data.frame(input_table, row.names = NULL), 
             justify = c("left", "left"), 
             split.tables = Inf)

# Load additional R-packages
addPacks <- c("ShortRead", "BiocGenerics", "Biostrings")
addPacksLoaded <- suppressMessages(
  sapply(addPacks, require, character.only = TRUE))
if(!all(addPacksLoaded)){
  pandoc.table(data.frame(
    "R-Packages" = names(addPacksLoaded), 
    "Loaded" = addPacksLoaded, 
    row.names = NULL))
  stop("Check dependancies.")
}

# Read sequence file
if(seqType == "fasta"){
  seqPointer <- ShortRead::readFasta(args$seqFile)
}else{
  seqPointer <- ShortRead::readFastq(args$seqFile)
}

seqs <- ShortRead::sread(seqPointer)
names(seqs) <- ShortRead::id(seqPointer)

factorSeqs <- factor(as.character(seqs))

key <- data.frame(
  "readNames" = names(factorSeqs),
  "seqID" = paste0(args$seqName, as.integer(factorSeqs))
)

consolidatedSeqs <- DNAStringSet(levels(factorSeqs))
names(consolidatedSeqs) <- paste0(args$seqName, 1:length(levels(factorSeqs)))

if(!is.null(args$output)){
  if(args$compress & !grepl(".gz", args$output)){
    args$output <- paste0(args$output, ".gz")
  }
}

# Write output and key files
# Output
if(is.null(args$output)){
  pandoc.table(
    data.frame(
      "seqID" = names(consolidatedSeqs),
      "sequence" = as.character(consolidatedSeqs),
      row.names = NULL),
    style = "simple")
}else{
  writeXStringSet(
    consolidatedSeqs, 
    filepath = args$output, 
    format = "fasta",
    compress = args$compress)
}

# Key file
if(!is.null(args$keyFile)){
  if(keyType == "csv"){
    write.csv(key, file = args$keyFile, row.names = FALSE, quote = FALSE)
  }else if(keyType == "tsv"){
    write.table(
      key, file = args$keyFile, sep = "\t", row.names = FALSE, quote = FALSE)
  }else if(keyType == "rds"){
    saveRDS(key, file = args$keyFile)
  }else if(keyType == "RData"){
    save(key, file = args$keyFile)
  }
}
q()
