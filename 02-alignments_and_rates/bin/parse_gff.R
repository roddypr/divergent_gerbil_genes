#!/usr/bin/R

# ---------------------------------------------------------------------------- #

library("optparse")

option_list <- list(
  make_option(c("-a", "--annotation_gff"), type="character", default=NULL,
              help="GFF file with transcript locations", metavar="character"),
  make_option(c("-o", "--out_csv"), type="character", default="out.csv",
              help="Output [default= %default]", metavar="character"))

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

# ---------------------------------------------------------------------------- #
# Input
## get the location of each transcript
tx      <- read.table(opt$annotation_gff, sep= "\t")

# ---------------------------------------------------------------------------- #
# Get trascript ID
tx_id   <- strsplit(as.character(tx$V9), split = ";")
tx_id   <- sapply(tx_id, function(x) x[grep("ID=", x)])
tx$tx   <- gsub("ID=transcript:", "", tx_id)

# Get gene name
tx_name <- strsplit(as.character(tx$V9), split = ";")
tx_name <- sapply(tx_name, function(x) x[grep("Name=", x)])
tx$name <- gsub("Name=", "", tx_name)

# Match trancript location with measurement table, based on transcript ID
colnames(tx)[c(1,4,5,7)] <- c('chr', 'start', 'end', 'strand')

tx <- tx[, c(1,4,5,7,10,11)]

# ---------------------------------------------------------------------------- #
# Write out
write.csv(tx, file = opt$out_csv, row.names = FALSE)

# ---------------------------------------------------------------------------- #
