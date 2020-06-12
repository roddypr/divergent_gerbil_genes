#!/usr/bin/R

# ---------------------------------------------------------------------------- #

library("optparse")

option_list <- list(
  make_option(c("-a", "--annotation_csv"), type="character", default=NULL,
              help="csv ID locations", metavar="character"),
  make_option(c("-i", "--id"), type="character", default=NULL,
              help="Dictionary with gene IDs", metavar="character"),
  make_option(c("-o", "--out_csv"), type="character", default="out.csv",
              help="Output [default= %default]", metavar="character"))

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

# ---------------------------------------------------------------------------- #
# Input
## get the location of each transcript
tx  <- read.csv(opt$annotation_csv)

## Get the table with the ID of each gene
ids           <- read.table(opt$id, sep = " ")
colnames(ids) <- c("species", "gene", "tx", "protein", "name")

# ---------------------------------------------------------------------------- #
# Get trascript ID
matches    <- match(tx$tx, ids$tx)
tx$gene    <- ids$gene[matches]
tx$protein <- ids$protein[matches]

# ---------------------------------------------------------------------------- #
tx <- tx[, c(5, 1:4, 6:8)]
# ---------------------------------------------------------------------------- #
# Write out
write.csv(tx, file = opt$out_csv, row.names = FALSE)

# ---------------------------------------------------------------------------- #
