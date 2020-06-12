#!/usr/bin/R

library(tidyverse)
library("optparse")

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="File produced by muridae_tree_to_branch_length.R", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="out.txt",
              help="File name for wide output [default= %default]", metavar="character"))

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

in_path  <- opt$input
out_path <- opt$output
# in_path  <- "tmp/biopp/gerbil_measurements.txt"
# out_path <- "tmp/biopp/gerbil_measurements_wide.txt"

dn_ds <- read.table(in_path, header=TRUE)

# Gather the rates for all species into one column,
#    Then spread based on all the variables
dn_ds %>%
  gather("species", "value", -orthogroup, -rate_category, -mutation_category) %>%
  unite("id", c(rate_category, mutation_category, species)) %>%
  spread(id, value) -> dn_ds_wide

write.table(dn_ds_wide,
            file      = out_path,
            sep       = "\t",
            row.names = FALSE,
            quote     = FALSE)
