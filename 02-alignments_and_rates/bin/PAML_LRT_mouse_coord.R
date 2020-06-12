#!/usr/bin/R

# ---------------------------------------------------------------------------- #

library("optparse")

option_list <- list(
  make_option(c("-r", "--reference_csv"), type="character", default=NULL,
              help="csv ID locations", metavar="character"),
  make_option(c("-d", "--rates"), type="character", default=NULL,
              help="Table with rates measurements", metavar="character"),
  make_option(c("-o", "--out_csv"), type="character", default="out.csv",
              help="Output [default= %default]", metavar="character"))

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

# ---------------------------------------------------------------------------- #
# Input
## get the location of each transcript
mouse_reference <- read.csv(opt$reference_csv)

rates <- read.csv(opt$rates)
# rates <- read.table("input/muridae_dnds_per_orthogroup_wide.txt", header=TRUE)

# ---------------------------------------------------------------------------- #

mouse_reference <- mouse_reference[!is.na(mouse_reference$orthogroup), ]
mouse_reference <- mouse_reference[mouse_reference$orthogroup %in% rates$orthogroup,]
# ---------------------------------------------------------------------------- #

rates_per_mouse <- rates[match(mouse_reference$orthogroup, rates$orthogroup), c(2:ncol(rates))]
rates_per_mouse <- cbind(mouse_reference, rates_per_mouse)

# ---------------------------------------------------------------------------- #
# Write out
write.csv(rates_per_mouse, file = opt$out_csv, row.names = FALSE)

# ---------------------------------------------------------------------------- #
