#!/usr/bin/R

# ---------------------------------------------------------------------------- #

library("optparse")

option_list <- list(
  make_option(c("-a", "--annotation_csv"), type="character", default=NULL,
              help="csv ID locations", metavar="character"),
  make_option(c("-m", "--orthogroup_ma"), type="character", default=NULL,
              help="Orthogroup list", metavar="character"),
  make_option(c("-s", "--species"), type="character", default="mus_musculus",
              help="Orthogroup list", metavar="character"),
  make_option(c("-o", "--out_csv"), type="character", default="out.csv",
              help="Output [default= %default]", metavar="character"))

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

# ---------------------------------------------------------------------------- #
# Input
## get the location of each transcript
tx  <- read.csv(opt$annotation_csv)

# Which species?
species <- opt$species

# Load orthogroups
ortho           <- read.table(opt$orthogroup_ma,
                              header=TRUE,
                              sep="\t", row.names = 1, stringsAsFactors=FALSE)

# Rename species
colnames(ortho) <- gsub(".pep", "", colnames(ortho))
stopifnot(species %in% colnames(ortho))
# Separate genes for each species in each orthgroup
# Count the number of genes for each species for each orthogroup
ortho <- as.matrix(ortho)
ortho_mouse <- ortho[,species]
ortho_mouse <- lapply(ortho_mouse, function(x) unlist(strsplit(x, split = ", ")))

ortho_mouse_names  <- lapply(1:length(ortho_mouse),
                             function(i) rep(names(ortho_mouse)[i],
                                             length(ortho_mouse[[i]])))
ortho_mouse_names  <- unlist(ortho_mouse_names)
ortho_mouse_vector <- unlist(ortho_mouse)

# ---------------------------------------------------------------------------- #

tx$orthogroup <- ortho_mouse_names[match(tx$protein, ortho_mouse_vector)]

# ---------------------------------------------------------------------------- #
# Write out
write.csv(tx, file = opt$out_csv, row.names = FALSE)

# ---------------------------------------------------------------------------- #
