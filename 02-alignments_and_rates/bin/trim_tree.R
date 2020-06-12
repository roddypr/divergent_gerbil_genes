#!/usr/bin/R

library("optparse")

option_list <- list(
  make_option(c("-i", "--index"), type="character", default=NULL,
              help="Table with a column for orthogroup ID and one column per species where '0' is absent and '1' is present", metavar="character"),
  make_option(c("-t", "--tree"), type="character", default=NULL,
              help="Tree to be trimmed (only topology matters)", metavar="character"),
  make_option(c("-o", "--out_dir"), type="character", default="trimmed_trees",
              help="Directory to output trimmed trees [default= %default]", metavar="character"))

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

# Inputs
orthogroup_index <- opt$index
species_tree     <- opt$tree
# Outpath
outdir         <- opt$out_dir

dir.create(outdir)

# ---------------------------------------------------------------------------- #
# Load orthogroup index
orthogroup_index <- read.table(orthogroup_index, header= TRUE)
stopifnot(colnames(orthogroup_index)[1] == "orthogroup")

# ---------------------------------------------------------------------------- #

# Get species names
orthogroup_index_ma <- orthogroup_index[, -1]
species             <- colnames(orthogroup_index)
orthogroup_index_ma <- as.matrix(orthogroup_index_ma)

# For each orthogroup, get which species are not present
stopifnot(orthogroup_index_ma %in% c(0,1))
orthogroup_index_na   <- t(apply(orthogroup_index_ma,
                                 1,
                                 function(x) x == 0))
species_remove        <- apply(orthogroup_index_na,
                               1,
                               function(x) names(which(x)))
names(species_remove) <- orthogroup_index$orthogroup

# ---------------------------------------------------------------------------- #

# Read species tree
library(ape)
tree           <- read.tree(species_tree)
tree$tip.label <- gsub(".pep", "", tree$tip.label)

# Make sure tree and matrix have same species
stopifnot(tree$tip.label %in% colnames(orthogroup_index_ma))
stopifnot(colnames(orthogroup_index_ma) %in% tree$tip.label)

# For each orthogroup, create a tree where we have removed the chosen species
for (i in 1:length(species_remove)) {

  orthogroup      <- names(species_remove)[i]
  orthogroup_tree <- tree
  if (length(species_remove[[i]]) > 0) {
    orthogroup_tree <- drop.tip(orthogroup_tree, species_remove[[i]])
  }
  # Write out
  out_path <- paste(outdir, "/", orthogroup, ".newick", sep = "")
  write.tree(orthogroup_tree, file = out_path)

}

# ---------------------------------------------------------------------------- #
# ENDS
# ---------------------------------------------------------------------------- #
