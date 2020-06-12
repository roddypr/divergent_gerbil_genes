#!/usr/bin/R

# ---------------------------------------------------------------------------- #
# Parses trees with rates of sequence change measured by bppml+mapnh
# ---------------------------------------------------------------------------- #

library("optparse")

option_list <- list(
  make_option(c("-i", "--index"), type="character", default=NULL,
              help="Table with a column for orthogroup ID and one column per species where '0' is absent and '1' is present", metavar="character"),
  make_option(c("-d", "--in_dir"), type="character", default=NULL,
              help="Directory each orthogroup produced by testNh, each named according to 'orthogroup.dnd'", metavar="character"),
  make_option(c("-o", "--out_table"), type="character", default="out.txt",
              help="File name for output table [default= %default]", metavar="character"))

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

library(ape)

# Input directory
in_directory <- opt$in_dir
in_index     <- opt$index
# Output file
out_path     <- opt$out_table

# # Input directory
# in_directory <- "results/dnds"
# in_index     <- "results/successful_dnds_index"
# # Output file
# out_path <- "tmp/biopp/gerbil_measurements.txt"

# ---------------------------------------------------------------------------- #
# Fixed parameters: clade classification

# Focus clades
gerbils <- c("psammomys_obesus",
             "meriones_unguiculatus")
murinae <- c("mus_musculus",
             "rattus_norvegicus")

both_clades <- c(gerbils, murinae)


# ---------------------------------------------------------------------------- #
# Index
orthogroup_index <- read.table(in_index, header= TRUE)

#    Name of orthogroup
tree_paths <- list.files(in_directory)
tree_paths <- tree_paths[grep(".dnd", tree_paths)]

# Orthogroup names
orthogroups <- strsplit(tree_paths, split = "\\.")
orthogroups <- sapply(orthogroups, function(x) x[1])

# Check
stopifnot(orthogroup_index$orthogroup %in% orthogroups)
stopifnot(orthogroups %in% orthogroup_index$orthogroup)

# Get rate category (dn or ds) and mutation category (WW, WS, SW, SS)
measurement   <- strsplit(tree_paths, split = "\\.")
measurement   <- sapply(measurement, function(x) x[2])
measurement   <- gsub("counts", "", measurement)

measurement   <- strsplit(measurement, split = "_")
rate_category <- sapply(measurement, function(x) x[2])
mut_category  <- sapply(measurement, function(x) x[4])
mut_category  <- gsub("->", "", mut_category)

# Tree file
rates_matrix           <- matrix(NA,
                                 ncol = length(both_clades),
                                 nrow = length(tree_paths))
colnames(rates_matrix) <- both_clades

for (i in 1:length(tree_paths)) {
  # Read tree
  tree_path <- file.path(in_directory, tree_paths[i])
  tree      <- read.tree(file = tree_path)
  # tree$tip.label
  # str(tree)
  # plot(tree, edge.width = 2, label.offset = 0.2)
  # nodelabels()
  # tiplabels()
  # tree$tip.label

  # subset part of tree with
  gerbils_in_tree <- gerbils[gerbils %in% tree$tip.label]
  murinae_in_tree <- murinae[murinae %in% tree$tip.label]
  stopifnot(length(gerbils_in_tree) > 0)
  stopifnot(length(murinae_in_tree) > 0)

  # Common ancestor between gerbiles and murinae
  common_ancestor_node <- mrca(tree)[gerbils_in_tree[1], murinae_in_tree[1]]

  # Matrix: distance between each node
  dist_matrix   <- dist.nodes(tree)

  # Get the distance between the common ancestor node and each tip
  muridae        <- c(gerbils_in_tree, murinae_in_tree)
  muridae_nodes  <- match(muridae, tree$tip.label)
  distances      <- sapply(muridae_nodes,
        function(i) dist_matrix[common_ancestor_node, i])

  match_col <- match(muridae, colnames(rates_matrix))
  rates_matrix[i,match_col] <- distances
}

results_df <- data.frame(orthogroup        = orthogroups,
                         rate_category     = rate_category,
                         mutation_category = mut_category)

results_df <- cbind(results_df, as.data.frame(rates_matrix))

# ---------------------------------------------------------------------------- #

write.table(results_df,
            file = out_path,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# ---------------------------------------------------------------------------- #
# ENDS
# ---------------------------------------------------------------------------- #
