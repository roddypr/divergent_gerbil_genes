#!/usr/bin/R

# Pseudo-code:

# For each gene:
#   1 - Load alignment to find out the number ID of the muridae species
#   2 - Load tree to find the IDs of the internal nodes
#   3 - Load "dN & dS for each branch" table, parsed from null model paml output
#   4 - Get the values for the specific nodes and make a very long line of results

correct_short_names <- function(species_names) {
  correct_names <- c("mus_musculus",
                     "rattus_norvegicus",
                     "meriones_unguiculatus",
                     "psammomys_obesus")
  short_version <- c("M_musculus",
                     "R_norvegicus",
                     "M_unguiculatus",
                     "P_obesus")

  correction <- correct_names[match(species_names, short_version)]
  if (any(!is.na(correction))) {
    to_correct <- which(!is.na(correction))
    species_names[to_correct] <- correction[to_correct]
  }
  return(species_names)
}

# Functions to parse in each type of data
species_name_from_alignment <- function(gene, alignment_path) {
  my_path       <- paste(gene, "fa", sep = ".")
  my_path       <- file.path(alignment_path, my_path)
  dna_alignment <- readDNAStringSet(my_path)
  species_names <- names(dna_alignment)
  # Correct species name (e.g. M_musculus -> mus_musculus)
  species_names <- correct_short_names(species_names)
  return(species_names)
}

read_tree <- function(gene, tree_path) {
  my_path       <- paste(gene, "newick", sep = ".")
  my_path       <- file.path(tree_path, my_path)

  tree          <- read.tree(my_path)
  tree$tip.label <- gsub("\\#.+", "", tree$tip.label)

  tree$tip.label <- correct_short_names(tree$tip.label)
  return(tree)
}

read_branch_lengths <- function(gene, branch_length_paths) {
  my_path <- file.path(branch_length_paths, gene)
  df      <- read.table(my_path, header = FALSE)
  colnames(df) <- c("branch", "t", "N", "S", "dN/dS", "dN", "dS", "N*dN", "S*dS")

  branch  <- strsplit(as.character(df$branch), split = "\\.\\.")
  branch  <- lapply(branch, as.numeric)
  branch  <- lapply(branch, sort)
  df$from <- sapply(branch, function(x) x[1])
  df$to   <- sapply(branch, function(x) x[2])

  return(df)
}

df2matrix <- function(df, column_name) {

  node_number <- max(c(df$from, df$to))
  results_ma <- matrix(NA, nrow = node_number, ncol = node_number)

  colnames(results_ma) <- 1:ncol(results_ma)
  rownames(results_ma) <- 1:nrow(results_ma)

  results_ma[cbind(df$from, df$to)] <- df[, column_name]

  return(results_ma)

}

rename_colnames <- function(meas_matrix, muridae_id) {
    stopifnot(colnames(meas_matrix) == rownames(meas_matrix))
    muridae_id <- muridae_id[!is.na(muridae_id)]
    meas_col   <- colnames(meas_matrix)
    meas_col[match(muridae_id, meas_col)] <- names(muridae_id)
    colnames(meas_matrix) <- meas_col
    rownames(meas_matrix) <- meas_col
    return(meas_matrix)
}

parse_distances <- function(meas_matrix, tree, gerbils, murinae) {

  gerbils_tree <- gerbils[gerbils %in% tree$tip.label]
  murinae_tree <- murinae[murinae %in% tree$tip.label]

  # MRCA node between all murids
  murid_mrca <- mrca(tree)[gerbils_tree[1], murinae_tree[1]]

  if (length(gerbils_tree) == 2) {
    # MRCA node between gerbils
    gerbil_mrca       <- mrca(tree)[gerbils_tree[1], gerbils_tree[2]]
    gerbilmrca_to_pob <- meas_matrix["psammomys_obesus", gerbil_mrca]
    gerbilmrca_to_mun <- meas_matrix["meriones_unguiculatus", gerbil_mrca]

    # MRCA node between gerbils --- MRCA node between all murids
    muridmrca_to_gerbilmrca <- meas_matrix[murid_mrca, gerbil_mrca]

    # MRCA node to each gerbil
    muridmrca_to_pob <- gerbilmrca_to_pob + muridmrca_to_gerbilmrca
    muridmrca_to_mun <- gerbilmrca_to_mun + muridmrca_to_gerbilmrca


  } else if ("psammomys_obesus" %in% tree$tip.label) {

    # MRCA node between gerbils
    gerbilmrca_to_pob <- NA
    gerbilmrca_to_mun <- NA

    # MRCA node between gerbils --- MRCA node between all murids
    muridmrca_to_gerbilmrca <- NA

    # MRCA node to each gerbil
    muridmrca_to_pob <- meas_matrix["psammomys_obesus", murid_mrca]
    muridmrca_to_mun <- NA


  } else if ("meriones_unguiculatus" %in% tree$tip.label) {

    # MRCA node between gerbils
    gerbilmrca_to_pob <- NA
    gerbilmrca_to_mun <- NA

    # MRCA node between gerbils --- MRCA node between all murids
    muridmrca_to_gerbilmrca <- NA

    # MRCA node to each gerbil
    muridmrca_to_pob <- NA
    muridmrca_to_mun <- meas_matrix["meriones_unguiculatus", murid_mrca]

  } else {
    stop()
  }

  if (length(murinae_tree) == 2) {
    # MRCA node between murins
    murine_mrca       <- mrca(tree)[murinae_tree[1], murinae_tree[2]]
    murinemrca_to_mmu <- meas_matrix["mus_musculus", murine_mrca]
    murinemrca_to_rno <- meas_matrix["rattus_norvegicus", murine_mrca]

    # MRCA node between murines --- MRCA node between all murids
    muridmrca_to_murinemrca <- meas_matrix[murid_mrca, murine_mrca]

    # MRCA node to each murine
    muridmrca_to_mmu <- murinemrca_to_mmu + muridmrca_to_murinemrca
    muridmrca_to_rno <- murinemrca_to_rno + muridmrca_to_murinemrca


  } else if ("mus_musculus" %in% tree$tip.label) {

    # MRCA node between murines
    murinemrca_to_mmu <- NA
    murinemrca_to_rno <- NA

    # MRCA node between murines --- MRCA node between all murids
    muridmrca_to_murinemrca <- NA

    # MRCA node to each murine
    muridmrca_to_mmu <- meas_matrix["mus_musculus", murid_mrca]
    muridmrca_to_rno <- NA


  } else if ("rattus_norvegicus" %in% tree$tip.label) {

    # MRCA node between murines
    murinemrca_to_mmu <- NA
    murinemrca_to_rno <- NA

    # MRCA node between murines --- MRCA node between all murids
    muridmrca_to_murinemrca <- NA

    # MRCA node to each murine
    muridmrca_to_mmu <- NA
    muridmrca_to_rno <- meas_matrix["rattus_norvegicus", murid_mrca]

  } else {
    stop()
  }

  return(data.frame(
    gerbilmrca_to_pob       = gerbilmrca_to_pob,
    gerbilmrca_to_mun       = gerbilmrca_to_mun,
    muridmrca_to_gerbilmrca = muridmrca_to_gerbilmrca,
    muridmrca_to_pob        = muridmrca_to_pob,
    muridmrca_to_mun        = muridmrca_to_mun,
    murinemrca_to_mmu       = murinemrca_to_mmu,
    murinemrca_to_rno       = murinemrca_to_rno,
    muridmrca_to_murinemrca = muridmrca_to_murinemrca,
    muridmrca_to_mmu        = muridmrca_to_mmu,
    muridmrca_to_rno        = muridmrca_to_rno
  ))
}

# ---------------------------------------------------------------------------- #

library(Biostrings)
library(ape)

gerbils <- c("psammomys_obesus",
             "meriones_unguiculatus")
murinae <- c("mus_musculus",
             "rattus_norvegicus")

muridae_species <- c("mus_musculus",
                     "rattus_norvegicus",
                     "meriones_unguiculatus",
                     "psammomys_obesus")

out_file <- "results/null_branch_lengths.csv"

# # Load index (orthogroup to gene name)
ortho_index           <- read.csv("tmp/paml_successful_gene_orthogroup.csv")

# List genes
alignment_path      <- "results/masked_alignments"
tree_path           <- "tmp/trimmed_trees"
branch_length_paths <- "tmp/branch_lengths"

alignments <- list.files(alignment_path)
trees      <- list.files(tree_path)
branches   <- list.files(branch_length_paths)

# Check the datasets have the same genes
getName <- function(my_vec) {
  sort(gsub("\\..+", "", my_vec))
}

stopifnot(getName(alignments) %in% getName(trees))
stopifnot(getName(branches) %in% getName(alignments))

stopifnot(ortho_index$orthogroup %in% getName(alignments))

# ortho_index   <- ortho_index[ortho_index$gene != "Macf1", ]

all_distances <- list()

for (i in 1:nrow(ortho_index)) {
  gene       <- as.character(ortho_index$name[i])
  orthogroup <- as.character(ortho_index$orthogroup[i])

  # Load alignment to find out the number ID of the muridae species
  species_names <- species_name_from_alignment(orthogroup, alignment_path)
  muridae_id    <- match(muridae_species, species_names)
  names(muridae_id) <- muridae_species
  # Make sure all murids are in the tree
  stopifnot(!all(is.na(muridae_id)))

  # Load tree to find the IDs of the internal nodes
  tree <- read_tree(orthogroup, tree_path)

  # Make sure tree and alignment have same species
  stopifnot(sort(tree$tip.label) == sort(species_names))

  # Double-check if muridae form a monophyletic group
  stopifnot(any(gerbils %in% tree$tip.label))
  stopifnot(any(murinae %in% tree$tip.label))

  stopifnot(is.monophyletic(tree, muridae_species))

  # Find the internal nodes linking the murids
  tree_muridae_species <- muridae_species[muridae_species %in% tree$tip.label]

  # Rename the Murid species in the trees
  murid_tree_match       <- match(muridae_species, tree$tip.label)
  murid_tree_match_no_na <- murid_tree_match[!is.na(murid_tree_match)]
  muridae_id_no_na       <- muridae_id[!is.na(murid_tree_match)]
  # tree$tip.label[murid_tree_match_no_na] <- muridae_id_no_na

  # Load "dN & dS for each branch" table, parsed from null model paml output
  branch_lengths <- read_branch_lengths(orthogroup, branch_length_paths)

  ds_matrix <- df2matrix(branch_lengths, column_name = "dS")
  ds_matrix <- rename_colnames(ds_matrix, muridae_id)
  dn_matrix <- df2matrix(branch_lengths, column_name = "dN")
  dn_matrix <- rename_colnames(dn_matrix, muridae_id)

  ds_distances <- parse_distances(meas_matrix = ds_matrix, tree, gerbils, murinae)
  dn_distances <- parse_distances(meas_matrix = dn_matrix, tree, gerbils, murinae)

  colnames(ds_distances) <- paste("dS_", colnames(ds_distances), sep = "")
  colnames(dn_distances) <- paste("dN_", colnames(dn_distances), sep = "")

  all_distances[[i]] <- data.frame(
    orthogroup = orthogroup,
    gene = gene,
    N = branch_lengths$N[1],
    S = branch_lengths$S[1],
    ds_distances,
    dn_distances)

}

all_distances_df <- do.call(rbind, all_distances)

write.csv(all_distances_df, file = out_file, row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------- #
