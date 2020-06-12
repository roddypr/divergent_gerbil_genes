#!/usr/bin/R

# ---------------------------------------------------------------------------- #
# For each gene lister in input/longest_tx_ids:
#    make a new fasta file, where the sequences for all species are printed

# ---------------------------------------------------------------------------- #

# Input
ortho             <- "tmp/filtered_orthogroups.txt"

gene_to_tx_to_pep <- "input/longest_tx_ids"

cds_species_path  <- "tmp/annotations/cds"

OUTDIR            <- "tmp/gene_sequences_per_orthogroup/"

# ---------------------------------------------------------------------------- #
# Load pre-filtered orthogroups
ortho <- read.table(ortho,
                    header=TRUE,
                    sep="\t", row.names = 1, stringsAsFactors=FALSE)

# ---------------------------------------------------------------------------- #
# Load table with Gene -> transcript -> CDS Names
gene_to_tx_to_pep <- read.table(gene_to_tx_to_pep,
                                sep = ' ',
                                stringsAsFactors=FALSE)
colnames(gene_to_tx_to_pep) <- c("species", "gene", "tx", "pep", "symbol")

# ---------------------------------------------------------------------------- #
# Make sure gene version numbers are removed from the datasets
gene_to_tx_to_pep$gene <- gsub("\\.[0-9]+", "", gene_to_tx_to_pep$gene)
gene_to_tx_to_pep$tx   <- gsub("\\.[0-9]+", "", gene_to_tx_to_pep$tx)
gene_to_tx_to_pep$pep  <- gsub("\\.[0-9]+", "", gene_to_tx_to_pep$pep)

ortho <- apply(ortho, 2, function(x) gsub("\\.[0-9]", "", x))

stopifnot(apply(ortho, 1, function(x) !any(grepl("\\.[0-9]", x))))

# ---------------------------------------------------------------------------- #
# Make ortho table of CDS ID (Tx) for the orthogroups in the rat

ortho_tx <- apply(ortho, 2,
      function(x) gene_to_tx_to_pep$tx[match(x, gene_to_tx_to_pep$pep)])

rownames(ortho_tx) <- rownames(ortho)

# ---------------------------------------------------------------------------- #
## Load the FASTA sequence for the CDS of each rodent species
#     (including gerbils)

library(Biostrings)
fa_list <- list()
for (species in colnames(ortho)) {
  species_path <- paste(cds_species_path, "/", species, ".fa", sep="")
  species_fa   <- readDNAStringSet(species_path)
  fa_list[[species]] <- species_fa
}
fa_list <- DNAStringSetList(fa_list)

# ---------------------------------------------------------------------------- #
## For each species, loop through the orthogroups, and append the gene sequence
# to the orthogroup fasta file

orthogroup_list <- list()
dir.create(OUTDIR)
stopifnot(colnames(ortho_tx) %in% names(fa_list))


for (species in colnames(ortho_tx)) {

  species_fa    <- fa_list[[species]]
  species_ortho <- ortho_tx[,species]

  for (i in 1:nrow(ortho_tx)) {

    orthogroup_id <- rownames(ortho_tx)[i]

    if (!is.na(species_ortho[i])) {
      new_file      <- paste(OUTDIR, orthogroup_id, ".fa", sep="")
      new_fa        <- species_fa[species_ortho[i]]
      names(new_fa) <- species
      writeXStringSet(new_fa, filepath=new_file, append=TRUE)

      # Add to list to make INDEX
      orthogroup_list[[orthogroup_id]] <- c(orthogroup_list[[orthogroup_id]],
                                            species)
    }
  }
}


# Make INDEX
index_ma <- matrix(nrow = length(orthogroup_list),
                   ncol = length(unique(sort(unlist(orthogroup_list)))))
colnames(index_ma) <- unique(sort(unlist(orthogroup_list)))
for (i in 1:nrow(index_ma)) {
  index_ma[i, orthogroup_list[[i]]] <- 1
}

index_ma[is.na(index_ma)] <- 0

index_df <- cbind(orthogroup = names(orthogroup_list),
                  data.frame(index_ma))

write.table(index_df,
            file      = paste(OUTDIR, "/INDEX", sep = ""),
            quote     = FALSE,
            sep       = "\t",
            row.names = FALSE)
