#!/usr/bin/R

# Input
in_ortho      <- "input/orthogroups"

# Output
out_ortho <- "tmp/filtered_orthogroups.txt"
out_count <- "tmp/filtered_orthogroups.count.txt"

# ---------------------------------------------------------------------------- #
# Load orthogroups
ortho           <- read.table(in_ortho,
                              header=TRUE,
                              sep="\t", row.names = 1, stringsAsFactors=FALSE)

# Rename species
colnames(ortho) <- gsub(".pep", "", colnames(ortho))

# Separate genes for each species in each orthgroup
# Count the number of genes for each species for each orthogroup
ortho       <- as.matrix(ortho)
gene_counts <- apply(ortho, 1,
                     function(r) sapply(strsplit(r, split = ", "), length))
gene_counts <- t(gene_counts)

# ---------------------------------------------------------------------------- #
# Filter: mouse, rat and human need to have 1 sequence each
murine_human  <- c("rattus_norvegicus", "mus_musculus", "homo_sapiens")

murine_human_count <- gene_counts[, murine_human]
murine_human_ok    <- apply(murine_human_count, 1,
                            function(x) all(x == 1))

# ---------------------------------------------------------------------------- #
# Filter: at least one gerbil has 1 sequence; neither has two
gerbil <- c("psammomys_obesus", "meriones_unguiculatus")

gerbil_count <- gene_counts[, gerbil]
gerbil_ok    <- apply(gerbil_count, 1,
                      function(x) all(x < 2) & any(x == 1))

# ---------------------------------------------------------------------------- #
# Other species: no more than two species can fewer or more than 1 sequence
other_species <- colnames(gene_counts)
other_species <- other_species[!other_species %in% c(gerbil, murine_human)]

species_count <- length(other_species)

other_species_count <- gene_counts[, other_species]
other_species_ok    <- apply(other_species_count, 1,
                             function(x) sum(x == 1) >= (species_count - 2))

# ---------------------------------------------------------------------------- #
# Join all filters

total_filter <- murine_human_ok & gerbil_ok & other_species_ok

# ---------------------------------------------------------------------------- #

# If any species has a multiple sequences for the same orthology group, mask the species
ortho[gene_counts != 1] <- ""

# ---------------------------------------------------------------------------- #
# Apply filter to matrix of orthogroup sequences
ortho_pass                   <- ortho[total_filter, ]
ortho_pass[ortho_pass == ""] <- NA

# Apply filter to matrix of sequence counts
gene_count_filtered <- gene_counts[total_filter, ]

# ---------------------------------------------------------------------------- #

## Write out
write.table(gene_count_filtered,
            file = out_count, sep="\t", quote=FALSE)

write.table(ortho_pass, file = out_ortho, sep="\t", quote=FALSE)

# ---------------------------------------------------------------------------- #
