#!/usr/bin/R

# Input

## Functions
df2matrix <- function(df) {
  ma           <- as.matrix(df[, -1])
  rownames(ma) <- df$orthogroup
  stopifnot(apply(ma, 2, is.numeric))
  return(ma)

}

## READ IN ALIGNMENT NAMES
orthogroup_seq <- "tmp/alignments/hmm_clean"
fa_files       <- list.files(orthogroup_seq)
fa_files       <- fa_files[grep(".fa", fa_files)]

## READ IN GENE SIZE IN ALIGNMENTS
in_sizes       <- "tmp/gene_sizes/gene_sizes_edited"
gene_sizes     <- read.table(in_sizes, header = TRUE)
gene_sizes     <- df2matrix(gene_sizes)

## READ IN GAP SIZE IN ALIGNMENTS
in_gaps_per    <- "tmp/gene_sizes/gap_percentages_edited"
gaps_per       <- read.table(in_gaps_per, header = TRUE)
gaps_per       <- df2matrix(gaps_per)

## Ingroup and outgroup species
gerbil       <- c("psammomys_obesus", "meriones_unguiculatus")
murine_human <- c("rattus_norvegicus", "mus_musculus", "homo_sapiens")

## Output directory
out_alignment_dir <- "results/masked_alignments"
dir.create(out_alignment_dir, recursive = TRUE)

out_file <- "results/masked_alignments_index"

# ---------------------------------------------------------------------------- #
# first eliminate alignments where the longest gene was smaller than 400 bp.
max_gene_size    <- apply(gene_sizes, 1, function(r) max(r, na.rm = TRUE))
max_gene_size_ok <- max_gene_size >= 400

gaps_per <- gaps_per[max_gene_size_ok, ]

# ---------------------------------------------------------------------------- #
# Removal summary

summary_class <- c("original",
                   "removed_400_and_shorter",
                   "longer_400_remaining",
                   "genes_with_gaps")
summary_numbers <- c(length(max_gene_size),
                     sum(!max_gene_size_ok),
                     sum(max_gene_size_ok),
           sum(apply(gaps_per, 1, function(x) any(x > 40, na.rm=TRUE))))

# ---------------------------------------------------------------------------- #
# masked any sequence for which gaps represented more than 70% of the non-gap
# size of the longest gene in the alignment.

gaps_per[gaps_per > 70] <- NA

# ---------------------------------------------------------------------------- #
# Filter: mouse, rat and human need to have 1 sequence each

murine_human_per <- gaps_per[, murine_human]
murine_human_ok  <- apply(murine_human_per, 1, function(r) all(!is.na(r)))

# ---------------------------------------------------------------------------- #
# Filter: at least one gerbil has 1 sequence

gerbil_per <- gaps_per[, gerbil]
gerbil_ok  <- apply(gerbil_per, 1, function(r) sum(!is.na(r)) > 0)

# For stats only:
gerbil_1_masked  <- apply(gerbil_per, 1, function(r) sum(is.na(r)) == 1)

# ---------------------------------------------------------------------------- #
# Other species: no more than two species can fewer or more than 1 sequence

special_species     <- c(gerbil, murine_human)
other_species_per   <- gaps_per[, !colnames(gaps_per) %in% special_species]

species_count       <- ncol(other_species_per)

other_species_ok    <- apply(other_species_per, 1,
                             function(r) sum(!is.na(r)) >= (species_count - 2))

# ---------------------------------------------------------------------------- #
# Join all filters

total_filter <- murine_human_ok & gerbil_ok & other_species_ok

gaps_per <- gaps_per[total_filter, ]

# ---------------------------------------------------------------------------- #
# Removal summary

summary_class   <- c(summary_class,
                     "Murine_human_have_masked",
                     "Gerbils_have_more_than_1_masked",
                     "Other_species_more_than_2_masked",
                     "Total_masking_issues",
                     "OK_but_gerbils_1_masked",
                     "Final_remaining")

summary_numbers <- c(summary_numbers,
                     sum(!murine_human_ok),
                     sum(!gerbil_ok),
                     sum(!other_species_ok),
                     sum(!total_filter),
                     sum(gerbil_1_masked),
                     nrow(gaps_per))

summary_table <- data.frame(class            = summary_class,
                            alignment_number = summary_numbers)

write.table(summary_table,
            file      = "results/orthogroup_filter_summary.txt",
            row.names = FALSE,
            quote     = FALSE)

# ---------------------------------------------------------------------------- #
## Apply the masking to the alignments
library(Biostrings)

fa_files_orthogroup <- gsub(".fa", "", fa_files)

to_keep <- fa_files_orthogroup %in% rownames(gaps_per)

fa_files_orthogroup <- fa_files_orthogroup[to_keep]
fa_files            <- fa_files[to_keep]

# Same order for fa_files and matrix
matrix_file_match  <- match(rownames(gaps_per), fa_files_orthogroup)
fa_files            <- fa_files[matrix_file_match]
fa_files_orthogroup <- fa_files_orthogroup[matrix_file_match]


# For each alignment
for (i in 1:length(fa_files)) {
  # Read in fasta file
  fa_name      <- fa_files[i]
  ortho_name   <- gsub("\\.fa", "", fa_name)
  fa_file      <- file.path(orthogroup_seq, fa_name)
  ortho_seq    <- readDNAStringSet(fa_file)

  # non-masked species
  gaps_per_i         <- gaps_per[i,]
  non_masked_species <- names(gaps_per_i)[!is.na(gaps_per_i)]

  # Include only non-masked species in orthogroup
  ortho_seq <- ortho_seq[names(ortho_seq) %in% non_masked_species]

  # Write out
  alignment_path <- file.path(out_alignment_dir,
                              paste(ortho_name, ".fa", sep = ""))
  writeXStringSet(ortho_seq, file=alignment_path)
}

# ---------------------------------------------------------------------------- #
# Make index
index_ma <- gaps_per
index_ma[!is.na(index_ma)] <- 1
index_ma[is.na(index_ma)]  <- 0

stopifnot(names(table(index_ma)) == c(0, 1))

index_df <- data.frame(orthogroup=rownames(index_ma),
                       data.frame(index_ma))

## Write out
write.table(index_df,
            file = out_file, sep="\t", quote=FALSE, row.names=FALSE)

# ---------------------------------------------------------------------------- #
