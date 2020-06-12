#!/usr/bin/Rscript

library(stringr)
library(Biostrings)

#Create function for identifying outgroups
get_outgroup_species <- function(species_vec) {
	outgroup_alternatives <- c("homo_sapiens", "oryctolagus_cuniculus")

	for(i in 1:length(outgroup_alternatives)) {
		if (outgroup_alternatives[i] %in% species_vec) {
			outgroup_species <- outgroup_alternatives[i]
			break()
		}
	}

	if (is.null(outgroup_species)) {
		stop("No outgroups found")
	}
	return(outgroup_species)
}

check_sequence_lengths <- function(split_alignment, file_name) {

  seq_lengths <- sapply(split_alignment, length)

  if (any(seq_lengths != seq_lengths[1])) {
    stop(paste("Alignment has sequences of different lengths:", file_name))
  }
}

#Create function to divide DNA alignment into codons
#  (chunks of three nucleotides, including gaps)

divide_into_codons <- function(my_seq) {
	start <- seq(1,nchar(my_seq),3)
	stop  <- pmin(start + 2, nchar(my_seq))
	return(str_sub(my_seq, start, stop))
}

score_per_sequence <- function(outgroup, ingroup, score_matrix) {
  stopifnot(length(ingroup) == length(outgroup))
  score_per_aa   <- score_matrix[cbind(ingroup, outgroup)]
  sequence_score <- sum(score_per_aa)
  return(sequence_score)

}

calculate_score <- function(split_alignment, outgroup, score_matrix) {

  # Get name of ingroup species
  species         <- names(split_alignment)
  ingroup_species <- species[species != outgroup_species]

  # Get outgroup sequence
  outgroup_seq   <- split_alignment[[outgroup]]

  # For each ingroup, calculate score
  ingroup_scores <- lapply(ingroup_species,
    function(x)
        score_per_sequence(outgroup = outgroup_seq,
                           ingroup  = split_alignment[[x]],
                           score_matrix))

  ingroup_scores        <- unlist(ingroup_scores)
  names(ingroup_scores) <- ingroup_species
  return(ingroup_scores)

}

# ---------------------------------------------------------------------------- #
# Protein dissimilarity coefficient matrices

sneath_matrix  <- read.table("Sneaths_index.txt", header = TRUE)
epstein_matrix <- read.table("Epsteins_coefficient.txt", header = TRUE)

# ---------------------------------------------------------------------------- #

input_directory  <- "protein_alignments"
output_directory <- "results"

alignment_files  <- list.files(input_directory)

if (!dir.exists(input_directory)) {
  stop(paste("Output directory does not exist:\n", output_directory, sep = ""))
}

if (!dir.exists(output_directory)) {
  stop(paste("Output directory does not exist:\n", output_directory, sep = ""))
}

results_list <- list()

for (gene_file in alignment_files) {

  gene_path <- file.path(input_directory, gene_file)
  gene_name <- gsub(".fasta", "", gene_file)

  aa_alignment  <- readAAStringSet(gene_path)

  # Identify outgroups
  outgroup_species <- get_outgroup_species(names(aa_alignment))

  # Split
  aa_split <- strsplit(as.character(aa_alignment), split = "")

  # Check all sequences have the same size
  check_sequence_lengths(split_alignment = aa_split, file_name = gene_file)

  # Sneath score
  sneath_score <- calculate_score(split_alignment = aa_split,
                                  outgroup        = outgroup_species,
                                  score_matrix    = sneath_matrix)

  # Epstein score
  epstein_score <- calculate_score(split_alignment = aa_split,
                                   outgroup        = outgroup_species,
                                   score_matrix    = epstein_matrix)

  # Rename
  names(sneath_score)  <- paste("sneath_", names(sneath_score), sep = "")
  names(epstein_score) <- paste("epstein_", names(epstein_score), sep = "")

  # Store in list
  results_list[[gene_name]] <- c(sneath_score, epstein_score)

}


# Check all column names
total_colnames <- unique(unlist(lapply(results_list, names)))
results_list   <- lapply(results_list, function(x) x[total_colnames])

# Now we need to turn the lists to data frames
results_df <- do.call(rbind, results_list)
results_df <- data.frame(gene = names(results_list), results_df)

out_file <- file.path(output_directory, "protein-dissimilarity.csv")
write.csv(results_df, file = out_file, row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------- #
