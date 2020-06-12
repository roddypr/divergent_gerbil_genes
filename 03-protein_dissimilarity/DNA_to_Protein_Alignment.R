#!/usr/bin/Rscript

#Script for parsing out gaps in DNA alignment while keeping integrity of entire alignment
#Script to then translate parsed DNA alignment and retain only sites that are conserved in two selected outgroups

library(stringr)
library(Biostrings)

#Create function for identifying outgroups
get_outgroup_species <- function(species_vec) {
	outgroup_alternatives <- list(
		c("cavia_porcellus", "homo_sapiens"),
	        c("homo_sapiens", "chinchilla_lanigera"),
	        c("homo_sapiens", "ictidomys_tridecemlineatus"),
		c("homo_sapiens", "octodon_degus"),
		c("homo_sapiens", "peromyscus_maniculatus"),
		c("oryctolagus_cuniculus", "ictidomys_tridecemlineatus"),
		c("homo_sapiens", "jaculus_jaculus"))

	for(i in 1:length(outgroup_alternatives)) {
		if (all(outgroup_alternatives[[i]] %in% species_vec)) {
			outgroup_species <- outgroup_alternatives[[i]]
			break()
		}
	}

	if (is.null(outgroup_species)) {
		stop("No outgroups found")
	}
	return(outgroup_species)
}

#Create function to divide DNA alignment into codons (chunks of three nucleotides, including gaps)

divide_into_codons <- function(my_seq) {
	start <- seq(1,nchar(my_seq),3)
	stop  <- pmin(start + 2, nchar(my_seq))
	return(str_sub(my_seq, start, stop))
}

#Create function to remove all codons that contain gaps
#Ensure this deletion happens stimutaneously for all species (i.e. ensure alignment still is in place)


remove_ns <- function(dna_alignment) {
	require(stringr)
	dna_alignment    <- as.character(dna_alignment)
	alignment_codons <- lapply(dna_alignment, divide_into_codons)
	gap_indeces      <- lapply(alignment_codons, function(x) grep("N", x))
	gap_indeces      <- sort(unique(unlist(gap_indeces)))

	if (length(gap_indeces > 0)) {
		no_gap_alignment <- lapply(alignment_codons, function(x) x[-gap_indeces])
	} else {
		no_gap_alignment <- alignment_codons
	}
	no_gap_alignment <- lapply(no_gap_alignment, function(x) paste(x, collapse = ""))
	return(DNAStringSet(unlist(no_gap_alignment)))
}

remove_gaps <- function(dna_alignment) {
	require(stringr)
	dna_alignment    <- as.character(dna_alignment)
	alignment_codons <- lapply(dna_alignment, divide_into_codons)
	gap_indeces      <- lapply(alignment_codons, function(x) grep("\\-", x))
	gap_indeces      <- sort(unique(unlist(gap_indeces)))

	if (length(gap_indeces > 0)) {
		no_gap_alignment <- lapply(alignment_codons, function(x) x[-gap_indeces])
	} else {
		no_gap_alignment <- alignment_codons
	}
	no_gap_alignment <- lapply(no_gap_alignment, function(x) paste(x, collapse = ""))
	return(DNAStringSet(unlist(no_gap_alignment)))
}

remove_stop_codons <- function(aa_alignment) {
  require(stringr)
  aa_alignment_split <- strsplit(as.character(aa_alignment), split = "")
  stop_indeces       <- lapply(aa_alignment_split, function(x) grep("\\*", x))
  stop_indeces       <- sort(unique(unlist(stop_indeces)))

  if (length(stop_indeces > 0)) {
    no_stop_alignment <- lapply(aa_alignment_split, function(x) x[-stop_indeces])
  } else {
    no_stop_alignment <- aa_alignment_split
  }

  no_stop_alignment <- lapply(no_stop_alignment, function(x) paste(x, collapse = ""))
  return(AAStringSet(unlist(no_stop_alignment)))
}

# ------------------------------------------------------------------------------------------ #

input_directory  <- "DNA_alignments"
output_directory <- "protein_alignments"

alignment_files  <- list.files(input_directory)

if (!dir.exists(input_directory)) {
	stop(paste("Output directory does not exist:\n", output_directory, sep = ""))
}

if (!dir.exists(output_directory)) {
	stop(paste("Output directory does not exist:\n", output_directory, sep = ""))
}

for (gene_file in alignment_files) {

	gene_path <- file.path(input_directory, gene_file)
	gene_name <- gsub(".fa", "", gene_file)

	dna_alignment <- readDNAStringSet(gene_path)
	dna_alignment <- remove_ns(dna_alignment)
	dna_alignment <- remove_gaps(dna_alignment)
 	aa_alignment  <- translate(dna_alignment)
 	aa_alignment  <- remove_stop_codons(aa_alignment)

	# Identify outgroups
	outgroup_species   <- get_outgroup_species(names(dna_alignment))

	# Identify the positions where both outgroups have the same AA
	outgroup_alignment <- aa_alignment[outgroup_species]

	# split into list of character vectors
	outgroup_alignment <- strsplit(as.character(outgroup_alignment), split = "")

	conserved_sites <- which(outgroup_alignment[[1]] == outgroup_alignment[[2]])

	# go back to alignment and subset to conserved sites
	aa_alignment_split <- strsplit(as.character(aa_alignment), split = "")

	aa_alignment_subset <- lapply(aa_alignment_split, function(seq) seq[conserved_sites])
	aa_alignment_subset <- lapply(aa_alignment_subset, function(seq) paste(seq, collapse = ""))
	aa_alignment_subset <- unlist(aa_alignment_subset)

	aa_alignment_subset <- AAStringSet(aa_alignment_subset)
	out_file <- file.path(output_directory, paste(gene_name, ".fasta", sep = ""))
 	writeXStringSet(aa_alignment_subset, file = out_file)

}
