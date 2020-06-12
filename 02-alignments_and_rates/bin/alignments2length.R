#!/usr/bin/R

# Reads alignments from directory
# For each alignment:
  # Removes gaps for each gene
  # Measures the number of characters in the gene
# Also makes some figures

library(Biostrings)

# Input
orthogroup_seq <- "tmp/alignments/hmm_clean"
out_sizes      <- "tmp/gene_sizes/gene_sizes"
out_gaps       <- "tmp/gene_sizes/gap_sizes"
out_gaps_per   <- "tmp/gene_sizes/gap_percentages"
out_fig        <- "tmp/gene_sizes/fig"

dir.create(out_fig, recursive=TRUE)

# ---------------------------------------------------------------------------- #
# List each fasta file
fa_files <- list.files(orthogroup_seq)
fa_files <- fa_files[grep("(\\.fasta|\\.fa)", fa_files)]


ortho_widths <- list()
ortho_gaps   <- list()

for (i in 1:length(fa_files)) {
  ## Read in alignment for specific orthogroup
  fa_name      <- gsub("\\.fa(sta)*", "", fa_files[i])
  fa_file      <- file.path(orthogroup_seq, fa_files[i])
  ortho_seq    <- readDNAStringSet(fa_file)

  ## Remove gaps and measure length of each sequence
  no_gaps      <- gsub("-", "", ortho_seq)
  ortho_widths[[fa_name]] <- width(no_gaps)

  ## Measure length of gaps for each sequence
  only_gaps    <- gsub("[^-]", "", ortho_seq)
  ortho_gaps[[fa_name]]   <- width(only_gaps)
}

# ---------------------------------------------------------------------------- #

ortho_species <- lapply(ortho_widths, names)
ortho_species <- unique(unlist(ortho_species))

# ---------------------------------------------------------------------------- #
# make matrix
list2matrix <- function(ortho_list, ortho_species) {
  ortho_matrix <- matrix(-9,
                         nrow = length(ortho_list),
                         ncol = length(ortho_species))

  colnames(ortho_matrix) <- ortho_species
  rownames(ortho_matrix) <- names(ortho_list)

  # populate matrix with gene sizes
  for (i in 1:length(ortho_list)) {
    seq_widths <- ortho_list[[i]]
    ortho_matrix[i, match(names(seq_widths), ortho_species)] <- seq_widths
  }
  ortho_matrix[ortho_matrix == -9] <- NA
  return(ortho_matrix)
}

ortho_widths_matrix <- list2matrix(ortho_list    = ortho_widths,
                                   ortho_species = ortho_species)

ortho_gaps_matrix   <- list2matrix(ortho_list    = ortho_gaps,
                                   ortho_species = ortho_species)

# ---------------------------------------------------------------------------- #
## make summary plots

# Histogram of gene width for each species
for (i in 1:ncol(ortho_widths_matrix)) {
  species <- colnames(ortho_widths_matrix)[i]
  to_plot <- ortho_widths_matrix[, i]
  to_plot <- to_plot[!is.na(to_plot)]
  pdf(paste(out_fig, "/", species, ".pdf", sep = ""))
  hist(to_plot,
       main=paste("Single-copy gene length distribution in", species),
       xlab="Gene length (bp)",
       xlim=c(0, 2000),
       breaks=500)
  dev.off()
}


ortho_summ <- data.frame(
  orthogroup = rownames(ortho_widths_matrix),
  mean_size  = apply(ortho_widths_matrix, 1, function(x) mean(x, na.rm=TRUE)),
  sd_size    = apply(ortho_widths_matrix, 1, function(x) sd(x, na.rm=TRUE))
)

pdf(paste(out_fig, "/orthogroup_mean_length.pdf", sep = ""))

hist(ortho_summ$mean_size,
     main=paste("Gene length per orthogroup (xlim at 2000 bp)"),
     xlab="Mean gene length (bp)",
     xlim=c(0, 5000),
     breaks=500,
     cex.axis = 1.3,
     cex.lab = 1.7)
dev.off()

pdf(paste(out_fig, "/orthogroup_sd_length.pdf", sep = ""))

hist(ortho_summ$sd_size,
     main=paste("Gene length per orthogroup (xlim at 2000 bp)"),
     xlab="Standard deviation of gene length (bp)",
     xlim=c(0, 2000),
     breaks=500,
     cex.axis = 1.3,
     cex.lab = 1.7)
dev.off()

# ---------------------------------------------------------------------------- #
# Maximum gene size per orthogroup
ortho_summ$max <- apply(ortho_widths_matrix, 1, function(x) max(x, na.rm=TRUE))
ortho_summ$min <- apply(ortho_widths_matrix, 1, function(x) min(x, na.rm=TRUE))

pdf(paste(out_fig, "/orthogroup_max_length.pdf", sep = ""))
hist(ortho_summ$max,
     main=paste("Length of longest gene per orthogroup (xlim at 2000 bp)"),
     xlab="Maximum gene length in orthogroup (bp)",
     xlim=c(0, 2000),
     breaks=500,
     cex.axis = 1.3,
     cex.lab = 1.7)
dev.off()

pdf(paste(out_fig, "/orthogroup_min_length.pdf", sep = ""))
hist(100*ortho_summ$min/ortho_summ$max,
     main=paste("Length of smallest gene per orthogroup"),
     xlab="Gene length as percentage of longest gene in orthogroup",
     xlim=c(0, 100),
     breaks=50,
     cex.axis = 1.3,
     cex.lab = 1.7)
dev.off()

# ---------------------------------------------------------------------------- #

# Gaps as a percentage of longest gene in orthogroup
gaps_per <- lapply(1:length(ortho_summ$max),
                   function(i) 100*ortho_gaps_matrix[i,]/ortho_summ$max[i])
gaps_per <- do.call(rbind, gaps_per)

## Plot largest gap per species
ortho_summ$max_gap <- apply(gaps_per, 1, function(x) max(x, na.rm=TRUE))

pdf(paste(out_fig, "/orthogroup_max_gap.pdf", sep = ""))
hist(ortho_summ$max_gap,
     main=paste("Gene with most gaps per orthogroup"),
     xlab="Total gap length (% of longest gene)",
     xlim=c(0, 100),
     breaks=50,
     cex.axis = 1.3,
     cex.lab = 1.7)
dev.off()

## Plot largest gap per species

for (i in 1:ncol(gaps_per)) {

  species     <- colnames(gaps_per)[i]

  dir.create(file.path(out_fig, "gaps_per_species"), recursive=TRUE)
  pdf(paste(out_fig, "/gaps_per_species/", species, ".pdf", sep = ""))
  hist(gaps_per[,i],
       main=paste("Gaps in alignment in", species),
       xlab="Alignment gaps as percentage of longest gene in orthogroup",
       xlim=c(0, 100),
       breaks=50,
       cex.axis = 1.3,
       cex.lab = 1.7)
  dev.off()
}

speciesVsSpecies <- function(gaps_per, species1, species2) {
  dir.create(file.path(out_fig, "pairwise_gaps"), recursive=TRUE)
  pdf(paste(out_fig, "/pairwise_gaps/", species1, "_", species2, ".pdf", sep = ""))
  plot(gaps_per[,species2] ~ gaps_per[,species1],
       xlab = paste("Alignment gaps in", species1, "(%)"),
       ylab = paste("Alignment gaps in", species2, "(%)"),
       xlim=c(0, 100),
       ylim=c(0, 100),
       cex.axis = 1.3,
       cex.lab = 1.7)
  dev.off()
}

pairwise_combinations <- combn(ncol(gaps_per), 2)
for (i in 1:ncol(pairwise_combinations)) {
  speciesVsSpecies(gaps_per,
                   species1 = colnames(gaps_per)[pairwise_combinations[1,i]],
                   species2 = colnames(gaps_per)[pairwise_combinations[2,i]])
}

# ---------------------------------------------------------------------------- #
# Write out

makeDf <- function(ma, orthogroups) {
  df <- data.frame(ma)
  df <- data.frame(orthogroups = orthogroups, ma)
  return(df)
}

write.table(makeDf(ortho_widths_matrix,
                   orthogroups = rownames(ortho_widths_matrix)),
            file = out_sizes, sep="\t", quote=FALSE, row.names=FALSE)
write.table(makeDf(ortho_gaps_matrix,
                   orthogroups = rownames(ortho_widths_matrix)),
            file = out_gaps, sep="\t", quote=FALSE, row.names=FALSE)
write.table(makeDf(gaps_per,
                   orthogroups = rownames(ortho_widths_matrix)),
            file = out_gaps_per, sep="\t", quote=FALSE, row.names=FALSE)

# ---------------------------------------------------------------------------- #
