#!/usr/bin/R

# Input
in_ortho      <- "input/orthogroups"
in_gene_names <- "input/longest_tx_ids"
t2d_genes     <- "tmp/t2d_genes"

# Output
out_ortho <- "tmp/filtered_orthogroups.txt"

# Filter threshold
maximum_indel <- 2

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
gene_names <- read.table(in_gene_names, sep = " ")

colnames(gene_names) <- c("species", "gene", "tx", "protein", "gene_name")

gene_names$gene_name <- gsub("-[0-9]+$", "", gene_names$gene_name)
gene_names_mouse     <- gene_names[gene_names$species == "mus_musculus", ]

t2d_genes            <- as.character(read.table(t2d_genes)$V1)
t2d_genes_prot       <- gene_names_mouse$protein[match(t2d_genes, gene_names_mouse$gene_name)]
t2d_genes_prot       <- as.character(t2d_genes_prot)

stopifnot(!is.na(t2d_genes_prot))

# ---------------------------------------------------------------------------- #
# Mouse groups
mouse_ortho <- ortho[, 'mus_musculus']
mouse_ortho <- strsplit(mouse_ortho, split = ", ")

t2d_in_ortho <-
sapply(t2d_genes_prot, function(t2d_gene)
  sapply(mouse_ortho, function(x) t2d_gene %in% x))

stopifnot(colnames(t2d_in_ortho) == t2d_genes_prot)
colnames(t2d_in_ortho) <- t2d_genes

t2d_in_ortho <- t2d_in_ortho[apply(t2d_in_ortho, 1, any), ]
t2d_genes_per_orthogroup <- apply(t2d_in_ortho, 1, sum)

## Present multiple times:
which_2 <- colnames(t2d_in_ortho)[t2d_in_ortho[t2d_genes_per_orthogroup > 1, ]]

t2d_orthogroup_which <- apply(t2d_in_ortho, 2, which)

t2d_orthogroup        <- rownames(t2d_in_ortho)[t2d_orthogroup_which]
names(t2d_orthogroup) <- names(t2d_orthogroup_which)


# ---------------------------------------------------------------------------- #
# Find which genes in Orthogroup list are part of the T2D list

td2_gene_counts <- gene_counts[match(t2d_orthogroup, rownames(gene_counts)), ]

chosen_species <- c("homo_sapiens",
                    "rattus_norvegicus",
                    "mus_musculus",
                    "psammomys_obesus",
                    "meriones_unguiculatus")

all_fine <- apply(td2_gene_counts ,1 , function(x) all(x[chosen_species] == 1))
all_fine_ortho <- names(all_fine[all_fine])
all_fine_gene <- names(t2d_orthogroup[t2d_orthogroup %in% all_fine_ortho])

# For each T2D genes, which have a missing sequence for one of the ingroup or
#    outgroup species

one_is_zero       <- apply(td2_gene_counts,
                           1,
                           function(x) any(x[chosen_species] == 0))

one_is_zero_ortho <- names(one_is_zero[one_is_zero])
one_is_zero_gene  <- names(t2d_orthogroup[t2d_orthogroup %in% one_is_zero_ortho])

table_of_zeros <- td2_gene_counts[one_is_zero_ortho, ]
zero_count     <- apply(table_of_zeros, 1, function(x) sum(x == 0))
inzero_count   <- apply(table_of_zeros[, chosen_species],
                        1,
                        function(x) sum(x == 0))

inzero_which   <- apply(table_of_zeros[, chosen_species],
                        1,
                        function(x) paste(names(x[x == 0]), collapse = ","))

## non-usable because of zero
names(t2d_orthogroup[t2d_orthogroup %in% names(zero_count[zero_count >= 5])])
# [1] "Irs3"  "Zfp57" "Ddit3" "Irs2"

write.table(data.frame(gene = one_is_zero_gene,
                       orthogroup = one_is_zero_ortho,
                       zero_count = zero_count
                       ),
            file="tmp/t2d_one2one_ortho.txt",
            row.names=FALSE, col.names=FALSE, quote=FALSE)

length(all_fine_gene)
# 80

# ---------------------------------------------------------------------------- #


ingroup_duplicated <- apply(td2_gene_counts ,1 , function(x) any(x[chosen_species] > 1))
ingroup_duplicated_ortho <- names(ingroup_duplicated[ingroup_duplicated])
ingroup_duplicated_gene <- names(t2d_orthogroup[t2d_orthogroup %in% ingroup_duplicated_ortho])

table_of_duplicateds <- td2_gene_counts[ingroup_duplicated_ortho, ]
duplicated_count <- apply(table_of_duplicateds, 1, function(x) sum(x > 1))
induplicated_count <- apply(table_of_duplicateds[, chosen_species], 1, function(x) sum(x > 1))
induplicated_which <- apply(table_of_duplicateds[, chosen_species], 1, function(x)
  paste(names(x[x > 1]), collapse = ","))


problem_df <- rbind(
data.frame(gene = all_fine_gene, orthogroup = all_fine_ortho,
           issue = "none", number_ingroup_affected = NA,
           number_species_affected = NA,
           which_ingroup_affected = NA),
data.frame(gene = c(one_is_zero_gene, ingroup_duplicated_gene),
          orthogroup = c(one_is_zero_ortho, ingroup_duplicated_ortho),
          issue = c(rep("deletion", length(zero_count)),
                    rep("expansion", length(duplicated_count))),
          number_ingroup_affected = c(inzero_count, induplicated_count),
          number_species_affected = c(zero_count, duplicated_count),
          which_ingroup_affected = c(inzero_which, induplicated_which)))


write.csv(problem_df, file="results/t2d_orthology_summary.csv")

# # 80
# # ---------------------------------------------------------------------------- #
#
# # Filter
#    # OK if for all species gene count is one
#    # OK if gene count for maximum of one species is not one
#       # but first remove that gene
#
# ortho[gene_counts != 1] <- ""
#
# ok_gene_counts <- apply(gene_counts, 1,
#       function(r) sum(r == 1) >= (length(r) - maximum_indel))
#
# ortho_pass     <- ortho[ok_gene_counts, ]
#
# ortho_pass[ortho_pass == ""] <- NA
#
#
# # ---------------------------------------------------------------------------- #
#
# write.table(ortho_pass, file = out_ortho, sep="\t", quote=FALSE)
#
# # ---------------------------------------------------------------------------- #
