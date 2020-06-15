# 1-to-1 orthologs between rodents and human

We want to find 1-to-1 orthologous protein-coding genes between the genomes of two gerbils, 10 other rodents and *H. sapiens*.

## Input

The directory `input` three folders, one for the cds fasta files for each species (`genus_species.cds.fa.gz`), one for the peptide fasta files (`genus_species.pep.fa.gz`) and one for the gff files (`genus_species.gff3.gz`), as downloaded from ensembl.

```sh

# Add the cds fasta file for all species in the format genus_species.cds.fa.gz
mkdir -p input/annotations/cds
ls input/annotations/cds | head -n2
# cavia_porcellus.cds.fa.gz
# chinchilla_lanigera.cds.fa.gz

mkdir -p input/annotations/pep
ls input/annotations/pep | head -n2
# cavia_porcellus.pep.fa.gz
# chinchilla_lanigera.pep.fa.gz

# Add the gff3 file for all species in the format genus_species.gff3.gz
mkdir -p input/annotations/gff
ls input/annotations/gff | head -n2
# cavia_porcellus.gff3.gz
# chinchilla_lanigera.gff3.gz

```

## Process

### Parse the annotation for each species

We start by getting the longest transcript per gene. Note that in `longest_transcript_per_gene.nf`, the fat sand rat is ignored, as its gff file is formatted in a different way.

```sh

conda activate nextflow
nextflow run \
  -w tmp/longest_tx bin/longest_transcript_per_gene.nf \
  1> tmp/longest_tx.out 2> tmp/longest_tx.err &

```

We now get the summary of the gff file for the fat sand rat.

```sh

zcat input/annotations/gff/psammomys_obesus.gff3.gz \
  | awk '$3 == "mRNA" {print $9}' > tmp/psammomys_obesus.mrna

conda activate r
Rscript bin/psammomys_obesus_gff_summary.R

cat tmp/psammomys_obesus.index >> results/longest_tx_fa/INDEX

zcat input/annotations/pep/psammomys_obesus.pep.fa.gz \
  | ruby -pe 'gsub(/ .+/, "")' \
  > results/longest_tx_fa/psammomys_obesus.pep.fa

grep -c ">" results/longest_tx_fa/*
# results/longest_tx_fa/cavia_porcellus.pep.fa:18095
# results/longest_tx_fa/chinchilla_lanigera.pep.fa:17809
# results/longest_tx_fa/dipodomys_ordii.pep.fa:16911
# results/longest_tx_fa/fukomys_damarensis.pep.fa:17730
# results/longest_tx_fa/homo_sapiens.pep.fa:19951
# results/longest_tx_fa/ictidomys_tridecemlineatus.pep.fa:18474
# results/longest_tx_fa/jaculus_jaculus.pep.fa:17845
# results/longest_tx_fa/meriones_unguiculatus.pep.fa:21116
# results/longest_tx_fa/mus_musculus.pep.fa:22026
# results/longest_tx_fa/nannospalax_galili.pep.fa:18647
# results/longest_tx_fa/octodon_degus.pep.fa:19982
# results/longest_tx_fa/psammomys_obesus.pep.fa:21807
# results/longest_tx_fa/rattus_norvegicus.pep.fa:22250

zgrep -c ">" input/annotations/pep/*
# input/annotations/pep/cavia_porcellus.pep.fa.gz:25582
# input/annotations/pep/chinchilla_lanigera.pep.fa.gz:25983
# input/annotations/pep/dipodomys_ordii.pep.fa.gz:23962
# input/annotations/pep/fukomys_damarensis.pep.fa.gz:23413
# input/annotations/pep/homo_sapiens.pep.fa.gz:109095
# input/annotations/pep/ictidomys_tridecemlineatus.pep.fa.gz:25958
# input/annotations/pep/jaculus_jaculus.pep.fa.gz:25180
# input/annotations/pep/meriones_unguiculatus.pep.fa.gz:31063
# input/annotations/pep/mus_musculus.pep.fa.gz:66760
# input/annotations/pep/nannospalax_galili.pep.fa.gz:26872
# input/annotations/pep/octodon_degus.pep.fa.gz:27047
# input/annotations/pep/psammomys_obesus.pep.fa.gz:21807
# input/annotations/pep/rattus_norvegicus.pep.fa.gz:29107

```

### Species tree

A species tree is in `input/2019-02-11-tree.newick`.

### Run OrthoFinder

```sh
mv results/longest_tx_fa/INDEX results/longest_tx_ids

conda activate orthofinder
orthofinder -t 30 \
  -f results/longest_tx_fa \
  -S diamond \
  -s input/2019-02-11-tree.newick \
  1> tmp/orthofinder.out 2> tmp/orthofinder.err

conda deactivate

```

Move orthofinder output to the `results` directory.

```sh

mv results/longest_tx_fa/Results_May13 results/orthofinder

```

I copied the following from the file `Statistics_Overall.csv`.

| Category | Number |
| ---------| ------ |
| Number of genes | 252643 |
| Number of genes in | 241338 |
| Number of unassigned genes | 11305 |
| Percentage of genes in orthogroups | 95.5 |
| Percentage of unassigned gene | 4.5 |
| Number of orthogroups | 18321 |
| Number of species-specific orthogroups | 26 |
| Number of genes in species-specific orthogroups | 112 |
| Percentage of genes in species-specific orthogroups | 0.0 |
| Mean orthogroup size | 13.2 |
| Median orthogroup size | 13.0 |
| G50 (assigned genes) | 13 |
| G50 (all genes) | 13 |
| O50 (assigned genes) | 7030 |
| O50 (all genes) | 7465 |
| Number of orthogroups with all species present | 10559 |
| Number of single-copy orthogroups | 7053 |


Although there are 10559 orthogroups with genes from all species, there is a large number of orthogroups with missing a gene in a single species. We need to make sure that there isn't a single species that is contributing to this result.

| Number of species in orthogroup	| Number of orthogroups
| --- | ---
| 1 | 26 |
| 2 | 743 |
| 3 | 346 |
| 4 | 334 |
| 5 | 273 |
| 6 | 232 |
| 7 | 239 |
| 8 | 286 |
| 9 | 418 |
| 10 | 559 |
| 11 | 1081 |
| 12 | 3225 |
| 13 | 10559 |

## Summary figures

```r

gene_count <- "results/orthofinder/Orthogroups.GeneCount.csv"
gene_count <- read.table(file = gene_count, header= T)

# m:..:n:0
one_gene_missing <- apply(gene_count[,-ncol(gene_count)],
                          1,
                          function(r) sum(r == 0) == 1)
sum(one_gene_missing)
# 3272

one_gene_missing <- gene_count[one_gene_missing,-ncol(gene_count)]
which_missing    <- apply(one_gene_missing, 1, function(r) which(r == 0))
which_missing    <- table(colnames(one_gene_missing)[which_missing])
which_missing    <- data.frame(which_missing)
colnames(which_missing) <- c("Species", "Frequency")

which_missing$Species <- gsub(".pep", "", which_missing$Species)
which_missing$Species <- gsub("_", " ", which_missing$Species)
library(Hmisc)
which_missing$Species <- capitalize(which_missing$Species)

pdf("results/m_n_0_missing.pdf")
library(ggplot2)
ggplot(which_missing) +
  geom_bar(stat="identity", aes(x = Species, y = Frequency)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()



# 1:..:1:0
one_gene_missing <- apply(gene_count[,-ncol(gene_count)],
                          1,
                          function(r) sum(r == 0) == 1 &
                          sum(r == 1) == (ncol(gene_count) - 2))
sum(one_gene_missing)
# 2369
one_gene_missing <- gene_count[one_gene_missing,-ncol(gene_count)]
which_missing    <- apply(one_gene_missing, 1, function(r) which(r == 0))
which_missing    <- table(colnames(one_gene_missing)[which_missing])
which_missing    <- data.frame(which_missing)
colnames(which_missing) <- c("Species", "Frequency")

which_missing$Species <- gsub(".pep", "", which_missing$Species)
which_missing$Species <- gsub("_", " ", which_missing$Species)
library(Hmisc)
which_missing$Species <- capitalize(which_missing$Species)


pdf("results/1_1_0_missing.pdf")
library(ggplot2)
ggplot(which_missing) +
  geom_bar(stat="identity", aes(x = Species, y = Frequency)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()

```
