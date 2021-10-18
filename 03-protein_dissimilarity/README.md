# Dissimilarity measurements

The scripts in this directory use edit and translate the DNA alignment for each group of orthologous protein-coding genes between rodent species, and then measure two dissimilarity measures between each sequence and three outgroups.

## Edit each alignment

We get each alignment from the output of `../02-alignments_and_rates`. Each alignment should be in a file called `gene_id.fa` (where `gene_id` is an arbitrary unique name).

```sh

mkdir -p DNA_alignments

# Then copy the alignments to the directory

```

In the following script, we translate the alignments and edit them, keeping only the positions that are identitical in three outgroup species. The output is stored in `protein_alignments`.

```sh

mkdir -p protein_alignments

Rscript DNA_to_Protein_Alignment.R

```

## Measure the dissimilarity values for each alignment

In the following script, we compare the amino-acid residues of each rodent to those of *H. sapiens* and score them according to the matrices in `Sneaths_index.txt` and `Epsteins_coefficient.txt`. The output is written to `results/protein-dissimilarity.csv`

```sh

mkdir -p results

Rscript protein_dissimilarity.R

```
