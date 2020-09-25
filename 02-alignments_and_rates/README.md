# Align orthogroups and measure substitution rates

The scripts in this directory perform the following steps:
* align the sequences of each group of orthologous genes
* filter the resulting alignments
* measure the substitution rates per mutational category with bioppml
* run branch test and branch-site test with PAML and godon
* retrieve the measurements for gerbil-murine divergence

## Input

We first get the groups of orthologous sequences:

```sh


ln -sfr ../01-orthogroups/results/longest_tx_ids input

ln -sfr \
  ../01-orthogroups/results/orthofinder/Orthogroups.csv \
  input/orthogroups

# grep mus_musculus input/longest_tx_ids > tmp/mouse_names

```

We then get and de-compress the CDS sequences for the rodents and the human.


```sh

mkdir -p input/annotations/gff
mkdir -p input/annotations/cds
mkdir -p input/annotations/pep

mkdir -p tmp

ANN_DIR=../input_data/2019-01-17-input/rodent_genome_annotations/
# This directory includes one folder per species with a cds fasta file,
#    and a pep fasta and a gff file (as downloaded from ensembl)

## Remove unwanted species
ls ${ANN_DIR} \
  | grep -v "ochotona_princeps" | grep -v "oryctolagus_cuniculus" \
  | cut -f 5 -d '/' > tmp/file_list

while read p; do
  ln -sfr ${ANN_DIR}/$p/*gff3.gz input/annotations/gff/$p.gff3.gz
  ln -sfr ${ANN_DIR}/$p/*cds.all.fa.gz input/annotations/cds/$p.cds.fa.gz
  ln -sfr ${ANN_DIR}/$p/*pep.all.fa.gz input/annotations/pep/$p.pep.fa.gz
done < tmp/file_list

```

The fat sand rat genome is not formatted in a standard way:

```sh

cat \
  ../input_data/2018-06-26-input/sand_rat_gene_annotation/Sandrat_geneV2.1.gff \
  | bgzip -c > input/annotations/gff/psammomys_obesus.gff3.gz

cat \
  ../input_data/2018-06-26-input/sand_rat_gene_annotation/Sandrat_geneV2.1.cds.fa \
  | bgzip -c > input/annotations/cds/psammomys_obesus.cds.fa.gz

cat \
  ../input_data/2018-06-26-input/sand_rat_gene_annotation/Sandrat_geneV2.1.pep.fa \
  | bgzip -c > input/annotations/pep/psammomys_obesus.pep.fa.gz


```

## Get the input tree

```sh

ln -sfr ../01-orthogroups/input/2019-02-11-tree.newick input/tree.newick

```


## Input for whole sequence alignment

### Filter orthogroups

We use an Rscript to filter the orthogroups. We allow orthogroups where there has been a lineage-specific increase or decrease in copy number in a maximum of 2 species. We use a second script to print a fasta file for each outgroup and to make a list of all the outgroups.

```sh

# Decompress the fasta file

mkdir -p tmp/annotations/cds

ls input/annotations/cds* \
   | grep -v "psammomys_obesus" \
   | cut -f 5 -d "/" | cut -f 1 -d "." > tmp/species_not_psammomys_obesus

while read p; do
  zcat input/annotations/cds/$p.cds.fa.gz \
    | ruby -pe 'gsub(/\.[0-9]+.*/, "")' \
    > tmp/annotations/cds/${p}.fa
done < tmp/species_not_psammomys_obesus

zcat input/annotations/cds/psammomys_obesus.cds.fa.gz \
  | ruby -pe 'gsub(/ .*/, "")' \
  > tmp/annotations/cds/psammomys_obesus.fa


conda activate r

Rscript bin/filter_orthogroups.R 2> tmp/filter_orthogroups.R.err

Rscript bin/make_alignment_input.R 2> tmp/make_alignment_input.R.err

conda deactivate

```

### Align

```sh

wget https://bioweb.supagro.inra.fr/macse/releases/macse_v2.03.jar

conda activate biobase
mkdir -p tmp/alignments/raw
ls tmp/gene_sequences_per_orthogroup \
  | cut -f 1 -d '.' \
  | grep -v INDEX \
  > tmp/genes

cat tmp/genes \
  | parallel -j 55 java -jar macse_v2.03.jar \
  -prog alignSequences \
  -seq tmp/gene_sequences_per_orthogroup/{}.fa \
  -out_NT tmp/alignments/raw/{}_nuc.aln \
  -out_AA tmp/alignments/raw/{}_aa.aln \
  1> tmp/alignments/raw/alignment.out \
  2> tmp/alignments/raw/alignment.err

```

## Clean alignment

We use HmmCleaner.pl to remove 'nonhomologous' sequences from each MSA. A difficult is that the alignment output from MACSE includes non-standard characters that have first to be removed from the alignment.

```sh

# First download macse_v2.03.jar

mkdir -p tmp/alignments/clean

cat tmp/genes \
  | parallel -j 40 java -jar macse_v2.03.jar \
  -prog exportAlignment \
  -align tmp/alignments/raw/{}_nuc.aln \
  -codonForInternalStop NNN \
  -codonForExternalFS --- \
  -codonForInternalFS --- \
  -out_NT tmp/alignments/clean/{}_nuc.aln \
  -out_AA tmp/alignments/clean/{}_aa.aln \
  1> tmp/alignments/clean/clean.out \
  2> tmp/alignments/clean/clean.err &

```

We can then clean the alignment, and re-create the nucleotide sequence of the cleaned alignment.

```sh

conda deactivate
conda deactivate

# First install HmmCleaner.pl and add to PATH

# Clean
cat tmp/genes \
  | parallel -j 10 \
  'HmmCleaner.pl tmp/alignments/clean/{}_aa.aln' \
  1> tmp/alignments/clean/clean.out \
  2> tmp/alignments/clean/clean.err &

# Amino acid to nucleotide
cat tmp/genes \
  | parallel -j 40 java -jar macse_v2.03.jar \
  -prog reportMaskAA2NT \
  -align_AA tmp/alignments/clean/{}_aa_hmm.fasta \
  -align tmp/alignments/clean/{}_nuc.aln \
  -min_NT_to_keep_seq 30 \
  -mask_AA $ \
  -min_seq_to_keep_site 4 \
  -min_percent_NT_at_ends 0.3 \
  -dist_isolate_AA 3 \
  -min_homology_to_keep_seq 0.3 \
  -min_internal_homology_to_keep_seq 0.5 \
  -out_NT tmp/alignments/clean/{}_nuc_hmm.fasta \
  1> tmp/alignments/clean/nuc_hmm.err &

```

There were 166 alignments that failed:

```sh

sort tmp/genes > tmp/gene_sequences_per_orthogroup_gene

ls tmp/alignments/clean/*_nuc_hmm.fasta | cut -f 4 -d "/" \
  | cut -f 1 -d "_" | sort > tmp/successful_alignments

diff tmp/successful_alignments tmp/gene_sequences_per_orthogroup_gene \
  | grep ">" | wc -l
# 166

```

Write to results directory:

```sh

mkdir -p tmp/alignments/hmm_clean

cat tmp/genes \
  | parallel -j 5 mv tmp/alignments/clean/{}_nuc_hmm.fasta \
    tmp/alignments/hmm_clean/{}.fa

```

We have successfully aligned 11,220 genes.

### Post alignment filter

#### Alignment quality

We first measure, for each alignment, the length of the sequence and the total number of gaps for each species. The following script saves three matrices in `tmp/gene_sizes`, and a number of figures:

```sh

conda activate r

Rscript bin/alignments2length.R

conda deactivate

cat tmp/gene_sizes/gap_percentages \
  | ruby -pe 'sub(/_nuc_hmm/, "")' \
  > tmp/gene_sizes/gap_percentages_edited

cat tmp/gene_sizes/gap_sizes \
  | ruby -pe 'sub(/_nuc_hmm/, "")' \
  > tmp/gene_sizes/gap_sizes_edited

cat tmp/gene_sizes/gene_sizes \
  | ruby -pe 'sub(/_nuc_hmm/, "")' \
  > tmp/gene_sizes/gene_sizes_edited

```

We then need to apply a filter to the alignments, removing sequences with too many gaps.

```sh

conda activate r
Rscript bin/filter_alignments.R
conda deactivate

```

The last step creates a directory, `results/masked_alignments`, which stores the masked alignments, and an index `results/masked_alignments_index`.

## Run bppml and mapnh

Because some species are missing a gene from the alignment files, we trim the tree for each alignment to include only the species present in the alignment.

```sh

mkdir -p tmp/trimmed_trees
conda activate r-pegas
Rscript bin/trim_tree.R \
  --index results/masked_alignments_index \
  --tree input/tree.newick \
  --out_dir tmp/trimmed_trees \
  2> tmp/trimmed_trees/tree.err
conda deactivate

```

We want to run `bppml` to estimate the branch length of each alignment. The programme requires a manual input of how many species there are in the alignment, so that's our first step:

```sh

mkdir -p tmp/biopp

parallel 'grep -c ">" results/masked_alignments/*.fa \
  | grep ":{}" \
  | cut -f 1 -d ":" \
  | cut -f 3 -d "/" \
  | cut -f 1 -d "." \
  > tmp/biopp/gene{} ' ::: {2..15}

#  cat tmp/biopp/gene* | wc -l
# 10554
```

After classing each alignment by number of species, we run `bppml` to estimate the branch length of each alignment:

```sh

# Download and install bppSuite
# cp -r ~/bppSuite .

mkdir -p tmp/biopp/ml

for i in {2..15}; do

  parallel -j 60 './bppSuite/bppml \
      GENE={2} \
      NODE_NUMBER={1} \
      INPUT_ALIGNMENT=results/masked_alignments/{2}.fa \
      INPUT_TREE=tmp/trimmed_trees/{2}.newick \
      OUTPUT_DIR=tmp/biopp/ml \
      param=bin/YN98.bpp \
      1> tmp/biopp/ml/{2}.out \
      2> tmp/biopp/ml/{2}.err' ::: `echo "(2 * ${i}) - 3" | bc` :::: tmp/biopp/gene${i}

done

ls tmp/biopp/ml/*dnd | wc -l
# 10014, as expected

```

We can now run `mapnh` on the estimated branch lengths to estimate dS and dN for each mutational category.

```sh

mkdir -p tmp/biopp/dnds

scp -r \
 susu:~/analyses/gerbil_gc/analysis/2019-02-02-high_gc_region_alignment/TestNH .

cat tmp/biopp/gene* | parallel -j 60 './TestNH/mapnh \
  GENE={} \
  INPUT_ALIGNMENT=results/masked_alignments/{}.fa \
  ML_DIR=tmp/biopp/ml/ \
  OUT_DIR=tmp/biopp/dnds/ \
  param=bin/mapnh_dnds.bpp \
  1> tmp/biopp/dnds/{}.out \
  2> tmp/biopp/dnds/{}.err'

conda deactivate

ls tmp/biopp/dnds/*_dS_X_W\-\>S.dnd \
   | cut -f 1 -d '.' \
   | cut -f 4 -d '/' \
   | sort > tmp/biopp/dnds_worked
cat tmp/biopp/gene* | sort > tmp/biopp/all_genes

diff tmp/biopp/dnds_worked tmp/biopp/all_genes

mkdir results/dnds
while read p; do
  cp tmp/biopp/dnds/${p}.counts_* results/dnds
done < tmp/biopp/dnds_worked

```

```r

dnds_worked <- read.table("tmp/biopp/dnds_worked")
dnds_worked <- as.character(dnds_worked$V1)

alignments_index <- read.table("results/masked_alignments_index", header = TRUE)

successful_dnds_index <- alignments_index[match(dnds_worked, alignments_index$orthogroup), ]

write.table(successful_dnds_index,
            file = 'results/successful_dnds_index',
            row.names = FALSE,
            quote = FALSE)

```

### Retrieving dN and dS per lineage

We then need to retrieve the branch lengths of each tree. This gets stored in the table `results/muridae_dnds_per_orthogroup.txt`.

```sh

conda activate r-pegas

Rscript bin/muridae_tree_to_branch_length.R \
  --index results/successful_dnds_index \
  --in_dir results/dnds \
  --out_table results/muridae_dnds_per_orthogroup.txt

conda deactivate

```

In the following script, we get the output from the last step and turn it into a wide format.


```sh

conda activate r

Rscript bin/dn_ds_2_long.R \
  --input results/muridae_dnds_per_orthogroup.txt \
  --output results/muridae_dnds_per_orthogroup_wide.txt

conda deactivate

```

## Give mouse coordinates


For each orthogroup in the final set, we need to get the reference gene.

```sh

ln -sfr ../01-orthogroups/results/orthofinder/Orthogroups.csv \
  input/original_orthogroups_index

zcat input/annotations/gff/mus_musculus.gff3.gz \
  | awk '$3 == "mRNA" { print }' > tmp/mouse_mrna.gff

conda activate r
Rscript bin/parse_gff.R -a tmp/mouse_mrna.gff -o tmp/mouse_tx_locations.csv
Rscript bin/add_gene_aa.R \
  -a tmp/mouse_tx_locations.csv \
  -i input/longest_tx_ids \
  -o tmp/mouse_gene_tx_aa_locations.csv

Rscript bin/add_orthogroup.R \
  -a tmp/mouse_gene_tx_aa_locations.csv \
  -m input/original_orthogroups_index \
  -o tmp/mouse_reference.csv

## Give mouse coordinates
Rscript bin/dn_ds_measurements_mouse_coord.R \
  --reference_csv tmp/mouse_reference.csv \
  --rates results/muridae_dnds_per_orthogroup_wide.txt \
  --out_csv results/2020-03-03-gerbil_measurements_mouse_coord.csv

cp tmp/mouse_reference.csv results/orthogroup2gene.csv

```

## Tidy up the results

```sh

tar -czf results/masked_alignments.tar.gz results/masked_alignments

# rm -rf results/masked_alignments

tar -czf results/dnds.tar.gz results/dnds
# rm -rf results/dnds

mv tmp/trimmed_trees results
ln -sfr results/trimmed_trees tmp

```

## Branch test in PAML

## Run branch test with codeml for each gene in the alignment

```sh

mkdir -p results/paml

tail -n +2 results/masked_alignments_index \
  | cut -f 1 \
  | parallel -j 40 'sh bin/run_codeml.sh {} \
  results/masked_alignments \
  tmp/trimmed_trees \
  results/paml'

# Some genes failed:
mkdir -p results/paml_failed_genes
grep "After deleting gaps. 0 sites" results/paml/* \
 | cut -f 3 -d "/" \
 | cut -f 1 -d "_" \
 | sort \
 | uniq \
  > results/paml_failed_genes/failed_genes

while read p; do
  mv results/paml/${p}_branch.txt results/paml_failed_genes
  mv results/paml/${p}_null.txt results/paml_failed_genes
done < results/paml_failed_genes/failed_genes

ls results/paml/*branch* | wc -l

python bin/parse_codeml.py results/paml \
  1> results/likelihoods.txt \
  2> results/likelihoods.txt.err

wc -l results/likelihoods.txt
#  21706 results/likelihoods.txt
wc -l results/likelihoods.txt.err
# 0 results/likelihoods.txt.err

```

Run a likelihood-ratio test for each alignment.

```sh

conda activate r
Rscript bin/paml-LRT.r

```

### Give gene name to branch test results


```sh

Rscript bin/PAML_LRT_mouse_coord.R \
  --reference_csv tmp/mouse_reference.csv \
  --rates results/log_ratio_test.csv \
  --out_csv results/2020-03-03-log_ratio_test-mouse_coord.csv

```

### Parse branch lengths for null model

We first get the branch length table for each gene.

```sh

grep -v -f results/paml_failed_genes/failed_genes results/masked_alignments_index \
 > results/paml_masked_alignments_index

```


```r

gene              <- read.csv("results/orthogroup2gene.csv")
genes_in_analysis <- read.table("results/paml_masked_alignments_index", header = TRUE)

gene <- gene[match(genes_in_analysis$orthogroup, gene$orthogroup), ]

write.csv(gene, "tmp/paml_successful_gene_orthogroup.csv")

```

```sh
mkdir -p tmp/branch_lengths

ls results/paml/*null.txt \
  | cut -f 3 -d "/" | cut -f 1 -d "_" \
  | parallel 'grep -A 50 "dN & dS for each branch" results/paml/{}_null.txt \
  | grep "\.\." \
  | cut -f 3- -d " " \
  > tmp/branch_lengths/{}'

conda activate r-pegas                                                                     

Rscript bin/null_branch_lengths.R

```

## Branch-site test

We use the branch-site test to identify sand rat proteins that show evidence of positive selection at specific residues. We use an implementation of this test that incorporates codon substitution rate variation and thereby controls for variation in the synonymous substitution rate caused by factors such as GC-biased gene conversion ([Davydov, Salamin and Robinson-Rechavi, 2019](https://academic.oup.com/mbe/article/36/6/1316/5371074)). This model is implemented in the programme Godon.

```sh


wget https://bitbucket.org/Davydov/godon/downloads/godon-master-linux-gnu-x86_64 -O godon

chmod +x godon
# branch: master, revision: ff28c19a52864162342c0756577eb656999439f4, build time: 2020-06-13_15:14:44

mkdir -p results/godon_branch_site
tail -n +2 results/masked_alignments_index \
  | cut -f 1 \
  | parallel -j 6 'bash bin/run_godon_branch_site.sh {} \
  results/masked_alignments \
  tmp/trimmed_trees \
  results/godon_branch_site'

```


Then we parse the output (D), defined as the difference in likelihoods between the null model (no positive selection at any residue) and the alternative model (positive selection at some residues).

```sh

ls results/godon_branch_site | cut -f 1 -d "_" > tmp/godon_list

while read p; do
  echo $p >> tmp/godon_branch_site_d
  grep "Final D" results/godon_branch_site/${p}_bs_codon_gamma.txt>> tmp/godon_branch_site_d
done < tmp/godon_list

cat tmp/godon_branch_site_d | ruby -pe 'gsub(/Final D=/, "")' \
  | paste - - >  tmp/godon_branch_site_d_parsed

```

We use a likelihood ratio test (LTR) to calculate the P-value for each branch-site test.


```r

# Performs likelihood-ratio test on Godon output,
  # comparing branch and null models

# ---------------------------------------------------------------------------- #

tests           <- read.table("tmp/godon_branch_site_d_parsed")
colnames(tests) <- c("id", "D")

# ---------------------------------------------------------------------------- #
# Degrees of freedom
tests$df <- 1

# P-value
tests$p.val <- pchisq(tests$D, tests$df, lower.tail=FALSE)

```

Note that multiple testing correction is not straightforward in this data set because there is a huge excess of tests where the P-value is 1. We applied the 1-sided robust FDR correction method ([Pounds and Cheng, 2006](https://academic.oup.com/bioinformatics/article/22/16/1979/208386)) for tests where `P < 0.99`:


```r

fdr_correction <- function(p) {

  require(prot2D)

  # Only run where p < 0.99
  p[p < 0.99] <- robust.fdr(p[p < 0.99])$q

  return(p)

}

tests$q_fdr <- fdr_correction(tests$p.val)

```
