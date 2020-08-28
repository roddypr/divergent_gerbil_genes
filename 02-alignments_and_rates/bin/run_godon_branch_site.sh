#!/bin/bash

. ${HOME}/analyses/gerbil_gc/software/anaconda2/etc/profile.d/conda.sh

conda activate biobase

# ---------------------------------------------------------------------------- #
# Setup

# The start of the ctl file looks like this:
#  seqfile = in.fa      * sequence data filename
# treefile = in.newick  * tree structure file name
#  outfile = output.txt * main result file name

wd=$PWD
gene=$1
alignment_dir=$2
tree_dir=$3

# Output directory
results_dir=$4

gene_dir=tmp/branch_site_godon/$gene

echo "Running godon branch-site model for ${gene}"

# ---------------------------------------------------------------------------- #
# Run null model

## Make directory with gene name
mkdir -p $gene_dir

## Copy in gene sequence
cp ${alignment_dir}/${gene}.fa ${gene_dir}/in.fa

## Copy in tree
cat ${tree_dir}/${gene}.newick \
  | ruby -pe 'gsub(/\:(1|0\.5)/, "")' \
  | ruby -pe 'gsub(/\:[0-9]\.*[0-9]*/, "")' \
  | ruby -pe 'gsub(/psammomys_obesus/, "psammomys_obesus#1")' \
  | ruby -pe 'gsub(/meriones_unguiculatus/, "meriones_unguiculatus#1")' \
  | ruby -pe 'gsub(/\(psammomys_obesus\#1,meriones_unguiculatus\#1\)/,\
                   "(psammomys_obesus#1,meriones_unguiculatus#1)#1")' \
  > ${gene_dir}/in.newick

## run godon from directory
cd ${gene_dir}

# Optimise branch lengths with M0
$wd/godon M0 -p 4 --out-tree M0_tree.nwk in.fa in.newick > m0.out

# Run branch-site test, adding add the codon gamma rate variation
$wd/godon test BS -p 4 \
  --ncat-codon-rate 4 \
  in.fa \
  M0_tree.nwk \
  -j bs_codon_gamma.json \
  > bs_codon_gamma.out

cd $wd

mv ${gene_dir}/bs_codon_gamma.out ${results_dir}/${gene}_bs_codon_gamma.txt

# ---------------------------------------------------------------------------- #

echo "Finished running Godon branch site model for ${gene}"

# ---------------------------------------------------------------------------- #
# DONE
# ---------------------------------------------------------------------------- #
