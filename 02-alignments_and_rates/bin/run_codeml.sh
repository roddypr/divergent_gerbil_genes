#!/bin/sh

. ${HOME}/analyses/gerbil_gc/software/anaconda2/etc/profile.d/conda.sh

conda activate paml

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

directory_null=tmp/run_codeml/null/$gene
directory_branch=tmp/run_codeml/branch/$gene

echo "Running PAML for ${gene}"

# ---------------------------------------------------------------------------- #
# Run null model

## Make directory with gene name
mkdir -p $directory_null

## Copy in gene sequence
cp ${alignment_dir}/${gene}.fa ${directory_null}/in.fa

## Copy in tree
cat ${tree_dir}/${gene}.newick \
  | ruby -pe 'gsub(/\:(1|0\.5)/, "")' \
  | ruby -pe 'gsub(/\:[0-9]\.*[0-9]*/, "")' \
  | ruby -pe 'gsub(/psammomys_obesus/, "psammomys_obesus#1")' \
  | ruby -pe 'gsub(/meriones_unguiculatus/, "meriones_unguiculatus#1")' \
  | ruby -pe 'gsub(/\(psammomys_obesus\#1,meriones_unguiculatus\#1\)/,\
                   "(psammomys_obesus#1,meriones_unguiculatus#1)#1")' \
  > ${directory_null}/in.newick

## run codeml from directory
cd ${directory_null}

ln -sf ${wd}/bin/codeml-for-null-model.ctl .

codeml codeml-for-null-model.ctl \
  1> codeml-for-null-model.out \
  2> codeml-for-null-model.err

cd ${wd}

mv ${directory_null}/output.txt ${results_dir}/${gene}_null.txt

# ---------------------------------------------------------------------------- #
# Run alternative model

## Make directory with gene name
mkdir -p $directory_branch

## Copy in gene sequence
ln -srf ${directory_null}/in.fa ${directory_branch}/in.fa

## Copy in tree
ln -srf ${directory_null}/in.newick ${directory_branch}/in.newick

## run codeml from directory
cd ${directory_branch}

ln -sf ${wd}/bin/codeml-for-branch-model.ctl .

codeml codeml-for-branch-model.ctl \
  1> codeml-for-branch-model.out \
  2> codeml-for-branch-model.err

cd ${wd}

mv ${directory_branch}/output.txt ${results_dir}/${gene}_branch.txt

# ---------------------------------------------------------------------------- #
# DONE
# ---------------------------------------------------------------------------- #
