#!/usr/bin/env python3
import os
import sys
import gffutils
from argparse import ArgumentParser
from collections import defaultdict

def get_longest_transcripts_in_gene(db):
    """
    return a dictionary where the keys are gene_ids and the values are the
    transcript with longest CDS
    """
    genes = {}
    for gene in db.features_of_type('gene'):
        if 'gene_id' in gene.attributes.keys():
            gene_id        = gene['gene_id'][0]
        elif 'gene' in gene.attributes.keys():
            gene_id        = gene['gene'][0]
        ## Make sure that biotype=protein_coding
        if 'biotype' in gene.attributes.keys():
            if gene['biotype'][0] == 'protein_coding':
                if (len(list(db.children(gene, featuretype='mRNA')))) > 0:
                    genes[gene_id] = gff_features(db, gene)
        elif 'gene_biotype' in gene.attributes.keys():
            if gene['gene_biotype'][0] == 'protein_coding':
                if (len(list(db.children(gene, featuretype='mRNA')))) > 0:
                    genes[gene_id] = gff_features(db, gene)
    return genes

def gff_features(db, gene):
    # Measure the length of each transcript
    ## Make sure you only get largest transcript for each gene
    tx_width = 0
    tx_width_max = 0
    tx_list = db.children(gene, featuretype='mRNA')
    for tx in tx_list:
        tx_width = db.children_bp(tx, child_featuretype='CDS')
        if 'Name' in tx.attributes.keys():
            gene_name = tx['Name'][0]
        else:
            gene_name = ""
        if tx_width > tx_width_max:
            tx_width_max  = tx_width
            max_trancript = tx
            max_gene      = gene_name
            # Get protein ID from the first CDS entry
            max_cds = list(db.children(tx, featuretype='CDS'))[0]['protein_id'][0]
    if tx_width_max > 0:
        tx_id          = max_trancript['transcript_id'][0]
        return [tx_id, max_cds, max_gene]
    # Print the gene, longest transcript mRNA, and all children
    # print(gene)
    # print(max_trancript)
    # for all_children in db.children(max_trancript, order_by="start"):
    #     print(all_children)

def get_gff_db(gff, db_file=":memory:"):
    """
    create a gffutils DB from a GFF file and will use an existing gffutils
    database if it already exists.
    """
    if db_file == ":memory:":
        db = gffutils.create_db(gff, dbfn=db_file, merge_strategy='create_unique')
    else:
        if not os.path.exists(db_file):
            db = gffutils.create_db(gff, dbfn=db_file, merge_strategy='create_unique')
        else:
            print("Reading db from", db_file, file=sys.stderr)
        db = gffutils.FeatureDB(db_file, keep_order=True)
    return db

if __name__ == "__main__":
    description = ("Get the longest protein-coding gene, which should be tagged in GFF as biotype=protein_coding")
    parser = ArgumentParser(description)
    parser.add_argument("gff", type=str, help="GFF file to calculate")
    parser.add_argument("--db_file", type=str, help="Optional tmp file in which to save DB")

    args = parser.parse_args()

    if args.db_file:
        db = get_gff_db(gff = args.gff, db_file = args.db_file)
    else:
        db = get_gff_db(gff = args.gff)

    gene_dict = get_longest_transcripts_in_gene(db)
    for gene, transcript_protein in gene_dict.items():
        print(gene, transcript_protein[0], transcript_protein[1], transcript_protein[2])
#
# db = get_gff_db(gff = 'test.gff')
#
# genes = {}
#
# for gene in db.features_of_type('gene'):
#     ## Make sure that biotype=protein_coding
#     if gene['biotype'][0] == 'protein_coding':
#         # Measure the length of each transcript
#         ## Make sure you only get largest transcript for each gene
#         tx_width = 0
#         tx_width_max = 0
#         tx_list = db.children(gene, featuretype='mRNA')
#
# for tx in tx_list:
#     tx_width = db.children_bp(tx, child_featuretype='CDS')
#     if tx_width > tx_width_max:
#         tx_width_max = tx_width
#         max_trancript = tx
#         genes[gene['gene_id'][0]] = max_trancript['transcript_id'][0]
