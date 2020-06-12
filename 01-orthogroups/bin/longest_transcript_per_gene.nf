#!/usr/bin/env nextflow

/*
  * Select one transcript for each protein-coding gene in a given number of
  * species
*/

// Inputs
params.gff    = "input/annotations/gff/*.gff3.gz"
params.pep_fa = "input/annotations/pep/*.pep.fa.gz"

params.output       = "results/longest_tx_fa"
params.output_index = "results/longest_tx_fa/INDEX"

gff = Channel
              .fromPath(params.gff)
              .filter { it.simpleName != "psammomys_obesus" }
              .map { it -> tuple(it.simpleName, it) }


pep_fa = Channel
              .fromPath(params.pep_fa)
              .filter { it.simpleName != "psammomys_obesus" }
              .map { it -> tuple(it.simpleName, it) }

gff
    .combine(pep_fa, by: 0)
    .set { gff_pep }

// Select the longest transcript per protein-coding gene.
process chooseLongestTranscript {

    conda '/home/zoo/zool2304/analyses/gerbil_gc/software/anaconda2/envs/biopython'

    maxForks 25

    input:
    set val(id), file('gff.gz'), file('fa.gz') from gff_pep

    output:
    set val(id), file('longest_tx'), file('fa.gz') into longest_tx_channel
    set val(id), file('longest_tx') into index_channel

    """
    python ${baseDir}/longest_transcript_per_gene.py gff.gz \
        > longest_tx

    """
}

process addSpecies {

    conda '/home/zoo/zool2304/analyses/gerbil_gc/software/anaconda2/envs/biopython'

    maxForks 25

    input:
    set val(id), file('longest_tx') from index_channel

    output:
    file 'species_longest_tx' into index_species_channel

    """

    cat longest_tx | awk -v var=\"${id}\" '{print var, \$0}' \
      > species_longest_tx

    """
}

index_species_channel
    .collectFile(name: params.output_index, newLine: true)

process getSequenceOfLongestTx {

    conda '/home/zoo/zool2304/analyses/gerbil_gc/software/anaconda2/envs/biopython'

    publishDir params.output, mode: 'move'

    maxForks 25

    input:
    set val(id), file('longest_tx'), file('fa.gz') from  longest_tx_channel

    output:
    file "${id}.pep.fa"

    """

    cut -f 3 ${longest_tx} -d " " > longest_tx.list

    # Parse fasta to remove spaces
    zcat fa.gz | ruby -pe 'gsub(/\\.\\d+ .+/, "")' > parsed_query.fa
    seqtk subseq -l 80 parsed_query.fa longest_tx.list > ${id}.pep.fa

    """
}


workflow.onComplete {
  println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
