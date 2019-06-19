#!/usr/bin/env nextflow

genome_dir = Channel
      .fromPath(params.genome_dir)
      .ifEmpty { exit 1, "${params.genome_dir} not found"}
reads_dir = Channel
      .fromPath(params.reads_dir)
      .ifEmpty { exit 1, "${params.reads_dir} not found"}
conditions = Channel
      .fromPath(params.conditions)
      .ifEmpty { exit 1, "${params.conditions} not found"}

/*
 * STEP 1 - Split chrs, run REAL on all reads & run isosegmenter on all chromosomes
 */
process input {
    // cpus threads
    tag "$reads_dir"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file(genome_dir) from genome_dir
    file(reads_dir) from reads_dir
    file(conditions) from conditions

    output:
    file("*") into reads

    script:
    """
    input.py --genome_dir $genome_dir --reads_dir $reads_dir --conditions $conditions
    """
}