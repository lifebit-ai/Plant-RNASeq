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
 * STEP 1 - Split chrs, run REAL on all reads
 */
process input {
    // cpus threads
    tag "$reads_dir"
    publishDir "${params.outdir}/input", mode: 'copy'

    input:
    file(genome_dir) from genome_dir
    file(reads_dir) from reads_dir
    file(conditions) from conditions

    output:
    file("aligned") into real_output
    file("chromosomes") into chrs

    script:
    //TODO: add extra params
    """
    input.py --genome_dir $genome_dir --reads_dir $reads_dir --conditions $conditions
    """
}

/*
 * STEP 2 - isoSegmenter - segment genomes into isochores i.e. split the genome into chromosomes
 */
process isoSegmenter {
    // cpus threads
    tag "$reads_dir"
    publishDir "${params.outdir}/input", mode: 'copy'

    input:
    file(chrs) from chrs

    output:
    file("*") into reads

    script:
    //TODO: add extra params
    """
    isosegmenter.py
    """
}

/*
 * STEP 3 - Compute the number of reads for each chromosome and isochore file
 */
// process reads {
//     // cpus threads
//     publishDir "${params.outdir}/reads", mode: 'copy'

//     input:
//     file(reads) from reads

//     output:
//     file("*") into expression

//     script:
//     """
//     reads.py
//     """
// }

// /*
//  * STEP 4 - compute the gene expression
//  */
// process expression {
//     // cpus threads
//     publishDir "${params.outdir}/expression", mode: 'copy'

//     input:
//     file(expression) from expression

//     output:
//     file("*") into results

//     script:
//     """
//     expression.py
//     """
// }