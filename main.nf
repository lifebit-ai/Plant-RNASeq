#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/plant-rnaseq
========================================================================================
 nf-core/plant-rnaseq Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/plant-rnaseq
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     nf-core/plant-rnaseq v${workflow.manifest.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run lifebit-ai/Plant-RNASeq --reads_folder /path/to/myFolder --fasta /path/to/my.fasta --conditions "4 6,9 10,13 14,17 18,21"

    Mandatory arguments:
      --reads_folder                Path to folder containing input FASTQ files
      --fasta                       Path to reference FASTA file used to align the reads against with real
      --conditions                  Integer to represent the number of conditions that the RNA data was grouped into intially eg 4

    Other options:
      --reads_prefix                Can be used to run the pipeline for a subset of the reads in the reads_folder
      --reads_extension             Extension of the FASTQ files in reads_folder (default = "fastq")
      --outdir                      The output directory where the results will be saved

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}
// Get as many processors as machine has
int threads = Runtime.getRuntime().availableProcessors()

// Validate inputs
fasta = Channel
      .fromPath(params.fasta)
      .ifEmpty { exit 1, "${params.fasta} not found"}
      //.map { file -> tuple(file.baseName, file) }
      .into { fasta_real; fasta_split_chr }

/*
 * Create a channel for input read files
 */
reads="${params.reads_folder}/*.${params.reads_extension}"
Channel
            .fromPath(reads)
            .map { file -> tuple(file.baseName, file) }
            .ifEmpty { exit 1, "${reads} was empty - no input files supplied" }
            .combine(fasta_real)
            .into { read_files_fastqc; read_files_real }

read_files_fastqc.subscribe{println "value: $it"}


// Header log info
log.info """=======================================================
lifebit-ai/plant-rnaseq v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'lifebit-ai/plant-rnaseq'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = reads
if(params.readPaths) summary['Reads'] = params.readPaths
if(params.reads) summary['Reads'] = params.reads
summary['Fasta Ref']    = params.fasta
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


/*
 * STEP 2 - REAL - align RNA data
 */
process real {
    cpus threads
    tag "$reads"
    publishDir "${params.outdir}/real", mode: 'copy'

    input:
    set val(name), file(reads), file(fasta) from read_files_real

    output:
    set val(name), file("*.aln") into real_output

    script:
    """
    real -T ${task.cpus} -p $reads -t $fasta -o ${name}.out
    awk '{ print >> \$9".aln" }' ${name}.out
    """
}


// emit individual SAM files with their prefix for no_reads process
real_output
   .flatMap { aln ->
        list_chrs = []
        aln[1].each { chrFile ->
                list_chrs << tuple(chrFile.baseName, aln[0], chrFile)
        }
        list_chrs
    }
    .set{aligned_reads_no_reads}


/*
 * STEP 3A - split_chr - split the fasta genome into seperate files for each chromosome using pyfaidx
 */
process split_chr {
   tag "$fasta"
   publishDir "${params.outdir}/split_chr", mode: 'copy'

   input:
   file fasta from fasta_split_chr

   output:
   file "*.fa" into chrs

   script:
   """
   faidx -x $fasta
   """
}


// emit individual FASTA files with their prefix for isoSegmenter
chrs.flatten()
    .map{ file -> tuple(file.baseName, file) }
    .set{chr}


/*
 * STEP 3B - isoSegmenter - segment genomes into isochores i.e. split the genome into chromosomes
 */
 process isosegmenter {

     tag "$chr"

     publishDir "${params.outdir}/isosegmenter", mode: 'copy'

     container 'bunop/isosegmenter:latest'

     input:
     set val(name), file(chr) from chr

     output:
     file "${name}.csv" into iso, iso_mk_gene_exp_input

     script:
     """
     isoSegmenter.py --infile $chr --window_size 100000 --outfile ${name}_wgap.csv
     awk 'BEGIN{FS=","}\$4!="gap"{print \$0}' ${name}_wgap.csv > ${name}.csv
     """
 }


// emit individual csv files with their prefix for no_reads
iso.flatten()
    .map{ file -> tuple(file.baseName, file) }
    .cross(aligned_reads_no_reads)
    .map { it ->
       [it[1][0], it[1][1], it[1][2], it[0][1]]
     }
    .set{ iso_no_reads }


/*
 * STEP 4A - no_reads - compute number of reads found in each isochore
 */
process no_reads {

    publishDir "${params.outdir}/no_reads", mode: 'copy'

    input:
    set val(isochore_name), val(sample_name), file(aligned_read), file(isochore) from iso_no_reads

    output:
    file "*.csv" into csv

    script:
    """
    sort -k 10,10n $aligned_read > ${aligned_read}.s
    noReads.py $isochore ${aligned_read}.s > ${isochore_name}_${sample_name}.csv
    """
}



/*
 * STEP 4B - combine outputs into one file in the correct format
 */
process mk_gene_exp_input {

    publishDir "${params.outdir}/mk_gene_exp_input", mode: 'copy'

    input:
    file iso from iso_mk_gene_exp_input.collect()
    file csv from csv.collect()

    output:
    file "gene_exp_input.csv" into gene_exp

    script:
    """
    mk_gene_exp_input.py
    """
}



/*
 * STEP 5 - gene_exp - compute the gene expression
 */
process gene_exp {
    tag "output_gene_expression.csv"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    file isochore_reads from gene_exp

    output:
    file "output_gene_expression.csv" into results

    script:
    """
    geneExp.py $isochore_reads $params.conditions >> output_gene_expression.csv
    """
}