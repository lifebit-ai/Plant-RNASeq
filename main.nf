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

// Configurable variables
params.name = false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false

multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")

// Validate inputs
fasta = Channel
      .fromPath(params.fasta)
      .ifEmpty { exit 1, "${params.fasta} not found"}
      //.map { file -> tuple(file.baseName, file) }
      .into { fasta_real; fasta_split_chr }

// AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Check workDir/outdir paths to be S3 buckets if running on AWSBatch
// related: https://github.com/nextflow-io/nextflow/issues/813
if( workflow.profile == 'awsbatch') {
    if(!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}

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
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

nf-core/plant-rnaseq v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'nf-core/plant-rnaseq'
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
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


def create_workflow_summary(summary) {

    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-plant-rnaseq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/plant-rnaseq Workflow Summary'
    section_href: 'https://github.com/nf-core/plant-rnaseq'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}


/*
 * Parse software version numbers
 */
// process get_software_versions {
//
//     output:
//     file 'software_versions_mqc.yaml' into software_versions_yaml
//
//     script:
//     """
//     echo $workflow.manifest.version > v_pipeline.txt
//     echo $workflow.nextflow.version > v_nextflow.txt
//     fastqc --version > v_fastqc.txt
//     multiqc --version > v_multiqc.txt
//     scrape_software_versions.py > software_versions_mqc.yaml
//     """
// }



/*
 * STEP 1 - FastQC
 */
// process fastqc {
//     tag "$reads"
//     publishDir "${params.outdir}/fastqc", mode: 'copy',
//         saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}
//
//     input:
//     set val(name), file(reads) from read_files_fastqc
//
//     output:
//     file "*_fastqc.{zip,html}" into fastqc_results
//
//     script:
//     """
//     fastqc -q $reads
//     """
// }



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



/*
 * STEP 2 - MultiQC
 */
// process multiqc {
//     publishDir "${params.outdir}/MultiQC", mode: 'copy'
//
//     input:
//     file multiqc_config
//     file ('fastqc/*') from fastqc_results.collect()
//     file ('software_versions/*') from software_versions_yaml
//     file workflow_summary from create_workflow_summary(summary)
//
//     output:
//     file "*multiqc_report.html" into multiqc_report
//     file "*_data"
//
//     script:
//     rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
//     rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
//     """
//     multiqc -f $rtitle $rfilename --config $multiqc_config .
//     """
// }



/*
 * STEP 3 - Output Description HTML
 */
// process output_documentation {
//     tag "$prefix"
//     publishDir "${params.outdir}/Documentation", mode: 'copy'
//
//     input:
//     file output_docs
//
//     output:
//     file "results_description.html"
//
//     script:
//     """
//     markdown_to_html.r $output_docs results_description.html
//     """
// }



/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/plant-rnaseq] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/plant-rnaseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/plant-rnaseq] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/plant-rnaseq] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[nf-core/plant-rnaseq] Pipeline Complete"

}
