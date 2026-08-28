#!/usr/bin/env nextflow

// Single Cell RNA, CITE-Seq pipeline for hashed samples

// Parameters used in this NextFlow pipeline are defined in scCITE_NF_CVID.config file

def helpMessage(){
    ///log.info nfcoreHeader()
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run SingleCellNF_HPC.nf nextflow.config
    Mandatory arguments:
      --fastqdir                    Path to fastq data (not required for SRA download)
      --chemistry                   Version of 10x chemistry, e.g. "--chemistry v2" or "--chemistry v3" or "--chemistry 5prime"
      --aligner                     Name of the tool to use for scRNA (pseudo-) alignment. Available are: "none", "cellranger", "alevin" or "kallisto". Default 'cellranger'.
      --library_type                Either: EXP, CITE, HASH, CITE_HASH or ATAC
      --renamefastq                 true/false, Set true when renaming is needed, not needed when Bcl2toFastq was performed or when similar datastructure is provided
    SRA download Options:
      --sradownload                 true/false, Set true to download datasets from SRA
      --sratxt                      Text file containing SRA ids (ids are line separated in .txt file), e.g. "sra_ids.txt"
      (--sradir                      Directory to store SRA files in)
    Renamefastq arguments:
      --renamedfastqstore           Directory to store renamed (and zipped) fastq files
      --extension1                  Extension of SRA file _1, options: "I1" for index, "R1" for read 1 or "R2" for read 2
      --extension2                  Extension of SRA file _2, options: "I1" for index, "R1" for read 1 or "R2" for read 2
      --extension3                  Extension of SRA file _3, options: "I1" for index, "R1" for read 1 or "R2" for read 2
    Bcl2toFastq Options:
      --bcl2fastq                   true/false, Set true when conversion is needed
      --barcodeMM                   mismatches, set at 0 if indexes differ with only 1 basepair
      --bcl2dir                     bcl2 files directory. Mandatory when bcl2fastq is true.
      --bcl2gz                      name of tar.gz bcl2 files
      --samplesheet                 An illumina samplesheet containing the 10X sample information and indexes
      --sampleindexcsv              A csv file with "Lane,Sample,Index" for each sample
      --untarbcl2                   true/false, If the bcl2 dir is tar file, put true.
    CellRanger Options:
      --cellranger_index            Supply a cellranger index to pipeline when performing cellranger quantification
      --cellranger_vdjindex         Supply a cellranger vdj index to pipeline when performing cellranger VDJ quantification
      --feature_ref                 Needed when for CITE-seq and Hashing
    CellRanger Aggregate Options:
      --aggregate                   true/false, Set true when aggregation is needed
      --libraries                   Set when available, will be created otherwise. CSV file with fastqs, sample, library type
    Alevin Options:
      --alevin_gene_map             A gene map used for Alevin
      --alevin_index                Supply an Alevin index to pipeline when performing cellranger quantification
      --barcode_whitelist           Custom file of whitelisted barcodes
    Kallisto/BusTools Options:
      --kallisto_gene_map           A gene map used for bustools correction of BUS files created by Kallisto
      --kallisto_index              Supply a Kallisto index to pipeline when performing kallisto pseudomapping
      --barcode_whitelist           Custom file of whitelisted barcodes
    Popscle Options:
      --freemuxlet                  true/false, Set true when running freemuxlet is needed
      --vcfref                      Path to a vcf reference
      --nsamples                    Numeric value for number of samples per single cell pool
    QC Options:
      --qcsample                    true/false, Set true when running QCreports is needed
      --annotationfile              Path to annotation file, including hashtags when needed
      --rmarkdown                   Path to markdown file used for qc
    Velocyto Options:
      --gtf                         GTF file, same to cellranger quantification index
    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --name                        Name for the pipeline run
    """.stripIndent()
}

// PREPARATION OF PIPELINE 
// Show help message when needed
if (params.help){
    helpMessage()
    exit 0
}

// Check if necessary inputs are present
// Check the aligner (kallisto, alevin, cellranger)
if ( params.aligner != "none" && params.aligner != 'cellranger' && params.aligner != 'alevin' && params.aligner != 'kallisto' ){
    exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'none','cellranger','alevin', 'kallisto'"
}

// Check the library tyep
if ( params.library_type != "EXP" && params.library_type != 'CITE' && params.library_type != 'CITE_TCR_BCR' && params.library_type != 'HASH' && params.library_type != 'CITE_HASH' && params.library_type != 'ATAC'  ){
    exit 1, "Invalid library type option: ${params.library_type}. Valid options: 'EXP','CITE','HASH', 'CITE_HASH', 'CITE_TCR_BCR','ATAC'"
}

// Check if cellranger index is supplied properly and convert chemistry
if( params.aligner == 'cellranger' ){
  Channel
  .fromPath(params.cellranger_index)
  .ifEmpty { exit 1, "Cellranger index not found: ${params.cellranger_index}" }
  if ( params.chemistry == 'v3' ){
    cellrangerchemistry = "SC3Pv3"
    } else if ( params.chemistry == 'v2' ){
      cellrangerchemistry = "SC3Pv2"
      } else if ( params.chemistry == "5prime"){
        cellrangerchemistry = "SC5P-R2"
        } else {
        log.info"""
        No chemistry provided, cellranger defaulting to autodetect
        """
        cellrangerchemistry = "auto"
      }
}

// Check if alevin index is supplied properly
if( params.aligner == 'alevin' ){
    Channel
        .fromPath(params.alevin_index)
        .ifEmpty { exit 1, "Alevin index not found: ${params.alevin_index}" }
    Channel
        .fromPath(params.alevin_gene_map)
        .ifEmpty { exit 1, "Alevin gene map not found: ${params.alevin_gene_map}" }
    if (params.chemistry == 'v3'){
        alevinchemistry = "chromiumV3"
    } else if (params.chemistry == 'v2'){
        alevinchemistry = "chromium"
    } else {
        exit 1, "Invalid chemistry type option: ${params.chemistry}. Valid options: 'v2','v3'"
    }
}


// Check if kallisto index is supplied properly
if( params.aligner == 'kallisto' ){
    Channel
        .fromPath(params.kallisto_index)
        .ifEmpty { exit 1, "Kallisto index not found: ${params.kallisto_index}" }
}

// Check if barcode whitelist is supplied properly (not needed for cell ranger)
// barcode white list not integrated yet
if( params.aligner != "cellranger" && params.barcode_whitelist ){
    Channel
        .fromPath(params.barcode_whitelist)
        .ifEmpty { exit 1, "Barcode whitelist is not found: ${params.barcode_whitelist}" }
}

//Check if feature reference file is provided when performing cellranger hashing / CITE-seq
if((params.library_type == "CITE" | params.library_type == "HASH" | params.library_type =="CITE_HASH" | params.library_type =="CITE_TCR_BCR") && params.aligner == 'cellranger'){
    Channel
        .fromPath(params.feature_ref)
        .ifEmpty { exit 1, "Feature reference is not found: ${params.feature_ref}" }
}

//Check if vdj reference file is provided when performing cellranger hashing / CITE-seq
if((params.library_type == "CITE_TCR_BCR") && params.aligner == 'cellranger'){
    Channel
        .fromPath(params.cellranger_vdjindex)
        .ifEmpty { exit 1, "Feature reference is not found: ${params.cellranger_vdjindex}" }
}

//Check if SRA ids are provided when SRA download is desired
if(params.sradownload){
    Channel
      .fromPath(params.sratxt)
      .ifEmpty { exit 1, "SRA ids text file not found: ${params.sratxt}" }

}

if (params.renamefastq && params.bcl2fastq) {
    exit 1, "Invalid, you can't rename and convert at the same time, as the naming is already correct by converting."
}

if (params.sradownload && !params.renamefastq){
  exit 1, "For SRA download, fastq rename is required"
}


//Function to extract ID from path
def getID( filename ) {
    regexpPE = /(.+)\/(.+)/
    (filename =~ regexpPE)[0][2]
}

//Function to extract ID
def getrunname( filename ) {
    regexpPE = /(.+)_(A|B)(.+)/
    (filename =~ regexpPE)[0][3]
}

//Set automatic outputdir when not provided
if (params.scratchdir){
    scratchDir = "${params.scratchdir}"
}
else {
    exit 1, "A scratchdir is required."
}

if (params.outdir){
    outputDir = "${params.outdir}"
}
else {
    outputDir = "${params.scratchdir}/${params.name}"
}

aggrdir = outputDir+"/"+params.name


// Inforamtion pipeline
log.info """\
 CVID scRNASEQ10x P I P E L I N E by Tine D'hamers
 Supervised by Celine Everaert
 =================================================
 TOBI lab
 =================================================
 your run: ${params.name}
 quantification algorithm: ${params.aligner}
 library type : ${params.library_type} + ${params.chemistry}
 output directory : ${outputDir}
 freemuxlet: ${params.freemuxlet}
 velocyto: ${params.velocyto}
 """



// Initiate  all input variables for processes at empty channel (because even inputs for unused channels need to be initiated)
runname = Channel.empty()
dirs_rename = Channel.empty()
fastqfiles_dirs = Channel.empty()
sampleIDs = Channel.empty()


// Revalue process inputs according to requested 'mode'
if( params.renamefastq == true){
    //When not correctly named fastq files, rename eg. when downloaded from SRA
    if (params.sradownload == true){
      log.info """\
      SRA download mode (+ rename mode)
      """
      dirs_rename = Channel.empty()
      fastqfiles_dirs = Channel.empty()
      targzfile = Channel.empty()

      sampleIDs = Channel
        .fromPath(params.sratxt)
        .splitText()
        .map{it -> it.trim()} //needed to trim "\n", otherwise error

    } else {
      log.info """\
      rename mode
      """
      //dir of fastqfiles to rename
      dirs_rename = Channel
                      .fromPath("${params.fastqdir}/*", checkIfExists:true, type:'dir')
                      .map{dir -> tuple(getID(dir), dir)}
      fastqfiles_dirs = Channel.empty()
      targzfile = Channel.empty()
      sampleIDs = Channel.empty()
    }
    if (params.library_type != "EXP"){
      exit 1, "Rename only possible for expression data (EXP library type)"
    }

//file_extensions
    if (params.extension1 && params.extension2 && params.extension3 ){
    file_extension1 = "${params.extension1}"
    file_extension2 = "${params.extension2}"
    file_extension3 = "${params.extension3}"
    } else {
    exit 1, "(One of the) file extension(s) not provided in config file"
    }

} else if (params.bcl2fastq == true){
    // When blc2 fastq conversion is still needed
    log.info """\
    bcl2fastq mode
    """
    Channel
        .fromPath(params.bcl2dir)
        .ifEmpty { exit 1, "Bcl2 directory not found: ${params.bcl2dir}"}

    if (params.untarbcl2) {
      Channel
          .fromPath("${params.bcl2dir}/*.tar.gz")
          .ifEmpty { exit 1, "Bcl2 gzip not found: ${params.runname}.tar.gz"}
      targzfile = Channel
          .fromPath("${params.bcl2dir}/*.tar.gz")
    } else {
      targzfile = Channel.fromPath(params.bcl2dir)
    }
    if ( !params.samplesheet && !params.sampleindexcsv){
        exit 1, "Invalid sample sheet for bcl2fastq. Either provide a samplesheet or a sampleindexcsv"
    }
    if ( params.sampleindexcsv ){
        Channel
            .fromPath(params.sampleindexcsv)
            .ifEmpty { exit 1, "sampleindex csv file not found: ${params.sampleindexcsv}"}
    }
    if ( params.samplesheet ){
        Channel
            .fromPath(params.samplesheet)
            .ifEmpty { exit 1, "sample sheet file not found: ${params.samplesheet}"}
    }
    dirs_rename = Channel.empty()
    fastqfiles_dirs = Channel.empty()
    sampleIDs = Channel.empty()
} else if (params.bcl2fastq == false && params.renamefastq == false){
    // When directories are filled with correctly named fastq files
    log.info """\
    quantification mode
    """
    fastqfiles_dirs = Channel
        .fromPath("${params.fastqdir}/*", checkIfExists:true, type:'dir')
        .map{dir -> tuple(getID(dir), dir)}

    dirs_rename = Channel.empty()
    targzfile = Channel.empty()
    sampleIDs = Channel.empty()

} else {
    dirs_rename = Channel.empty()
    fastqfiles_dirs = Channel.empty()
    targzfile = Channel.empty()
    sampleIDs = Channel.empty()
}


// Store runname
if (params.bcl2runname){
    runname = getrunname(params.bcl2runname)
} else {
    runname = params.name
}

// Check necessary files for freemuxlet
if (params.freemuxlet == true) {
    Channel
        .fromPath(params.vcfref)
        .ifEmpty { exit 1, "vcf reference file not found: ${params.vcfref}"}
    if (params.nsamples <= 1) {
        exit 1, "Invalid number of samples, you can't demultiplex 1 or no samples."
    }
}

// Check necessary files for quality control of the samples
if (params.qcsample == true) {
    Channel
        .fromPath(params.annotationfile)
        .ifEmpty { exit 1, "annotation file not found: ${params.annotationfile}"}
    Channel
        .fromPath(params.rmarkdown)
        .ifEmpty { exit 1, "markdown file for qc analysis not found: ${params.annotationfile}"}
}

// Check necessary files for velocyto
if (params.velocyto == true) {
    Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "gtf file not found: ${params.vcfref}"}
}

// Start PROCESSES
// Run CellRanger bcl to fastq conversion when needed for cite-seq/hashed data
process cellranger_mkfastq {

    // Output directory
    publishDir "$params.fastqdir", mode: 'copy', overwrite: true

    // Custom label of process (in log file)
    tag "${params.name}"

    // Load necessary module
    module 'purge:CellRanger/7.0.0:bcl2fastq2/2.20.0-GCC-11.2.0'
    label 'mid_mem'
    time '48h'

    input:
    val runname
    
    when: params.bcl2fastq

    output:
    // Store the path of the output files in a new variable
    path "fastqfiles/${runname}/*", type:"dir" into fastq_afterbcl_dirs

    script:
    """
    # Create directory for fastq files
    mkdir -p fastqfiles
    # Cellranger mkfastq workflow parameters
    #       --run: path of Illumina BCL run folder (outputfiles)
    #       --id: name of folder created by mkfastq
    #       --csv: information on samples (lane, sample and index) needed for demultiplexing
    #       --output-dir: generate fastq files in this directory
    #       --barcode_mismatches: number of mismatches allowed per index adapter
    #       --localcores: maximum cores requested for the process
    #       --localmem: maximum GB requested for this process

    cellranger mkfastq --id=mkfastq_${runname} \\
                    --run=${params.bcl2dir} \\
                    --csv=${params.sampleindexcsv} \\
                    --output-dir=fastqfiles \\
                    --barcode-mismatches=${params.barcodeMM} \\
                    --localcores=${task.cpus} \\
                    --localmem=70
    cd fastqfiles/${runname}
    for fastqfile in *.fastq.*; do
        filename=\${fastqfile%.fastq*}
        samplename=\${filename%%_*}
        mkdir -p "\$samplename"
        mv -i "\$fastqfile" "\$samplename"/
    done
    """
}

// Store the directory where fastq files are stored
if (params.bcl2fastq){
    fastq_afterbcl_dirs
    .flatten()
    .map{dir -> tuple(getID(dir), dir)}
    .into{fastqfiles_combined_cellranger; fastqfiles_combined_alevin; fastqfiles_combined_kallistobus; fastqfiles_combined_cellranger_atac}
}
else if (params.renamefastq){
    fastq_afterrename_dirs
    .flatten()
    .map{dir -> tuple(getID(dir), dir)}
    .into{fastqfiles_combined_cellranger; fastqfiles_combined_alevin; fastqfiles_combined_kallistobus; fastqfiles_combined_cellranger_atac}
}
else {
    fastqfiles_dirs.into{fastqfiles_combined_cellranger; fastqfiles_combined_alevin; fastqfiles_combined_kallistobus; fastqfiles_combined_cellranger_atac}
}

// Cell ranger count
// Generate count matrices and web summary reports 
process cellranger_count{

    tag "${samplename}"

    publishDir "$outputDir/$samplename/cellranger", mode: 'copy', overwrite: true

    module 'CellRanger/7.0.0'
    label 'big_mem'
    time '72h'
    scratch "$scratchDir"

    input:
    tuple val(samplename), path(fastqdir) from fastqfiles_combined_cellranger

    output:
    tuple val(samplename), path("${params.name}_${samplename}/outs", type:"dir") into countfiles

    when: params.aligner == "cellranger" && params.library_type != "ATAC"

    script:
    if (params.library_type == "CITE_HASH" )
        """
        # Create multi config CSV file necessary as input for cellranger_multi
            echo '[gene-expression]' > libraries.csv
            echo 'reference,${params.cellranger_index}' >> libraries.csv
            echo 'chemistry,${cellrangerchemistry}' >> libraries.csv
#           echo 'cmo-set,${params.cmo_ref}' >> libraries.csv
            echo '' >> libraries.csv
            echo '[feature]' >> libraries.csv
            echo 'reference,${params.feature_ref}' >> libraries.csv
            echo '' >> libraries.csv
            echo '[libraries]' >> libraries.csv
            echo 'fastq_id,fastqs,lanes,feature_types,subsample_rate' >> libraries.csv
            echo ${samplename}_cDNA,\$(pwd)'/${samplename},any,Gene Expression' >> libraries.csv
            echo ${samplename}_FB5P,\$(pwd)'/${samplename},any,Antibody Capture' >> libraries.csv
#           echo ${samplename}_FB5P,\$(pwd)'/${samplename},any,Multiplexing Capture' >> libraries.csv
#           echo ${samplename}_ADT,\$(pwd)'/${samplename},any,Antibody Capture' >> libraries.csv
#           echo ${samplename}_HTO,\$(pwd)'/${samplename},any,Antibody Capture' >> libraries.csv
#           echo '' >> libraries.csv
#           echo '[samples]' >> libraries.csv
#           echo 'sample_id,cmo_ids' >> libraries.csv
#           echo 'sample1,Hashtag1' >> libraries.csv
#           echo 'sample2,Hashtag2' >> libraries.csv
            cellranger multi \\
                --id=${params.name}_${samplename} \\
                --csv=libraries.csv \\
                --localcores=${task.cpus} \\
                --localmem=80 \\
        """
    else if (params.library_type == "CITE_TCR_BCR" )
            """
            echo '[gene-expression]' > libraries.csv
            echo 'reference,${params.cellranger_index}' >> libraries.csv
            echo 'chemistry,${cellrangerchemistry}' >> libraries.csv
            echo '' >> libraries.csv
            echo '[vdj]' >> libraries.csv
            echo 'reference,${params.cellranger_vdjindex}' >> libraries.csv
            echo '' >> libraries.csv
            echo '[feature]' >> libraries.csv
            echo 'reference,${params.feature_ref}' >> libraries.csv
            echo '' >> libraries.csv
            echo '[libraries]' >> libraries.csv
            echo 'fastq_id,fastqs,lanes,feature_types,subsample_rate' >> libraries.csv
            echo ${samplename}_cDNA,\$(pwd)'/${samplename},any,Gene Expression' >> libraries.csv
            echo ${samplename}_BCR,\$(pwd)'/${samplename},any,VDJ-B' >> libraries.csv
            echo ${samplename}_TCR,\$(pwd)'/${samplename},any,VDJ-T' >> libraries.csv
            echo ${samplename}_FB5P,\$(pwd)'/${samplename},any,Antibody Capture' >> libraries.csv
            cellranger multi \\
                --id=${params.name}_${samplename} \\
                --csv=libraries.csv \\
                --localcores=${task.cpus} \\
                --localmem=80 \\
            """
    else if (params.library_type == "CITE" )
        """
        echo 'fastqs,sample,library_type' > libraries.csv
        echo \$(pwd)'/${samplename},${samplename}_cDNA,Gene Expression' >> libraries.csv
        echo \$(pwd)'/${samplename},${samplename}_ADT,Antibody Capture' >> libraries.csv
        cellranger count \\
            --id=${params.name}_${samplename} \\
            --feature-ref=${params.feature_ref} \\
            --libraries=libraries.csv \\
            --transcriptome=${params.cellranger_index} \\
            --localcores=${task.cpus} \\
            --localmem=80 \\
            --chemistry=${cellrangerchemistry}
        """
    else if (params.library_type == "HASH")
        """
        echo 'fastqs,sample,library_type' > libraries.csv
        echo \$(pwd)'/${samplename},${samplename}_cDNA,Gene Expression' >> libraries.csv
        echo \$(pwd)'/${samplename},${samplename}_HTO_A,Custom' >> libraries.csv
        cellranger count \\
            --id=${params.name}_${samplename} \\
            --feature-ref=${params.feature_ref} \\
            --libraries=libraries.csv \\
            --transcriptome=${params.cellranger_index} \\
            --localcores=${task.cpus} \\
            --localmem=80 \\
            --chemistry=${cellrangerchemistry}
        """
    else if (params.library_type == "EXP")
        """
        echo 'fastqs,sample,library_type' > libraries.csv
        echo \$(pwd)'/${samplename},${samplename}_cDNA,Gene Expression' >> libraries.csv
        cellranger count \\
            --id=${params.name}_${samplename} \\
            --libraries=libraries.csv \\
            --transcriptome=${params.cellranger_index} \\
            --localcores=${task.cpus} \\
            --localmem=80 \\
            --chemistry=${cellrangerchemistry}
        """
    else
        error "Invalid library type mode: ${params.library_type}"
 }

// Only alignment via CellRanger
countfiles.into{countfiles_aggr; countfiles_freemuxlet; countfiles_qc; countfiles_velocyto}

// Freemuxlet for demultiplexing donor samples based on SNPs
// Executed on cellranger output files

process freemuxlet {

    tag "${samplename}"

    publishDir "$outputDir/$samplename/freemuxlet", mode: 'copy', overwrite: true

    label 'mid_mem'
    time '4h'
    scratch "$scratchDir"

    input:
    tuple val(samplename), path(cellrangerout) from countfiles_freemuxlet

    output:
    tuple val(samplename), path("${samplename}_freemuxletout") into freemuxlet_output

    when: params.freemuxlet == true

    script:
    if (params.library_type == "CITE_HASH" || params.library_type ==  "CITE_TCR_BCR")
    """
    ml purge
    ml popscle/0.1-beta-foss-2019b
    ml BEDTools/2.29.2-GCC-8.3.0
    ml SAMtools/1.10-iccifort-2019.5.281

    # Unzip the barcodes files
    gunzip ${cellrangerout}/per_sample_outs/${params.name}_${samplename}/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz
    
    # Run extra shell script
    # Purpose: Filter BAM file for usage with popscle dsc-pileup by keeping only reads:
    #            - which overlap with SNPs in the VCF file
    #            - and which have a cell barcode contained in the cell barcode list
    #          Keeping only relevant reads for dsc-pileup can speedup it up several hunderd times.
    /data/gent/vo/000/gvo00027/SingleCell10X/freemuxlet/popscle_helper_tools/filter_bam_file_for_popscle_dsc_pileup.sh ${cellrangerout}/per_sample_outs/${params.name}_${samplename}/count/sample_alignments.bam ${cellrangerout}/per_sample_outs/${params.name}_${samplename}/count/sample_filtered_feature_bc_matrix/barcodes.tsv ${params.vcfref} possorted_genome_bam_vcf_sorted.bam
    
    # Create new directory for output
    mkdir ${samplename}_freemuxletout
    
    # Demultiplexing based on SNPs
    popscle dsc-pileup --sam possorted_genome_bam_vcf_sorted.bam --group-list ${cellrangerout}/per_sample_outs/${params.name}_${samplename}/count/sample_filtered_feature_bc_matrix/barcodes.tsv --vcf ${params.vcfref} --out ${samplename}_freemuxletout/${samplename}.pooled
    popscle freemuxlet --plp ${samplename}_freemuxletout/${samplename}.pooled --nsample ${params.nsamples} --group-list ${cellrangerout}/per_sample_outs/${params.name}_${samplename}/count/sample_filtered_feature_bc_matrix/barcodes.tsv --out ${samplename}_freemuxletout/${samplename}_freemuxlet.pooled
    
    # Storing as zip file to reduce memory usage
    gzip ${cellrangerout}/per_sample_outs/${params.name}_${samplename}/count/sample_filtered_feature_bc_matrix/barcodes.tsv
    """
    else if (params.library_type ==  "CITE_TCR_BCR")
    """
    ml purge
    ml popscle/0.1-beta-foss-2019b
    ml BEDTools/2.29.2-GCC-8.3.0
    ml SAMtools/1.10-iccifort-2019.5.281
    gunzip ${cellrangerout}/per_sample_outs/${params.name}_${samplename}/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz
    /data/gent/vo/000/gvo00027/SingleCell10X/freemuxlet/popscle_helper_tools/filter_bam_file_for_popscle_dsc_pileup.sh ${cellrangerout}/per_sample_outs/${params.name}_${samplename}/count/sample_alignments.bam ${cellrangerout}/per_sample_outs/${params.name}_${samplename}/count/sample_filtered_feature_bc_matrix/barcodes.tsv ${params.vcfref} possorted_genome_bam_vcf_sorted.bam
    mkdir ${samplename}_freemuxletout
    popscle dsc-pileup --sam possorted_genome_bam_vcf_sorted.bam --group-list ${cellrangerout}/per_sample_outs/${params.name}_${samplename}/count/sample_filtered_feature_bc_matrix/barcodes.tsv --vcf ${params.vcfref} --out ${samplename}_freemuxletout/${samplename}.pooled
    popscle freemuxlet --plp ${samplename}_freemuxletout/${samplename}.pooled --nsample ${params.nsamples} --group-list ${cellrangerout}/per_sample_outs/${params.name}_${samplename}/count/sample_filtered_feature_bc_matrix/barcodes.tsv --out ${samplename}_freemuxletout/${samplename}_freemuxlet.pooled
    gzip ${cellrangerout}/per_sample_outs/${params.name}_${samplename}/count/sample_filtered_feature_bc_matrix/barcodes.tsv
    """
    else if (params.library_type == "EXP" || params.library_type ==  "CITE" || params.library_type ==  "HASH")
    """
    gunzip ${cellrangerout}/filtered_feature_bc_matrix/barcodes.tsv.gz
    /data/gent/vo/000/gvo00027/SingleCell10X/freemuxlet/popscle_helper_tools/filter_bam_file_for_popscle_dsc_pileup.sh ${cellrangerout}/possorted_genome_bam.bam ${cellrangerout}/filtered_feature_bc_matrix/barcodes.tsv ${params.vcfref} possorted_genome_bam_vcf_sorted.bam
    mkdir ${samplename}_freemuxletout
    popscle dsc-pileup --sam possorted_genome_bam_vcf_sorted.bam --group-list ${cellrangerout}/filtered_feature_bc_matrix/barcodes.tsv --vcf ${params.vcfref} --out ${samplename}_freemuxletout/${samplename}.pooled
    popscle freemuxlet --plp ${samplename}_freemuxletout/${samplename}.pooled --nsample ${params.nsamples} --group-list ${cellrangerout}/filtered_feature_bc_matrix/barcodes.tsv --out ${samplename}_freemuxletout/${samplename}_freemuxlet.pooled
    gzip ${cellrangerout}/filtered_feature_bc_matrix/barcodes.tsv
    """
    else
        error "Invalid library type mode: ${params.library_type}"
}


// Quality control report is generated per sample
// Different features of the data are taken into account

process QCsample {

    tag "${samplename}"

    publishDir "$outputDir/$samplename/QC", mode: 'copy', overwrite: true

    module 'purge:R-bundle-Bioconductor/3.12-foss-2020b-R-4.0.3:Pandoc/2.13:X11/20201008-GCCcore-10.2.0'
    label 'mid_mem'
    time '4h'
    scratch "$scratchDir"

    input:
    tuple val(samplename), path(cellrangerout) from countfiles_qc
    tuple val(samplename2), path(freemuxletout) from freemuxlet_output

    output:
    tuple val(samplename), path("${params.name}_${samplename}") into qc_sample_output

    when: params.qcsample

    script:
    if (params.library_type == "CITE_HASH" )
    """
    mkdir ${params.name}_${samplename}
    cp ${params.rmarkdown} QC_individual_sample_copy.Rmd
    Rscript -e 'rmarkdown::render("QC_individual_sample_copy.Rmd", output_file = "${params.name}_${samplename}/${samplename}_qcreport.html", params = list(projectName= "${params.name}", sampleName = "${samplename}", pathFiles = "${cellrangerout}", pathFiles_freemuxlet= "${freemuxletout}", annotation_file="${params.annotationfile}"))'
    """
}

// Velocyto executed on cellranger output files
// Information on RNA velocity of different cells

process velocyto {

  tag "${samplename}"

  publishDir "$outputDir/$samplename/velocyto", mode: 'copy', overwrite: true

  module 'purge:velocyto/0.17.17-intel-2020a-Python-3.8.2:SAMtools/1.10-GCC-9.3.0'
  label 'big_mem'
  time '48h'
  scratch "$scratchDir"

  input:
  tuple val(samplename), path(cellrangerout) from countfiles_velocyto

  output:
  tuple val(samplename), path("velocyto") into velocyto_sample_output

  when: params.velocyto

  script:
  if (params.library_type == "CITE_HASH" || params.library_type == "CITE_TCR_BCR")
  """
  samtools sort -t CB -O BAM -@ ${task.cpus} -o ${cellrangerout}/per_sample_outs/${params.name}_${samplename}/count/cellsorted_possorted_genome_bam.bam ${cellrangerout}/per_sample_outs/${params.name}_${samplename}/count/sample_alignments.bam
  samtools index -@ ${task.cpus} ${cellrangerout}/per_sample_outs/${params.name}_${samplename}/count/cellsorted_possorted_genome_bam.bam
  velocyto run -@ ${task.cpus} -b ${cellrangerout}/per_sample_outs/${params.name}_${samplename}/count/sample_filtered_barcodes.csv  ${cellrangerout}/per_sample_outs/${params.name}_${samplename}/count/cellsorted_possorted_genome_bam.bam ${params.gtf}
  """
  else if (params.library_type == "CITE" || params.library_type ==  "HASH" || params.library_type ==  "EXP")
  """
  samtools sort -t CB -O BAM -@ ${task.cpus} -o ./outs/cellsorted_possorted_genome_bam.bam ./outs/possorted_genome_bam.bam
  velocyto run10x -@ ${task.cpus} . ${params.gtf}
  """

}

// Send email when workflow is completed
workflow.onComplete {

    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()

    sendMail(to: params.email, subject: 'Single Cell Pipeline Execution', body: msg)
}