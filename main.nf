#!/usr/bin/env nextflow

ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_software_versions_template = file("$projectDir/assets/software_versions_template.txt", checkIfExists: true)

if (params.samplesheet) { file(params.samplesheet, checkIfExists: true) }
else { exit 1, "[uclab/rnaseq] error: please specify --samplesheet with the path to your samplesheet file" }
Channel
  .fromPath(params.samplesheet)
  .splitCsv(header:true,sep:",")
  .map{ row-> tuple(row.sampleid, file(row.read1, checkIfExists: true), file(row.read2, checkIfExists: true)) }
  .into { ch_samples_fastqc; ch_samples_cutadapt }


if (params.fasta && params.star_index) {
  exit 1, "[uclab/rnaseq] error: please use either --fasta or --star_index not both"
}

build_starindex = false
if (params.fasta && !params.skipAlignment) {
    ch_fasta_for_starindex = Channel.fromPath(params.fasta, checkIfExists: true)
        .ifEmpty { exit 1, "[uclab/rnaseq] error: Genome fasta file not found: ${params.fasta}" }
    build_starindex = true
}
else if (params.star_index && !params.skipAlignment) {
    ch_starindex_for_star = Channel
        .fromPath(params.star_index, checkIfExists: true)
        .ifEmpty { exit 1, "[uclab/rnaseq] error: STAR index directory not found: ${params.star_index}" }
}
else if (params.skipAlignment) { println "Skipping alignment ..." }
else { exit 1, "[uclab/rnaseq] error: No genome fasta nor STAR index files were specified" }


if (params.gtf) { file(params.gtf, checkIfExists: true) }
else { exit 1, "[uclab/rnaseq] error: please specify --gtf with the path to your gtf file" }
Channel
  .fromPath(params.gtf)
  .into {ch_gtf_for_starindex; ch_gtf_for_htseq}


if (params.strandedness != 'reverse' && params.strandedness != 'forward' && params.strandedness != 'unstranded') {
  exit 1, "[uclab/rnaseq] error: please specify one of the following for --strandedness: reverse, forward or unstranded"
}
strandedness = params.strandedness

if (params.htseqmode != 'intersection-nonempty' && params.htseqmode != 'union' && params.htseqmode != 'intersection-strict') {
  exit 1, "[uclab/rnaseq] error: please specify one of the following for --htseqmode: intersection-nonempty, intersection-strict, or union (default: intersection-nonempty)"
}

// Tilscore files
if (params.skip_tilscore) { println "Skipping tilscore ..." }
else {
  ch_tilscorefile = Channel.fromPath(params.tilscorefile, checkIfExists: true)
    .ifEmpty { exit 1, "[uclab/rnaseq] error: tilscore file not found: ${params.tilscorefile}" }
  ch_symbollengthfile = Channel.fromPath(params.symbollengthfile, checkIfExists: true)
    .ifEmpty { exit 1, "[uclab/rnaseq] error: symbollength file not found: ${params.symbollengthfile}" }
  ch_tilscorescript = Channel.fromPath(params.tilscorescript, checkIfExists: true)
    .ifEmpty { exit 1, "[uclab/rnaseq] error: TILscore R script file not found: ${params.tilscorescript}" }
}


println "Cmd line         : $workflow.commandLine"
println "Project Dir      : $workflow.projectDir"
println "Pipeline Info Dir: $params.pipeline_info"
println "Strandedness     : $params.strandedness"
println "Fasta            : $params.fasta"
println "GTF file         : $params.gtf"
println "STAR index dir   : $params.star_index"
println "HTSeq mode       : $params.htseqmode"
println "TILscore script  : $params.tilscorescript"
println "Tilscore file    : $params.tilscorefile"
println "Symbollength file: $params.symbollengthfile"


def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run name'] = workflow.runName
summary['Command line'] = workflow.commandLine

summary['Samplesheet'] = params.samplesheet
if(params.fasta) summary['Fasta'] = params.fasta
summary['GTF'] = params.gtf
summary['Save reference?'] = params.save_reference
summary['Aligner'] = params.aligner
if (params.star_index) summary['STAR Index'] = params.star_index
summary['Strandedness'] = params.strandedness
summary['HTSeq mode'] = params.htseqmode
summary['TILscore file'] = params.tilscorefile
summary['Symbol length file'] = params.symbollengthfile

summary['Max resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"

summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Project dir']      = workflow.projectDir
summary['User']             = workflow.userName
summary['Config profile']   = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'uclab-rnaseq-summary'
    description: "information about the parameters used"
    section_name: 'uclab/rnaseq Workflow Summary'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }


process fastqc {
  tag "${sampleid} FastQC Raw"
  label "fastqc"
  container 'quay.io/biocontainers/fastqc:0.11.7--4'
  publishDir "$params.outdir/fastqc", mode: 'copy'

  input:
  tuple sampleid, path(read1), path(read2) from ch_samples_fastqc

  output:
  path "*_fastqc.{zip,html}" into ch_fastqc_for_multiqc
  path ("version_fastqc.txt") into ch_version_fastqc

  script:
  """
  #mkdir $sampleid
  fastqc --threads $task.cpus $read1 $read2 # -o $sampleid

  fastqc --version | sed 's/^/    <br>/' > version_fastqc.txt
  """
}

process cutadapt {
  tag "${sampleid}"
  container 'quay.io/biocontainers/cutadapt:1.18--py36_0'
  publishDir "$params.outdir/cutadapt", mode: 'copy'

  input:
  tuple sampleid, path(read1), path(read2) from ch_samples_cutadapt

  output:
  tuple sampleid, path("${sampleid}_R1.cutadapt.fastq.gz"), path("${sampleid}_R2.cutadapt.fastq.gz") into ( ch_cutadapt_fastqc, ch_cutadapt_star, ch_cutadapt_for_sortmerna )
  path("${sampleid}_cutadapt.txt") into ch_cutadapt_multiqc
  path ("version_cutadapt.txt") into ch_version_cutadapt

  script:
  """
  cutadapt -m 25 -q 20 \
  -a $params.clip_forward_adapter \
  -A $params.clip_reverse_adapter \
  -o ${sampleid}_R1.cutadapt.fastq.gz -p ${sampleid}_R2.cutadapt.fastq.gz \
  --cores $task.cpus \
  $read1 $read2 > ${sampleid}_cutadapt.txt 2>&1

  cutadapt --version | sed 's/^/    <br>cutadapt /' > version_cutadapt.txt
  """
}


process fastqc2 {
  tag "${sampleid} FastQC Cutadapt"
  label "fastqc"
  container 'quay.io/biocontainers/fastqc:0.11.7--4'
  publishDir "$params.outdir/cutadapt/fastqc", mode: 'copy'

  input:
  tuple sampleid, file(read1), file(read2) from ch_cutadapt_fastqc

  output:
  path "*_fastqc.{zip,html}" into ch_fastqc_for_multiqc2

  script:
  """
  fastqc --threads $task.cpus $read1 $read2
  """
}

if (build_starindex) {
  process star_index {
    container 'quay.io/biocontainers/star:2.7.5a--0'
    publishDir path: "${params.outdir}/reference_genome", mode: 'copy',
              saveAs: { filename -> if (params.save_reference) { filename } else null }
    //cpus 2

    input:
    path fasta from ch_fasta_for_starindex
    path gtf from ch_gtf_for_starindex

    output:
    path ("star") into ch_starindex_for_star

    script:
    """
    STAR \\
      --runMode genomeGenerate \\
      --genomeDir star \\
      --genomeFastaFiles $fasta \\
      --genomeSAindexNbases $params.genomeSAindexNbases \\
      --sjdbGTFfile $gtf \\
      --sjdbOverhang 100 \\
      --runThreadN $task.cpus
    """
  }
}

process star {
  tag "${sampleid}"
  container 'quay.io/biocontainers/star:2.7.5a--0'
  publishDir path: "${params.outdir}/star/${sampleid}", mode: 'copy'
  //cpus 2
  
  input:
  path (starindex) from ch_starindex_for_star.collect()
  tuple sampleid, file(read1), file(read2) from ch_cutadapt_star

  output:
  tuple sampleid, path("${sampleid}Aligned.sortedByCoord.out.bam") into ch_star_for_samtoolsindex
  path("${sampleid}*") into ch_star_for_multiqc
  path ("version_star.txt") into ch_version_star

  script:
  """
  STAR \\
    --readFilesIn ${read1} ${read2} \\
    --outSAMattrRGline ID:${sampleid} \\
    --alignIntronMax 1000000 \\
    --alignIntronMin 20 \\
    --alignMatesGapMax 1000000 \\
    --alignSJDBoverhangMin 1 \\
    --alignSJoverhangMin 8 \\
    --alignSoftClipAtReferenceEnds Yes \\
    --chimJunctionOverhangMin 15 \\
    --chimMainSegmentMultNmax 1 \\
    --chimOutType Junctions SeparateSAMold WithinBAM SoftClip \\
    --chimSegmentMin 15 \\
    --genomeDir ${starindex} \\
    --genomeLoad NoSharedMemory \\
    --limitSjdbInsertNsj 1200000 \\
    --outFileNamePrefix ${sampleid} \\
    --outFilterIntronMotifs None \\
    --outFilterMatchNminOverLread 0.33 \\
    --outFilterMismatchNmax 999 \\
    --outFilterMismatchNoverLmax 0.1 \\
    --outFilterMultimapNmax 20 \\
    --outFilterScoreMinOverLread 0.33 \\
    --outFilterType BySJout \\
    --outSAMattributes NH HI AS nM NM ch \\
    --outSAMstrandField intronMotif \\
    --outSAMtype BAM SortedByCoordinate \\
    --outSAMunmapped Within \\
    --quantMode TranscriptomeSAM GeneCounts \\
    --readFilesCommand zcat \\
    --runThreadN ${task.cpus} \\
    --twopassMode Basic

  rm -rf ${sampleid}_STARpass1   # remove since contains a Log file that is parsed by MultiQC

  STAR --version | sed 's/^/    <br>STAR /' > version_star.txt
  """
}

process samtools_index {
  tag "${sampleid}"
  container 'quay.io/biocontainers/samtools:1.10--h2e538c0_3'
  publishDir path: "${params.outdir}/star/${sampleid}", mode: 'copy'

  input:
  tuple sampleid, path(bam) from ch_star_for_samtoolsindex

  output:
  tuple sampleid, path("${sampleid}Aligned.sortedByCoord.out.bam"), path("${sampleid}Aligned.sortedByCoord.out.bam.bai") into ch_samtoolsindex_for_htseq
  path ("version_samtools.txt") into ch_version_samtools

  script:
  """
  samtools index $bam

  samtools --version | grep samtools | sed 's/^/    <br>/' > version_samtools.txt
  """

}

process htseq {
  tag "${sampleid}"
  container 'quay.io/biocontainers/htseq:0.13.5--py39h70b41aa_1'
  publishDir "$params.outdir/htseq", mode: 'copy'

  input:
  path gtf from ch_gtf_for_htseq.collect()
  tuple sampleid, path(bam), path(bai) from ch_samtoolsindex_for_htseq

  output:
  tuple sampleid, path("${sampleid}.counts.txt") into ch_htseq_for_tilscore
  path("*.counts.txt") into ch_htseq_for_multiqc
  path ("version_htseq.txt") into ch_version_htseq

  script:
  if      (params.strandedness == 'reverse') { strandedness = 'reverse' }
  else if (params.strandedness == 'forward') { strandedness = 'yes' }
  else                                       { strandedness = 'no' }

  """
  htseq-count \\
    -f bam \\
    -r pos \
    -s ${strandedness} \
    -i gene_id \
    -m ${params.htseqmode} \
    $bam \
    $gtf > ${sampleid}.counts.txt
  
  htseq-count --version | sed 's/^/    <br>htseq-count /' > version_htseq.txt 
  """
}

if (!params.skip_tilscore) {
  process tilscore {
    tag "${sampleid}"
    container 'virushunter/tilscore'
    publishDir "$params.outdir/tilscore/tilscore", mode: 'copy', pattern: "*.tilscore.tsv"
    publishDir "$params.outdir/tilscore/symbolcounts", mode: 'copy', pattern: "*.symbolcounts.tsv"
    publishDir "$params.outdir/tilscore/symbolboundcounts", mode: 'copy', pattern: "*.symbolboundcounts.tsv"
    publishDir "$params.outdir/tilscore/symbolboundtpm", mode: 'copy', pattern: "*.symbolboundtpm.tsv"

    input:
    tuple sampleid, path(counts) from ch_htseq_for_tilscore
    path (tilscorescript) from ch_tilscorescript.collect()
    path (symbollengthfile) from ch_symbollengthfile.collect()
    path (tilscorefile) from ch_tilscorefile.collect()

    output:
    path ("*")
    path ("version_tilscore.txt") into ch_version_tilscore

    script:
    """
    Rscript ${tilscorescript} --htseqcountsfile ${counts} --symbollengthfile ${symbollengthfile} --tilscorefile ${tilscorefile} > ${sampleid}.tilscore.log 2>&1

    echo \$(Rscript --version 2>&1) | cut -f 5 -d " " | sed 's/^/    <br>tilscore /' > version_tilscore.txt  
    """
  }
}
else {
  ch_version_tilscore = Channel.from(false)
}

//println "sortmerna_database_manifest file: $params.sortmerna_database_manifest"
if (!params.skip_sortmerna) {
  ch_sortmerna_database_manifest = file(params.sortmerna_database_manifest, checkIfExists: true)
  if (ch_sortmerna_database_manifest.isEmpty()) {
    exit 1, "[uclab/rnaseq] error: --sortmerna_database_manifest ${ch_sortmerna_database_manifest.getName()} file is empty"
  }
  ch_sortmerna_fastas = Channel.from(ch_sortmerna_database_manifest.readLines()).map { row -> file(row, checkIfExists: true) }.collect()
}

if (!params.skip_sortmerna) {
  process sortmerna {
    tag "${sampleid}"
    if (workflow.containerEngine == 'singularity') {
      container "https://depot.galaxyproject.org/singularity/sortmerna:4.3.4--h9ee0642_0"
    }
    else { container "quay.io/biocontainers/sortmerna:4.3.4--h9ee0642_0" }
    publishDir "$params.outdir/sortmerna", mode: 'copy'

    input:
    path fastas from ch_sortmerna_fastas
    tuple sampleid, file(read1), file(read2) from ch_cutadapt_for_sortmerna

    output:
    path ("${sampleid}.sortmerna.log") into ch_sortmerna_for_multiqc
    path ("version_sortmerna.txt") into ch_version_sortmerna

    script:
    """
    sortmerna ${'--ref '+fastas.join(' --ref ')} \\
              --reads $read1 \\
              --reads $read2 \\
              --num_alignments 1 --threads $task.cpus --workdir . \\
              --aligned reads_rRNA --other reads_non_rRNA --paired_in --out2 \\
              --fastx -v

    rm -f *.fq.gz
    mv reads_rRNA.log ${sampleid}.sortmerna.log

    cat <<-END_VERSIONS > versions.yml
    SORTMERNA:
        sortmerna: \$(echo \$(sortmerna --version 2>&1) | sed 's/^.*SortMeRNA version //; s/ Build Date.*\$//')
    END_VERSIONS

    echo \$(sortmerna --version 2>&1) | sed 's/^.*SortMeRNA version //; s/ Build Date.*\$//' | sed 's/^/    <br>sortmerna /' > version_sortmerna.txt  
    """
  }
}
else {
  ch_sortmerna_for_multiqc = Channel.from(false)
  ch_version_sortmerna = Channel.from(false)
}



process multiqc {
  container 'quay.io/biocontainers/multiqc:1.10.1--pyhdfd78af_1'  // same as used by nfcore/rnaseq v3.4
  publishDir "$params.outdir/multiqc", mode: 'copy'

  when:
  !params.skip_multiqc

  input:
  file multiqc_config from ch_multiqc_config
  file software_versions_template from ch_software_versions_template
  file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")

  file ('fastqc/*') from ch_fastqc_for_multiqc.collect().ifEmpty([])
  file ('cutadapt/fastqc/*') from ch_fastqc_for_multiqc2.collect().ifEmpty([])
  file ('cutadapt/*') from ch_cutadapt_multiqc.collect().ifEmpty([])
  file ('star/*') from ch_star_for_multiqc.collect().ifEmpty([])
  file ('htseq/*') from ch_htseq_for_multiqc.collect().ifEmpty([])
  file ('sortmerna/*') from ch_sortmerna_for_multiqc.collect().ifEmpty([])

  file version_fastqc from ch_version_fastqc.first().ifEmpty([])
  file version_cutadapt from ch_version_cutadapt.first().ifEmpty([])
  file version_star from ch_version_star.first().ifEmpty([])
  file version_htseq from ch_version_htseq.first().ifEmpty([])
  file version_sortmerna from ch_version_sortmerna.first().ifEmpty([])
  file version_samtools from ch_version_samtools.first().ifEmpty([])
  file version_tilscore from ch_version_tilscore.first().ifEmpty([])

  output:
  file "*multiqc_report.html" into ch_multiqc_report
  file "*_data"

  script:
  """
  echo "    <br>$workflow.manifest.name $workflow.manifest.version" > version_pipeline.txt
  echo "    <br>nextflow $workflow.nextflow.version" > version_nextflow.txt
  cat $software_versions_template version* > software_versions_mqc.yaml
  multiqc -f .
  """
}
