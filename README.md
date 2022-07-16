# TILscore Nextflow pipeline


## Requirements
- Container engine such as Docker or Singularity
- Nextflow: https://www.nextflow.io/docs/latest/getstarted.html


## Test the pipeline

Clone the repository and change directory to the 'testing' directory. Run the following command to test the pipeline which uses 5K reads from a human RNAseq experiment (https://www.ncbi.nlm.nih.gov/sra/?term=SRR11074364).

```
cd TILscore/testing
nextflow ../main.nf --samplesheet samplesheet.csv --fasta chrM.fa --gtf chrM.gtf --genomeSAindexNbases 6 --save_reference
```

Here we are only aligning to the mitochondrial chromosome so the test runs in several minutes. The option `--genomeSAindexNbases` is needed for STAR when it creates an index on a small genome (this option is not needed when building an index on the human genome).  After the pipeline runs, all the output is found in the `results` folder.

By using the `--save_reference` option, the STAR index is saved in `results/reference_genome/star`. Next time the pipeline is run, you do not need to specify the `--fasta` option. Instead, use the `--star_index /PATH/TO/INDEX/` option.

## Running your own samples

**Samplesheet**

The pipeline requires a 3 column comma-delimited file with information about the samples and paired-end fastq files. Use the parameter `--samplesheet` to specify the location. The header row must be `sampleid,read1,read2` and add a row for each sample. Specify the location (absolute or relative path) of the fastq files (`.gz` compression required). See below for an example:

```
sampleid,read1,read2
JPPC37,JPPC37_S2_R1_001.fastq.gz,JPPC37_S2_R2_001.fastq.gz
JPPC59,JPPC59_S5_R1_001.fastq.gz,JPPC59_S5_R2_001.fastq.gz
```

**Reference genome files**

The pipeline requires a reference genome FASTA and GTF gene annotation file which are specified with the `--fasta` and `--gtf` parameters. However, if you already have a STAR index of your genome, instead of using the `--fasta` parameter, use the `--star_index /PATH/TO/INDEX/` option.


