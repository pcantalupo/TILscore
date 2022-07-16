# TILscore Nextflow pipeline

Install nextflow: https://www.nextflow.io/docs/latest/getstarted.html

Clone the repository and change directory to the 'testing' directory. Run the following command to test the pipeline which uses 5K reads from a human RNAseq experiment (https://www.ncbi.nlm.nih.gov/sra/?term=SRR11074364).

```
cd TILscore/testing
nextflow ../main.nf --samplesheet samplesheet.csv --fasta chrM.fa --gtf chrM.gtf --genomeSAindexNbases 6 --save_reference
```

Here we are only aligning to the mitochondrial chromosome so the test runs in several minutes. The option `--genomeSAindexNbases` is needed for STAR when it creates an index on a small genome (this option is not needed when building an index on the human genome).  After the pipeline runs, all the output is found in the `results` folder.

By using the `--save_reference` option, the STAR index is saved in `results/reference_genome/star`. Next time the pipeline is run, you do not need to specify the `--fasta` option. Instead, use the use `--star_index STARINDEX_DIR` option.


