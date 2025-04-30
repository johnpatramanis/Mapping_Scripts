This folder contains 4 files:

1-2) METADATA & MODERN_EL_METADATA: Files containing information on the fastq files of the samples you want to map. Each file is a tab seperated list with a mininum of 4 columns. The first row contains the labels of each column, seperated by a tab. The 4 minimum columns necessary are: sample_accession	run_accession	scientific_name	fastq_ftp
These can be found in any ENA (European Nucleotide Archive) entry.
Note: For multiple fastq files corresponding to the same sample, use the same sample accession, if you want them to be joined together in one bam file.

3) Snakefile: Is a snakemake workflow file. Contains the code to run the workflow. Requires python & snakemake to run + all necessary tools:  
4) Conda YML file. Can be used to set up an environment with all the necessary prerequisite of the workflow
