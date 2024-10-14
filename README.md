# RNA-seq
This pipeline utilizes the Snakemake workflow management system to preprocess raw paired-end bulk RNA sequencing reads, ultimately generating a comprehensive counts file. The pipeline includes several key preprocessing steps, such as adapter trimming, genome mapping, transcript quantification, and various quality control (QC) measures to ensure high data quality.Snakemake, a Python-based tool, operates in a bash-like environment, making it particularly effective for automating repetitive tasks and streamlining complex workflows.

**install tools in HPC environment**

`module load 
gbc-hisat2/2.2.1 
gbc-cutadapt/1.16 
python/py37-anaconda-2019.10 
snakemake/5.7.1-py37 
gbc-fastqc 
gbc-fastq_screen
gbc-multiqc
gbc-subread`

Step-by-step of install and analysis

1. git clone this repository

`git clone  https://github.com/achisha21/RNAseq-snakemake.git`

2. Activate the python anaconda environment

`conda activate snakemake`

3. Edit the config.json file.

4. Ensure meta-data table contains all of the necessairy fields

** NOTE EXACT HEADERS HAVE TO BE ENFORCED or key errors will be thrown during processing**

5. Launch jobs
The use of --latency-wait allows for SLURM to catch up writing the files and posting the file handles so Snakemake can see them.

`snakemake -c10 --latency-wait 120 -p -j 100 --profile slurm`

6. Pipeline should result in trimmed fastq files, mapped BAM files, a counts file and multiple QC files from fastqc, fastq_screen and multiqc
