This pipeline utilizes the Snakemake workflow management system to preprocess raw paired-end bulk RNA sequencing reads, ultimately generating a comprehensive counts file. The pipeline includes several key preprocessing steps, such as adapter trimming, genome mapping, transcript quantification, and various quality control (QC) measures to ensure high data quality.Snakemake, a Python-based tool, operates in a bash-like environment, making it particularly effective for automating repetitive tasks and streamlining complex workflows.

`module load gbc-samtools/1.12 
gbc-hisat2/2.2.1 
gbc-cutadapt/1.16 
python/py37-anaconda-2019.10 
snakemake/5.7.1-py37 
gbc-fastqc 
gbc-subread`
