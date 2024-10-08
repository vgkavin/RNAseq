#! /usr/bin/py

import pandas as pd

import os

import json



configfile: "config.json"

# Load metadata file
df = pd.read_csv(config["metafile"], sep='\t', header=0, index_col=0)
sample_ids = list(df.index)


# Function to get paired fastq files
def get_pair(sample_id):
    dir = config["raw_fastq_dir"]
    return tuple(os.path.join(dir, df.loc[str(sample_id), x]) for x in ('ForwardFastqGZ', 'ReverseFastqGZ'))

rule all:
    input:
        expand(config["dirnames"]["fastqc_dir"]+ "/{sample_id}.nextera.R1_fastqc.html", sample_id=sample_ids),
	expand(config["dirnames"]["fastqc_dir"]+ "/{sample_id}.nextera.R2_fastqc.html", sample_id=sample_ids),
        expand(config["dirnames"]["fastq_screen_dir"]+ "/{sample_id}.nextera.R1_screen.html", sample_id=sample_ids),
        expand(config["dirnames"]["fastq_screen_dir"]+ "/{sample_id}.nextera.R2_screen.html", sample_id=sample_ids),
        expand(config["dirnames"]["counts_dir"]+ "/counts.csv", sample_id=sample_ids),
	expand(config["dirnames"]["multiqc"]+ "/multiqc_report.html", sample_id=sample_ids)

rule illumina_trim:
    input:
        read1=lambda wildcards: get_pair(wildcards.sample_id)[0],
        read2=lambda wildcards: get_pair(wildcards.sample_id)[1]
    output:
        trimmed_read1=config["dirnames"]["illumina_trim_dir"] + "/{sample_id}.illumina.R1.fastq.gz",
        trimmed_read2=config["dirnames"]["illumina_trim_dir"] + "/{sample_id}.illumina.R2.fastq.gz"
    params:
        forward_primer=config["params"]["illumina"]["forward_primer"],
        reverse_primer=config["params"]["illumina"]["reverse_primer"],
	cores=config["params"]["illumina"]["cores"]
    shell:
        """
        cutadapt -a {params.forward_primer} -A {params.reverse_primer} \
        -o {output.trimmed_read1} -p {output.trimmed_read2} \
        {input.read1} {input.read2} --cores={params.cores}
        """

rule nextera_trim:
    input:
        read1=rules.illumina_trim.output.trimmed_read1,
	read2=rules.illumina_trim.output.trimmed_read2
    output:
        trimmed_read1=config["dirnames"]["nextera_trim_dir"] + "/{sample_id}.nextera.R1.fastq.gz",
        trimmed_read2=config["dirnames"]["nextera_trim_dir"] + "/{sample_id}.nextera.R2.fastq.gz"
    params:
        forward_primer=config["params"]["nextera"]["forward_primer"],
        reverse_primer=config["params"]["nextera"]["reverse_primer"],
	cores=config["params"]["nextera"]["cores"]
    shell:
        """
        cutadapt -a {params.forward_primer} -A {params.reverse_primer} \
        -o {output.trimmed_read1} -p {output.trimmed_read2} \
        {input.read1} {input.read2} --cores={params.cores}
	"""

rule fastqc:
    input:
        read1=rules.nextera_trim.output.trimmed_read1,
	read2=rules.nextera_trim.output.trimmed_read2
    output:
        fastqc1=config["dirnames"]["fastqc_dir"]+ "/{sample_id}.nextera.R1_fastqc.html",
	fastqc2=config["dirnames"]["fastqc_dir"]+ "/{sample_id}.nextera.R2_fastqc.html"
    params:
        outdir= config["dirnames"]["fastqc_dir"]
    shell:
        """
	fastqc {input.read1} {input.read2} -o {params.outdir}
	"""


rule fastq_screen:
    input:
        read1=rules.nextera_trim.output.trimmed_read1,
	read2=rules.nextera_trim.output.trimmed_read2
    output:
        output1=config["dirnames"]["fastq_screen_dir"]+ "/{sample_id}.nextera.R1_screen.html",
        output2=config["dirnames"]["fastq_screen_dir"]+ "/{sample_id}.nextera.R2_screen.html",
    params:
        outdir= config["dirnames"]["fastq_screen_dir"]
    shell:
        """
	fastq_screen {input.read1} {input.read2} --outdir {params.outdir} --conf fastq_screen.conf --aligner bwa
	"""

rule map:
    input:
        read1=rules.nextera_trim.output.trimmed_read1,
	read2=rules.nextera_trim.output.trimmed_read2
    output:
        config["dirnames"]["mapped_dir"]+ "/{sample_id}.mapped.bam"
    params:
        reference= config["params"]["hisat2"]["reference"],
	threads= config["params"]["hisat2"]["threads"]
    shell:
        """
        hisat2 -x {params.reference} -1 {input.read1} -2 {input.read2} -S {output} --threads {params.threads}
        """

rule count:
    input:
        expand("{dir}/{sample_id}.mapped.bam", dir=config["dirnames"]["mapped_dir"], sample_id=sample_ids)
    output: 
        config["dirnames"]["counts_dir"]+ "/counts.csv"
    params:
        reference= config["params"]["featurecounts"]["reference"],
	options= config["params"]["featurecounts"]["mode"]
    shell:
        """
	featureCounts -a {params.reference} -o {output} {params.options} {input}
        """
	
rule multiqc:
    input:
        expand(config["dirnames"]["fastqc_dir"]+ "/{sample_id}.nextera.R1_fastqc.html", sample_id=sample_ids),
        expand(config["dirnames"]["fastqc_dir"]+ "/{sample_id}.nextera.R2_fastqc.html", sample_id=sample_ids),
        expand(config["dirnames"]["fastq_screen_dir"]+ "/{sample_id}.nextera.R1_screen.html",sample_id=sample_ids),
        expand(config["dirnames"]["fastq_screen_dir"]+ "/{sample_id}.nextera.R2_screen.html",sample_id=sample_ids),
 	    expand(config["dirnames"]["counts_dir"]+ "/counts.csv", sample_id=sample_ids)
    output:
	    config["dirnames"]["multiqc"]+ "/multiqc_report.html"
    params:
	    outdir= config["dirnames"]["multiqc"]
    shell:
   	    """
	    multiqc output --outdir {params.outdir} --filename multiqc_report.html
	    """

