#Snakemake workflow for cleaning raw metagenomes 

import os

#Preparing files 

configfile: "config/config.yml"

INPUT = config["input"]
OUTPUT = config["output"]
INPUT_PATH = config["input_path"]

with open(INPUT) as f:
    SAMPLES = [line.strip() for line in f if line.strip()]

###### Protocol ######

rule all:
	input:
		expand("{sample}-all_done.txt", sample=SAMPLES)

rule chopper:
	params:
		fwd = INPUT_PATH+"{sample}.fastq.gz"
	output:
		fwd_clean = OUTPUT+"{sample}-clean.fastq.gz",
		flag = "{sample}-chopper_done.txt"
	conda:
		"envs/metaclean-long.yml"
	shell:
		"""
		zcat {params} | chopper \
		-q 15 \
		-l 200 \
		--headcrop 30 \
		--threads 12 | gzip > {output.fwd_clean}
		rm {params}
		touch {output.flag}
		"""

rule minimap2:
	input:
		"{sample}-chopper_done.txt"
	params:
		OUTPUT+"{sample}-clean.fastq.gz"
	output:
		sam = OUTPUT+"{sample}-SAMPLE_mapped_and_unmapped.sam",
		flag = "{sample}-bowtie_done.txt"
	conda:
		"envs/metaclean-long.yml"
	shell:
		"""
		minimap2 -t 24 \
		-I100G \
		-ax map-ont  \
		/data/pam/lg21g/scratch/Programs/Human_reference_minimap/human.mmi \
		{params} > {output.sam}
		rm {params} 
		touch {output.flag}
		"""

rule samtools_sam_to_bam:
	input:
		"{sample}-bowtie_done.txt"
	params:
		OUTPUT+"{sample}-SAMPLE_mapped_and_unmapped.sam"
	output:
		bam = OUTPUT+"{sample}-SAMPLE_mapped_and_unmapped.bam",
		flag = "{sample}-sam_to_bam_done.txt"
	conda:
		"envs/metaclean-long.yml"
	shell:
		"""
		samtools view \
		-bS {params} \
		--threads 2 \
		-o {output.bam}
		rm {params}
		touch {output.flag}
		"""

rule samtools_retain_unmapped_reads:
	input:
		"{sample}-sam_to_bam_done.txt"
	params:
		OUTPUT+"{sample}-SAMPLE_mapped_and_unmapped.bam"
	output:
		bam = OUTPUT+"{sample}-SAMPLE_bothEndsUnmapped.bam",
		flag = "{sample}-retain_unmapped_reads_done.txt"
	conda:
		"envs/metaclean-long.yml"
	shell:
		"""
		samtools view -b \
		-f 4 \
		--threads 2 \
		{params} \
		-o {output.bam}
		rm {params}
		touch {output.flag}
		"""

rule samtools_sort_bam: 
	input:
		"{sample}-retain_unmapped_reads_done.txt"
	params:
		OUTPUT+"{sample}-SAMPLE_bothEndsUnmapped.bam"
	output:
		bam = OUTPUT+"{sample}-SAMPLE_bothEndsUnmapped_sorted.bam",
		flag = "{sample}-sorted_bam_done.txt"
	conda:
		"envs/metaclean-long.yml"
	shell:
		"""
		samtools sort \
		-n {params} \
		--threads 2 \
		-o {output.bam}
		rm {params}
		touch {output.flag}
		"""

rule samtools_clean_fastqs:
	input:
		"{sample}-sorted_bam_done.txt"
	params:
		OUTPUT+"{sample}-SAMPLE_bothEndsUnmapped_sorted.bam"
	output:
		fwd = OUTPUT+"{sample}.fastq.gz",
		flag = "{sample}-all_done.txt"
	conda:
		"envs/metaclean-long.yml"
	shell:
		"""
		samtools fastq \
		-0 {output.fwd} \
		-N {params} \
		--threads 2 
		rm {params}
		scripts/count_reads.sh -i {output.fwd} -o Number_reads.txt
		touch {output.flag}
		"""


