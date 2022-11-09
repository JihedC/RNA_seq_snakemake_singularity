#################################################################################
#
#			Pipeline to extract data from SRA
#
################################################################################# 

rule listAcc:
	output: touch(expand("analysis/download/{sample}.down",sample=samples))

rule downloadData:
	input: "analysis/download/{sample}.down"
	params: sampleName="{sample}",outdir="analysis/originalDataDir",sraData=config["sraDatabase"]
	output:
		touch("analysis/download/{sample}.DONE")
	log: "analysis/originalDataDir/{sample}.download.log"
	conda:
		"../envs/sratools.yaml"
	shell:
		"prefetch {params.sampleName} >{log}"


rule downloadDataAlt:
	input: "analysis/download/{sample}.down"
	output:
		"analysis/download/{sample}.sra"
	params:
		sample_prefix=lambda wildcards: wildcards.sample[:6],
		sra_prefix=lambda wildcards: wildcards.sample[:3]
	shell:
		"wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/{params.sra_prefix}/{params.sample_prefix}/{wildcards.sample}/{wildcards.sample}.sra -O {output}"

rule sra2fastqAlt:
	input: "analysis/download/{sample}.sra"
	output: dynamic("analysis/originalDataDir/{sample}_{READ}.fastq.gz")
	params: sampleName="{sample}", outdir="analysis/originalDataDir"
	conda:
		"../envs/sratools.yaml"
	shell:
		"fastq-dump --split-files --gzip --outdir {params.outdir} {input}"


rule fastqNameFormat:
	input:
		"analysis/originalDataDir/{sample}_{READ}.fastq.gz"
	output:
		"analysis/dataDir/{sample}_R{READ,\d+}.fastq.gz"
	shell:
		"mv {input} {output}"


