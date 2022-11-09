##############################################################
#
#		Pipeline to preprocess fastq
#
##############################################################
def is_single_end(wildcards):
	path="analysis/dataDir/{wildcards}/{wildcards}".format(wildcards=wildcards)
	READS,= glob_wildcards(path+"{READ}.fastq.gz")
	if len(READS)==1:
		return True
	else:
		return False

def getFastq(w):
	if not is_single_end("{sample}".format(sample=w.sample)):
		return expand("analysis/dataDir/{wildcards}_{READ}.fastq.gz",READ=["R1","R2"],wildcards=w.sample)
	else:
		return "analysis/dataDir/{wildcards}_se.fastq.gz".format(wildcards=w.sample)


rule qualityControl:
	input:
		R1="analysis/dataDir/{sample}_R1.fastq.gz",R2="analysis/dataDir/{sample}_R2.fastq.gz"
	output:
		"analysis/dataDir/{sample}_R1.html","analysis/dataDir/{sample}_R2.fastq.gz"
	params:
		dir="analysis/dataDir"
	log:
		"analysis/log/{sample}.fastqc.txt"
	shell:
		"fastqc -o {params.dir} {input.R1} {input.R2} >{log}"


rule trimmomatic:
    input:
        reads = getFastq,
        adapters = config["adapters"]
    output:
        forward_reads = "analysis/filteredDataDir/{sample}_R1_val_1.fq.gz",
        reverse_reads = "analysis/filteredDataDir/{sample}_R2_val_2.fq.gz",
        forwardUnpaired = "analysis/filteredDataDir/{sample}_R1_val_1_unpaired.fq.gz",
        reverseUnpaired = "analysis/filteredDataDir/{sample}_R2_val_2_unpaired.fq.gz"
    message: "trimming {wildcards.sample} reads"
    log:
        "analysis/logs/trimmomatic/{sample}.log"
    params :
        seedMisMatches =            str(config['trimmomatic']['seedMisMatches']),
        palindromeClipTreshold =    str(config['trimmomatic']['palindromeClipTreshold']),
        simpleClipThreshhold =      str(config['trimmomatic']['simpleClipThreshold']),
        LeadMinTrimQual =           str(config['trimmomatic']['LeadMinTrimQual']),
        TrailMinTrimQual =          str(config['trimmomatic']['TrailMinTrimQual']),
        windowSize =                str(config['trimmomatic']['windowSize']),
        avgMinQual =                str(config['trimmomatic']['avgMinQual']),
        minReadLen =                str(config['trimmomatic']['minReadLength']),
        phred = 		            str(config["trimmomatic"]["phred"])
    threads: 10
    conda:
        "../envs/trimmomatic_env.yaml"
    shell:
        "trimmomatic PE {params.phred} -threads {threads} "
        "{input.reads} "
        "{output.forward_reads} "
        "{output.forwardUnpaired} "
        "{output.reverse_reads} "
        "{output.reverseUnpaired} "
        "ILLUMINACLIP:{input.adapters}:{params.seedMisMatches}:{params.palindromeClipTreshold}:{params.simpleClipThreshhold} "
        "LEADING:{params.LeadMinTrimQual} "
        "TRAILING:{params.TrailMinTrimQual} "
        "SLIDINGWINDOW:{params.windowSize}:{params.avgMinQual} "
        "MINLEN:{params.minReadLen} &>{log}"

rule adaptNames:
	input:
	        R1="analysis/filteredDataDir/{sample}_R1_val_1.fq.gz",
                R2="analysis/filteredDataDir/{sample}_R2_val_2.fq.gz"
	output:
                R1="analysis/finalDataDir/{sample}_R1.fastq.gz",
                R2="analysis/finalDataDir/{sample}_R2.fastq.gz"
	shell:
		"""
		cp {input.R1} {output.R1}
		cp {input.R2} {output.R2}
		"""