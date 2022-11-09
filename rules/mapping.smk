########################################################################################################
#
#					Pipeline para mapeo
#
#######################################################################################################

rule createGenomeIndex:
	input:
		fasta="{referenceDir}/reference.fa".format(referenceDir=config["reference"]["index"])
	output:
		touch("{indexDirectory}/index.DONE".format(indexDirectory=config["reference"]["index"]))
	threads:
		8
	params:
		mode="--runMode genomeGenerate",index=config["reference"]["index"]
	conda:
		"../envs/STAR.yaml"	
	shell:
		"STAR {params.mode} --runThreadN 8 --genomeDir {params.index} --genomeFastaFiles {input.fasta} --limitGenomeGenerateRAM 30000000000"

def is_single_end(wildcards):
	path="analysis/dataDir/{wildcards}/{wildcards}".format(wildcards=wildcards)
	READS,= glob_wildcards(path+"{READ}.fastq.gz")
	if len(READS)==1:
		return True
	else:
		return False

def getFastq(w):
	if not is_single_end("{sample}".format(sample=w.sample)):
		return expand("analysis/finalDataDir/{wildcards}_{READ}.fastq.gz",READ=["R1","R2"],wildcards=w.sample)
	else:
		return "analysis/finalDataDir/{wildcards}.fastq.gz".format(wildcards=w.sample)


rule mapping:
	input:
		fastq=getFastq,
		indexdone="{indexDirectory}/index.DONE".format(indexDirectory=config["reference"]["index"]),
		annotation= "{annotationDir}/annotation.gff".format(annotationDir=config["reference"]["annotation"])
	params:
		STAR="--outSAMtype BAM SortedByCoordinate --bamRemoveDuplicatesType UniqueIdentical --outWigType bedGraph --outWigNorm RPM --readFilesCommand zcat  --quantMode GeneCounts --twopassMode Basic",
		prefix="analysis/mappedDataDir/{sample}/{sample}",
		genome=config["reference"]["index"]
	threads:
		8
	log:
		"analysis/log/mapping/{sample}.log"
	conda:
		"../envs/STAR.yaml"
	output:
		"analysis/mappedDataDir/{sample}/{sample}Aligned.sortedByCoord.out.bam",
		"analysis/mappedDataDir/{sample}/{sample}ReadsPerGene.out.tab"
	shell:
		"STAR --runMode alignReads {params.STAR} --outFileNamePrefix {params.prefix} --runThreadN 8 --sjdbGTFfile {input.annotation} --genomeDir {params.genome} --readFilesIn {input.fastq} 2>{log}"

rule indexBAM:
	input:
		"analysis/mappedDataDir/{sample}/{sample}Aligned.sortedByCoord.out.bam"
	output:
		"analysis/mappedDataDir/{sample}/{sample}Aligned.sortedByCoord.out.bam.bai"
	shell:
		"samtools index {input}"


rule getNormalizedCoverage:
	input:
		BAM="analysis/mappedDataDir/{sample}/{sample}Aligned.sortedByCoord.out.bam",
		BAI="analysis/mappedDataDir/{sample}/{sample}Aligned.sortedByCoord.out.bam.bai"
	output:
		"analysis/toUCSC_norm/{sample}.normalized.bw"
	conda:
		"../envs/deeptools.yaml"
	params: norm="norm"
	shell:
		"bamCoverage --normalizeUsing RPKM --bam {input.BAM} -o {output}"



rule getCoverage:
	input:
		BAM="analysis/mappedDataDir/{sample}/{sample}Aligned.sortedByCoord.out.bam",
		BAI="analysis/mappedDataDir/{sample}/{sample}Aligned.sortedByCoord.out.bam.bai"
	output:
		"analysis/toUCSC/{sample}.bw"
	conda:
		"../envs/deeptools.yaml"
	shell:
		"bamCoverage  --bam {input.BAM} -o {output}"

rule publish:
	input:
		expand("analysis/toUCSC/{sample}.bw",sample=samples),expand("analysis/toUCSC_norm/{sample}.normalized.bw",sample=samples)
	output:
		touch("analysis/toUCSC/pasteToUCSC.txt")