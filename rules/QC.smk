#############################################################################################################
#
#
#					Metrics module
#
#
#############################################################################################################
rule fastqc_dataDir:
	input:
		files=expand("analysis/mappedDataDir/{sample}/{sample}Aligned.sortedByCoord.out.bam",sample=samples)
	output:
		touch("analysis/QC/dataDir/fastqc.DONE")
	params:
		dir="analysis/QC/dataDir",
		dirInput="analysis/dataDir"
	conda:
		"envs/fastqc.yaml"
	shell:
		"fastqc -o {params.dir} {params.dirInput}/*gz"

rule fastqc_filteredDataDir:
	input:
		files=expand("analysis/mappedDataDir/{sample}/{sample}Aligned.sortedByCoord.out.bam",sample=samples)
	output:
		touch("analysis/QC/filteredDataDir/fastqc.DONE")
	params:
		dir="analysis/QC/filteredDataDir",
		dirInput="analysis/filteredDataDir"
	conda:
		"envs/fastqc.yaml"
	shell:
		"fastqc -o {params.dir} {params.dirInput}/*gz"

rule fastqc_finalDataDir:
	input:
		files=expand("analysis/mappedDataDir/{sample}/{sample}Aligned.sortedByCoord.out.bam",sample=samples)
	params:
		dir="analysis/QC/finalDataDir",
		dirInput="analysis/finalDataDir"
	output:
		touch("analysis/QC/finalDataDir/fastqc.DONE")
	conda:
		"envs/fastqc.yaml"
	shell:
		"fastqc -o {params.dir} {params.dirInput}/*gz"

rule picardMetrics:
	input:
		BAM="analysis/mappedDataDir/{sample}/{sample}Aligned.sortedByCoord.out.bam",
		reference="{referenceDir}/reference.fa".format(referenceDir=config["reference"]["index"])
	output:
		touch("analysis/QC/{sample}_stats.DONE")
	params:
		dir="analysis/QC/{sample}"
	conda:
		"envs/picard.yaml"
	shell:
		"""
		picard CollectMultipleMetrics VALIDATION_STRINGENCY=SILENT I={input.BAM} R={input.reference} O={params.dir}
		"""

rule idxstats:
	input:
		BAM="analysis/mappedDataDir/{sample}/{sample}Aligned.sortedByCoord.out.bam",
		BAI="analysis/mappedDataDir/{sample}/{sample}Aligned.sortedByCoord.out.bam.bai"
	output:
		"analysis/QC/{sample}.idxstats.txt"
	shell:
		"""
		samtools idxstats {input.BAM} >{output}
		"""

rule deduplicatedStats:
	input:
		BAM="analysis/mappedDataDir/{sample}/{sample}Aligned.sortedByCoord.out.bam"
	output:
		BAM="analysis/mappedDataDir/{sample}/{sample}.dedup.bam",
		STATS="analysis/QC/{sample}.dedup.stats"
	shell:
		"""
		picard MarkDuplicates VALIDATION_STRINGENCY=SILENT I={input} O={output.BAM} M={output.STATS}
		"""

rule createRefFlat:
	input:
		"{annotationDir}/annotation.gff".format(annotationDir=config["reference"]["annotation"])
	output: 
		"{referenceDir}/ref_flat.txt".format(referenceDir=config["reference"]["index"])
	shell:
		"""
		wget {input} -O {output}
		"""
rule createRibosomal:
	input:
		annotation="{annotationDir}/annotation.gff".format(annotationDir=config["reference"]["annotation"]),
		chrom_sizes="{referenceDir}/chrom_sizes.txt".format(referenceDir=config["reference"]["index"])
	output:
		"{referenceDir}/ribosomal.txt".format(referenceDir=config["reference"]["index"])
	shell:
		"""
				
		# Sequence names and lengths. (Must be tab-delimited.)
		perl -lane 'print "\@SQ\tSN:$F[0]\tLN:$F[1]\tAS:mm10"' $chrom_sizes |grep -v _ >> {output}

		# Intervals for rRNA transcripts.
		grep 'gene_type "rRNA"' $genes | awk '$3 == "transcript"' | cut -f1,4,5,7,9 | perl -lane '/transcript_id "([^"]+)"/ or die "no transcript_id on $.";
	        print join "\t", (@F[0,1,2,3], $1)' | sort -k1V -k2n -k3n
>> {output}
		"""

rule RNAmetrics:
	input:
		BAM="analysis/mappedDataDir/{sample}/{sample}Aligned.sortedByCoord.out.bam",
		REF_FLAT="{referenceDir}/ref_flat.txt".format(referenceDir=config["reference"]["index"]),
		RIBOSOMAL="{referenceDir}/ribosomal.txt".format(referenceDir=config["reference"]["index"])
	output:
		"analysis/QC/{sample}.RNAstats.txt"
	params:
		 STRAND="NONE"
	shell:
		"""
		rnaseqc {input} {output}

		"""

rule makereport:
	input:
		"analysis/QC/dataDir/fastqc.DONE",
		"analysis/QC/filteredDataDir/fastqc.DONE",
		"analysis/QC/finalDataDir/fastqc.DONE",
		expand("analysis/QC/{sample}_stats.DONE",sample=samples),
		expand("analysis/QC/{sample}.idxstats.txt",sample=samples),
		expand("analysis/QC/{sample}.dedup.stats",sample=samples)

	output:
		"analysis/results/report/multiqc_report.html"
	conda:
		"envs/multiqc.yaml"		
	shell:
		"multiqc -o analysis/results/report analysis"