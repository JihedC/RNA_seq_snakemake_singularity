###########################################################################
#
#			Repetitive regions pipeline
#
###########################################################################

rule createCountsPerRepetitiveRegions:
	input:
		bamFiles=expand("analysis/mappedDataDir/{sample}/{sample}Aligned.sortedByCoord.out.bam",sample=samples),
		baiFiles=expand("analysis/mappedDataDir/{sample}/{sample}Aligned.sortedByCoord.out.bam.bai",sample=samples),
		annotation=config["REPCLASSES_ZIP_bed_URL"]
	output:
		"analysis/counts/all.countsPerRepetitiveRegions.csv"
	params:
		header="chr\\\tstart\\\tend\\\tID\\\t\\\tsize\\\tstrand\\\t"+"\\\t".join(samples)
	singularity:'docker://biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'	
	shell:
		"""
		echo {params.header}>{output}
		bedtools multicov -bams {input.bamFiles} -bed {input.annotation}>> {output}
		"""

rule diffExpressionRepetitiveClasses:
	input:
		counts="analysis/counts/all.countsPerRepetitiveRegions.csv",
		experiment=config["samples"],
		tableRepClasses=config["REPCLASSES_ZIP_TABLE_URL"]
	params:
		cond1="{condition1}",cond2="{condition2}"
	conda:
		"../envs/R.yaml"
	output:
		table="analysis/results/repClasses/{condition1}vs{condition2}.repClasses.csv",
		figures="analysis/results/repClasses/{condition1}vs{condition2}.figures.repClasses.pdf" 
	script:
		"../scripts/repClasses.R"