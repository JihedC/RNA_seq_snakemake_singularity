############################################################################################################
#
#				DEXSEQ pipeline
#
##########################################################################################################

rule prepareDEXSEQannotation:
	input:
		"{annotationDir}/annotation.gtf".format(annotationDir=config["reference"]["annotation"])
	output:
		"analysis/annotation/DEXSEQ.gff"
	conda:
		"../envs/DEXSEQ.yaml"
	shell:
		"python2 scripts/dexseq_prepare_annotation.py {input} {output}"

rule createCountsPerGenePerExonTable:
        input:
                bam="analysis/mappedDataDir/{sample}/{sample}Aligned.sortedByCoord.out.bam",
                annotation="analysis/annotation/DEXSEQ.gff"
	params:
		stranded=config["stranded"]
	conda:
		"../envs/DEXSEQ.yaml"
	output:
		"analysis/counts/{sample}.countsPerGeneperExon.csv"
	shell:
		"python2 scripts/dexseq_count.py -s {params.stranded} -r pos -f bam {input.annotation} {input.bam} {output}"

rule diffExpressionExon:
	input:
		counts= expand("analysis/counts/{sample}.countsPerGeneperExon.csv",sample=samples),
		experiment=config["samples"],
		GFF="analysis/annotation/DEXSEQ.gff"
	params:
		contrast1="{condition1}",contrast2="{condition2}"
	output:
		table="analysis/results/DEXSEQ/{condition1}vs{condition2}.csv",
		figures="analysis/results/DEXSEQ/{condition1}vs{condition2}.figures.pdf"
	conda: 
		"../envs/R.yaml"
	script:
		"../scripts/DEXSEQ.R"