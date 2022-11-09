def multimapping():
	if (config["multipleMapping"]=="yes"):
		return "--nonunique all"
	else:
		return ""

rule count_matrix:
	input:
		expand("analysis/mappedDataDir/{sample}/{sample}ReadsPerGene.out.tab", sample=samples)
	output:
		"analysis/counts/all2.tsv"
	shell:
		"awk 'NF > 1{ a[$1] = a[$1]\"\t\"$2} END {for( i in a ) print i a[i]}' {input} > {output}"

rule count_htseq_multi:
	input:
		bamFiles=expand("analysis/mappedDataDir/{sample}/{sample}Aligned.sortedByCoord.out.bam", sample=samples),
		annotation="{annotationDir}/annotation.gff".format(annotationDir=config["reference"]["annotation"])
	output:
		"analysis/counts/all.tsv"
	conda:
		"../envs/htseq.yaml"
	params:
		stranded=config["stranded"] ,multimap=multimapping(),header="Gene\\\t"+"\\\t".join(samples)
	shell:
		"""
		echo {params.header}>{output}
		
		htseq-count {params.multimap} -s {params.stranded} -f bam -r pos {input.bamFiles} {input.annotation} >> {output}
		"""


def get_contrast(wildcards):
	return config["contrasts"][wildcards.contrast]


rule diffExpresionGene:
	input:
		countsFile="analysis/counts/all.tsv",
		experiment=config["samples"]
	output:
		table="analysis/results/diffexp/{condition1}vs{condition2}.csv",
		plots="analysis/results/diffexp/{condition1}vs{condition2}.figures.pdf"
	params:
		contrast1="{condition1}",contrast2="{condition2}"
	log:
		"analysis/log/{condition1}vs{condition2}.diffexp.log"
	conda:
		"../envs/R.yaml"
	script:
		"../scripts/deseq2.R"
#For multiple files as input
#rule diffExpressionGeneMulti:
#	input:
#		countsFiles=expand("analysis/mappedDataDir/{sample}/{sample}Gene.out.tab", sample=samples)
#	output:
#		table="analysis/results/diffexp/{condition1}vs{condition2}.csv",
#		plots="analysis/results/diffexp/{condition1}vs{condition2}.figures.pdf"
#	params:
#		contrast1="{condition1}",contrast2="{condition2}"
#	log:
		#"analysis/log/{condition1}vs{condition2}.diffexp.log"
#	conda:
#		"../envs/R.yaml"
#	script:
#		"../scripts/deseq2multi.R"