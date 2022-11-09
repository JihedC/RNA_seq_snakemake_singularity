rule getIndex:
	output:
		touch("reference/centrifuge/{indexname}.DONE".format(indexname=config["centrifuge"]["index"]))
	conda:
		"../envs/centrifuge.yaml"
	params:
		url=config["centrifuge"]["url"],outputfile="reference/centrifuge/{indexname}.tar.gz".format(indexname=config["centrifuge"]["index"])
	shell:
		"""
		wget -q -O {params.outputfile} {params.url}
		tar -C reference/centrifuge -xvf {params.outputfile}
		"""

rule centrifuge_PE:
	input:
		indexDONE="reference/centrifuge/{indexname}.DONE".format(indexname=config["centrifuge"]["index"]),
		R1="analysis/finalDataDir/{sample}_R1.fastq.gz",
		R2="analysis/finalDataDir/{sample}_R2.fastq.gz"
	output:
		report="analysis/results/centrifuge/{sample}_report.txt",
		classification="analysis/results/centrifuge/{sample}_classification.txt"
	conda:
		"../envs/centrifuge.yaml"
	threads:
		config["threads"]
	params:
		index="reference/centrifuge/{indexname}".format(indexname=config["centrifuge"]["index"])
		
	shell:
		"centrifuge -p {threads} -x {params.index} -1 {input.R1} -2 {input.R2} --report-file {output.report} |cut -f 1,3 >  {output.classification}"


rule updateTaxonomy:
	output:
		"reference/centrifuge/taxonomy.tab"
	conda:
		"../envs/centrifuge.yaml"
	shell:
		"""
		ktUpdateTaxonomy.sh reference/centrifuge
		"""

rule plotKrona:
	input:
		data=expand("analysis/results/centrifuge/{sample}_classification.txt",sample=samples),
		taxonomy="reference/centrifuge/taxonomy.tab"
	output:
		"analysis/results/krona/report.html"
	log:
		"analysis/log/centrifuge/krona.log"
	conda:
		"../envs/centrifuge.yaml"
	shell:
		"""
		ktImportTaxonomy -tax reference/centrifuge -o {output} {input.data} >{log}
		"""
