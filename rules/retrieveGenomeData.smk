#################################################################################
#
#			Pipeline to extract data from SRA
#
################################################################################# 

rule download_genome:
	params: 
		fastaFile=config["GENOME_ZIP_FASTA_URL"]
	output:
		"{referenceDir}/reference.fa".format(referenceDir=config["reference"]["index"])
	shell:
		"curl {params.fastaFile} | gunzip -c > {output}"

rule download_gen_annotation:
	params:
		gffFile=config["GENOME_ZIP_GFF_URL"]
	output:
		"{annotationDir}/annotation.gff".format(annotationDir=config["reference"]["annotation"])
	shell:
		"curl {params.gffFile} | gunzip -c > {output}"




rule download_gen_annotation_gtf:
        params:
                gtfFile=config["GENOME_ZIP_GTF_URL"]
        output:
                "{annotationDir}/annotation.gtf".format(annotationDir=config["reference"]["annotation"])
        shell:
                "curl {params.gtfFile} | gunzip -c > {output}"