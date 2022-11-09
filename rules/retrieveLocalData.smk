rule retrieveLocalData:
	input:
		R1=lambda wildcards: tableSamples.loc[wildcards.sample]["R1"],
		R2=lambda wildcards: tableSamples.loc[wildcards.sample]["R2"]
	output:
		R1="analysis/dataDir/{sample}_R1.fastq.gz",R2="analysis/dataDir/{sample}_R2.fastq.gz"
	shell:
		"""
		cat {input.R1} > {output.R1}
		cat {input.R2} > {output.R2}
		"""
