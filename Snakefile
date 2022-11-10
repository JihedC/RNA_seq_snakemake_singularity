# Snakemake pipeline for the analysis of RNA-seq Illumina data


################## Import libraries ##################


import pandas as pd
import os
import sys
from subprocess import call
import itertools
from snakemake.utils import R

################## Configuration file ##################

configfile: "experiment_params.yaml"


tableSamples = pd.read_table(config["samples"],sep="\t",index_col=0)
samples=list(set(tableSamples.index))
samples=sorted(samples)
# If there are R1/R2 columns it is not necessary dowload data from GEO
if config["preprocess"]=="no":
    dataDir="analysis/finalDataDir"
else:
    dataDir="analysis/dataDir"

###########################################################################################################################################################
#
#						Loading data from SRA or local files
#
###########################################################################################################################################################

if 'R1' in tableSamples:
    include: "rules/retrieveLocalData.smk"
    #There is no fastq files
else:
    include: "rules/retrieveData.smk"


include: "rules/retrieveGenomeData.smk"
include: "rules/preprocess.smk"
##########################################################################################################################################################
#
#						Get final file names
#
##########################################################################################################################################################

conditions = tableSamples.condition
conditionsUnique=set(conditions)

def getComparisonNames(conditions):
    x=list()
    conditionList=list(itertools.combinations(conditionsUnique,2))
    for elem1,elem2 in conditionList:
        if config["deseq2"]=="yes":
            x.append("analysis/results/diffexp/{condition1}vs{condition2}.csv".format(condition1=elem1,condition2=elem2))
            x.append("analysis/results/diffexp/{condition1}vs{condition2}.figures.pdf".format(condition1=elem1,condition2=elem2))
        if config["DEXEQ"]=="yes":
            x.append("analysis/results/DEXSEQ/{condition1}vs{condition2}.csv".format(condition1=elem1,condition2=elem2))
            x.append("analysis/results/DEXSEQ/{condition1}vs{condition2}.figures.pdf".format(condition1=elem1,condition2=elem2))
        if config["repetitiveAnalysis"]=="yes":
            x.append("analysis/results/repClasses/{condition1}vs{condition2}.repClasses.csv".format(condition1=elem1,condition2=elem2))
            x.append("analysis/results/repClasses/{condition1}vs{condition2}.figures.repClasses.pdf".format(condition1=elem1,condition2=elem2))
        if config["QC"]=="yes":
            x.append("analysis/results/report/multiqc_report.html")
    return x
#########################################################################################################################################################
#
#							Snakemake rules
#
#########################################################################################################################################################

rule all:
        input:
                getComparisonNames(conditions),"analysis/toUCSC/pasteToUCSC.txt","analysis/results/krona/report.html"

def getSample(wildcards):
        return samples


include: "rules/mapping.smk"
include: "rules/deseq2.smk"
include: "rules/DEXSEQ.smk"
include: "rules/diffExpRepClasses.smk"
include: "rules/QC.smk"
include: "rules/centrifuge.smk"