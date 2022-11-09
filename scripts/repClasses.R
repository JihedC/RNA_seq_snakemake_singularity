library(DESeq2)
library(Glimma)
library("RColorBrewer")
library(pheatmap)
library(ggplot2)

file=snakemake@input[["counts"]]
fileExperiment=snakemake@input[["experiment"]]
fileTableRep=snakemake@input[["tableRepClasses"]]
annoFile=snakemake@input[["tableRepClasses"]]
control=snakemake@params[["cond1"]]
other=snakemake@params[["cond2"]]
outFile=snakemake@output[["table"]]
figuresFile=snakemake@output[["figures"]]
#It reads the table of region counts
table=read.table(file,header=TRUE)
table$rownames=table[,4]
#Read tehe experiment
designExp=read.table(fileExperiment,header=TRUE,row.names="sampleID")
designExp$condition=factor(designExp$condition)

#we summarize the data by repetitive classes
resumetable=aggregate(table[,7:(ncol(table)-1)],by=list(table$rownames),sum)
#clean environment
rm(table)
row.names(resumetable)=resumetable$Group.1
resumetable$Group.1=c()
names(resumetable)=row.names(designExp)
#Read tehe experiment
dds=DESeqDataSetFromMatrix(countData=resumetable,colData=designExp,design = ~condition)
dds$condition=relevel(dds$condition, ref=control)
dds=DESeq(dds)
res2=results(dds)

write.table(res2[order(res2$padj),],file=outFile,sep="\t")
pdf(figuresFile,paper="a4")
        
	plotMA(res2)
	rld=rlog(dds,blind=TRUE)
        plotPCA(rld,n=10000)

        sampleDists <- dist(t(assay(rld)))
        sampleDistMatrix <- as.matrix(sampleDists)
        rownames(sampleDistMatrix) <- rld$condition
        colnames(sampleDistMatrix) <- NULL
        colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

        pheatmap(sampleDistMatrix,
                clustering_distance_rows=sampleDists,
                clustering_distance_cols=sampleDists,
                col=colors)


dev.off()
dt.res2<- as.numeric(res2$padj<0.05)
annot=read.table(gzfile(annoFile),header=FALSE,sep="\t")
names(annot)=c("GeneID","repClass","repFamily","number")

glMDSPlot(dds, groups=designExp,path=dirname(outFile),folder=sub(".csv","",basename(outFile)),launch=FALSE)
glMDPlot(res2, status=dt.res2,path=dirname(outFile),folder=sub(".csv","",basename(outFile)),
        counts=counts(dds,normalized=TRUE),
        groups=dds$condition,
	anno=annot,
        transform=TRUE,
        samples=colnames(dds),
        launch=FALSE)


