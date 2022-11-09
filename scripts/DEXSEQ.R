library(DEXSeq)
GFF=snakemake@input[["GFF"]]
outFile=snakemake@output[["table"]]
outFiguresFile=snakemake@output[["figures"]]
countFiles = snakemake@input[["counts"]]
sampleTableFile = snakemake@input[["experiment"]]


sampleTable=read.table(sampleTableFile,header=TRUE,sep="\t")
row.names(sampleTable)=sampleTable$sampleID
dxd = DEXSeqDataSetFromHTSeq(
	countFiles,
	sampleData=sampleTable,
	design= ~ sampleID + exon + condition:exon,
	flattenedfile=GFF )

dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd)
dxd = testForDEU( dxd)
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")
dxr = DEXSeqResults( dxd )
#table( dxr2$padj < 0.1 )

#print(table(now = dxr2$padj < 0.1 ))

#dxr = DEXSeq(dxd)

#Plots
pdf(outFiguresFile)
	plotDispEsts(dxd)
	plotMA(dxr,cex=0.8)
dev.off()

write.csv(dxr[!is.na(dxr$padj) & dxr$padj<0.05,],file=outFile)
DEXSeqHTML(dxr,FDR=0.05,color=c("#FF000080", "#0000FF80"),path=dirname(outFile))
