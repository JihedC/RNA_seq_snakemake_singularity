library(ggplot2)
library(DESeq2)
library("RColorBrewer")
library(pheatmap)
library("biomaRt")
###############################################################################################################################
#
#                                             FUNCTIONS
#
###############################################################################################################################
volcanoplot2 <- function (res) {
  # Make a basic volcano plot
  with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot",xlim=c(-3,3)))
  
  # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
  with(subset(res, padj<.05 ), points(log2FoldChange, -log10(padj), pch=20, col="red"))
  
}

###############################################################################################################################
#
#                                             MAIN
#
###############################################################################################################################

files=snakemake@input[["countsFiles"]]
fileExperiment=snakemake@input[["experiment"]]
#contrast=snakemake@params[["contrasts"]]
plotsFile=snakemake@output[["plots"]]
#TODO hacerlo para varios contrastes
control=snakemake@params[["contrast1"]]
print(control)
other=snakemake@params[["contrast2"]]

diffExpTable=snakemake@output[["table"]]
#diffExpNormTable=snakemake@output$normTable

###############################################################################################################################
#
#                                             DESEQ2 ANALYSIS
#
###############################################################################################################################


#Read the experiment
designExp=read.table(fileExperiment,header=TRUE,row.names="sampleID")
designExp$fileName=files

designExp$condition=factor(designExp$condition)
dds=DESeqDataSetFromHTSeqCount(countData=countTable,colData=designExp,design = ~condition)
dds$condition=relevel(dds$condition, ref=control)
dds=DESeq(dds)
res2=results(dds,contrast=c("condition",other,control))
rld=rlog(dds,blind=TRUE)
###############################################################################################################################
#
#                                             Download annotation
#
 ###############################################################################################################################

pdf(plotsFile,paper="a4")
	plotMA(res2)

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


###############################################################################################################################
#
#                                             Download annotation
#
###############################################################################################################################

library('biomaRt')


#Download annotation

genes <- row.names(countTable)
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))

G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",'mgi_symbol', 'chromosome_name',
                                                          'start_position', 'end_position', 'strand',
                                                          "description"),values=genes,mart= mart)
annotation=aggregate(G_list[,"mgi_symbol"],by=list(G_list$ensembl_gene_id,G_list$chromosome_name,
                                                    G_list$start_position,G_list$end_position,
                                                    G_list$strand,G_list$description),function(x){paste(x,sep="_")})
annotation=aggregate(mgi_symbol~ensembl_gene_id+chromosome_name+start_position+end_position+strand+description,data=G_list,FUN=paste,collapse=",")
row.names(annotation)=annotation$ensembl_gene_id
annotation$name=unlist(annotation$mgi_symbol)


 ###############################################################################################################################
 #
 #                                             Create reports and tables
 #
 ###############################################################################################################################
library(Glimma)

annotation=annotation[row.names(res2),]
annotation$GeneID=annotation$ensembl_gene_id
dt.res2<- as.numeric(res2$padj<0.05)
finalTable=merge(annotation,as.data.frame(res2),by.x="row.names",by.y="row.names")

row.names(finalTable)=finalTable$Row.names

finalTable=merge(finalTable,table,by.x="row.names",by.y="row.names")
finalTableNorm=merge(finalTable,counts(dds,normalized=TRUE),by.x="row.names",by.y="row.names")

write.table(finalTable,file=diffExpTable,sep="\t",quote=FALSE)
#write.table(finalTable,file=diffExpNormTable,sep="\t",quote=FALSE)

glMDSPlot(dds, groups=designExp,path=dirname(diffExpTable),folder=sub(".csv","",basename(diffExpTable)),launch=FALSE)
glMDPlot(res2, status=dt.res2,path=dirname(diffExpTable),folder=sub(".csv","",basename(diffExpTable)), 
	counts=counts(dds,normalized=TRUE), 
	groups=dds$condition, 
	transform=TRUE,
	samples=colnames(dds), 
	anno=annotation[,c(9,7,6)],
	launch=FALSE)


