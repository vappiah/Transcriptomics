
################################################################
#   Differential expression analysis with DESeq2


# load libraries
library(DESeq2)
library(pheatmap)
library(RColorBrewer)


# set working directory
setwd("/home/kobina/deseq2")

# read count data
count_table <- read.csv('raw_counts.tsv',sep='\t',row.names=1)
head(count_table)
dim(count_table)

# read the sample information
sample_info <- read.csv('design.tsv',sep='\t',row.names=1)
sample_info
dim(sample_info)

dim(count_table)

# set factor levels
factors <- factor(sample_info$Group)
groups <- unique(sample_info$Group)

groups

groups <- rev(groups)

groups 

sample_info$Group <- factors

sample_info$Group


# create DESeq object
dds <- DESeqDataSetFromMatrix(countData=count_table, colData=sample_info, design= ~Group)

# set the reference for the Group factor
dds$Group <- relevel(dds$Group, ref="control")

# filter out low count genes
# keep genes with at least N counts >= 10, where N = size of smallest group

keep <- rowSums(counts(dds) >=10)>= min(table(sample_info$Group))
dds <- dds[keep,]


# perform statistical tests
dds <- DESeq(dds,test="Wald",sfType='poscount')

# get the deseq result
deseq_result <- results(dds)

deseq_result

deseq_result <- as.data.frame(deseq_result)

class(deseq_result)

head(deseq_result)

dim(deseq_result)


names(deseq_result)


deseq_result$GeneName <- row.names(deseq_result)
names(deseq_result)
head(deseq_result)


deseq_result <- subset(deseq_result, 
                      select=c("GeneName","padj","pvalue","lfcSE","stat","log2FoldChange","baseMean"))

names(deseq_result)

write.table(deseq_result, file="deseq_result.all.tsv", row.names=F,sep="\t")


# extract de genes with padj < 0.05 and log2Foldchange <=-1 or >=1

deg <- subset(deseq_result, padj<0.05 & abs(log2FoldChange) >=1)

dim(deg)
dim(deseq_result)

deg <- deg[order(deg$padj),]

head(deg)

write.table(deg, file="deseq_deg.tsv", row.names=F, sep="\t")




################################################################
# Gene expression data visualization


# plot dispersion estimates
plotDispEsts(dds, main="GSE203159 Dispersion Estimates")

# create histogram plot of p-values
hist(deseq_result$padj, breaks=seq(0,1,length=21), col = "grey", border = "white",
     xlab="", ylab="", ylim=c(0,8000), main="GSE203159 Frequencies of padj-values")


# volcano plot
#set colors
old.pal <- palette(c("#00BFFF", "#FF3030"))

#set margin size
par(mar=c(4,4,2,1), cex.main=1.5)

#set title
title=paste(groups[1], "vs", groups[2])

#plot values
plot(deseq_result$log2FoldChange, -log10(deseq_result$padj), main=title,
     xlab="log2FC", ylab="-log10(Padj)", pch=20, cex=0.5)

with( subset(deseq_result, padj <0.05 & abs(log2FoldChange) >=1),
      points(log2FoldChange, -log10(padj), pch=20,  col=(sign(log2FoldChange) +3)/2, cex=1))


legend("bottomleft", title=paste("Padj<", 0.05, sep=""), 
       legend=c("down","up"), pch=20, col=1:2)






# variance stabilizing transformation
vsd <- vst(dds,blind=FALSE)



# pca plot
# use transformed values to generate a pca plot
plotPCA(vsd,intgroup=c("Group"))

Heatmaps
#R Package:pheatmap

#Heatmap of log transformed normalized counts. We will use the top 10 genes.


#top 10 genes
normalized_counts <- counts(dds,normalized=T)
transformed_counts <- log2(normalized_counts+1)
top_hits=row.names(deg[1:10,])
top_hits
top_hits <-transformed_counts[top_hits,]
pheatmap(top_hits, cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=FALSE)

#add annotation
annot_info <- as.data.frame(colData(dds) [c('Group')])
annot_info
pheatmap(top_hits,cluster_rows=FALSE,show_rownames=TRUE,cluster_cols=FALSE,annotation_col=annot_info)
pheatmap(top_hits,cluster_rows=TRUE,show_rownames=TRUE,cluster_cols=TRUE,annotation_col=annot_info)




#Heatmap of sample-to-sample distance matrix (with clustering) based on the normalized counts.

#generate the distance matrix
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix)

#set a color scheme
colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
#generate the heatmap
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists, col=colors)


################################################################
#   General expression data visualization

dat <- log10(counts(dds, normalized = T) + 1) # extract normalized counts

# box-and-whisker plot
lbl <- "log10(raw counts + 1)"

ss_info <- read.csv('design.tsv',sep='\t')
ss <- ss_info[order(ss_info$Group),]$SampleID

group_ <- subset(ss_info, ss_info$Group==groups[1])$SampleID
group_ <- append(group_,subset(ss_info, ss_info$Group==groups[2])$SampleID)

status <- as.integer(allgroups==groups[1])+1

ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
boxplot(dat[,group_], boxwex=0.6, notch=T, main="GSE203159", ylab="lg(norm.counts)", outline=F, las=2, col=status)
legend("topleft", groups, fill=palette(), bty="n")


# MD plot
par(mar=c(4,4,2,1), cex.main=1.5)
plot(log10(deseq_result$baseMean), deseq_result$log2FoldChange, main=paste(groups[1], "vs", groups[2]),
     xlab="log10(mean of normalized counts)", ylab="log2FoldChange", pch=20, cex=0.5)
with(subset(deseq_result, padj<0.05 & abs(log2FoldChange) >= 1),
     points(log10(baseMean), log2FoldChange, pch=20, col=(sign(log2FoldChange) + 3)/2, cex=1))
legend("bottomleft", title=paste("Padj<", 0.05, sep=""), legend=c("down", "up"), pch=20,col=1:2)
abline(h=0)
old.pal <- palette(c("#00BFFF", "#FF3030"))
palette(old.pal) # restore palette


# UMAP plot (multi-dimensional scaling)
library(umap)
dat <- dat [!duplicated(dat ), ] # first remove duplicates
par(mar=c(3,3,2,6), xpd=TRUE, cex.main=1.5)
ump <- umap(t(dat), n_neighbors = 3, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=5", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.25,0), legend=groups, pch=20,
       col=1:length(groups), title="Group", pt.cex=1.5)


#The codes in this section were picked directly from the deseq2 tutorial
#http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
par(mfrow=c(2,2),mar=c(2,2,1,1))
ylim <- c(-2.5,2.5)
resGA <- results(dds, lfcThreshold=.5, altHypothesis="greaterAbs")
resLA <- results(dds, lfcThreshold=.5, altHypothesis="lessAbs")
resG <- results(dds, lfcThreshold=.5, altHypothesis="greater")
resL <- results(dds, lfcThreshold=.5, altHypothesis="less")
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()
plotMA(resLA, ylim=ylim); drawLines()
plotMA(resG, ylim=ylim); drawLines()
plotMA(resL, ylim=ylim); drawLines()



