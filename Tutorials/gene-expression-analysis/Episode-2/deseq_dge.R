#Load libraries

library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)



#Set the working directory
setwd('/home/kobina/deseq2')

#Load the count data
count_data <- read.csv('count_matrix.csv',header=TRUE,row.names=1)
colnames(count_data)
head(count_data)

#Load the sample information
sample_info <- read.csv('design.csv',header=TRUE,row.names=1)
colnames(sample_info)
head(sample_info)



#Set factor levels
sample_info$Treatment <- factor(sample_info$Treatment)
sample_info$Sequencing <- factor(sample_info$Sequencing)



#Create a deseq object and import the count data and sample information.
dds <- DESeqDataSetFromMatrix(countData = count_data,colData = sample_info, design = ~Sequencing + Treatment)


#Set the reference for the Treatment factor
dds$Treatment <- factor(dds$Treatment, levels = c("untreated","treated"))
#dds$Treatment <- relevel(dds$Treatment, ref = "untreated")

#Filter the genes
keep <- rowSums(counts(dds)>10) >= min(table(sample_info$Treatment))
dds <- dds[keep,]




#Perform the statistical test(s) to identify differentially expressed genes
dds <- DESeq(dds)
deseq_result <- results(dds)
deseq_result



#Change DESeq Object to R object(dataframe)
deseq_result <- as.data.frame(deseq_result)
class(deseq_result)

head(deseq_result)

#Order the result table by increasing p value
deseq_result_ordered <- deseq_result[order(deseq_result$pvalue),]
head(deseq_result_ordered)



#Make some queries



#Is FBgn0003360 gene differentially expressed?
deseq_result["FBgn0003360",]




#Is the Pasilla gene (ps, FBgn0261552) downregulated by the RNAi treatment?
deseq_result["FBgn0261552",]



#Extract the most differentially expressed genes due to the treatment.
#Select genes with a significant change in gene expression (adjusted p-value below 0.05)
#And log2fold change <1 and >1


#Step 1: filter based on p adjusted value
filtered <- deseq_result %>% filter(deseq_result$padj < 0.05)


#Step 2: filter based on fold changes. here we will use a threshold of 1
filtered <- filtered %>% filter(abs(filtered$log2FoldChange) > 1 )


dim(deseq_result)
dim(filtered)


#Make queries 



#Save the deseq result. We will save the both the original data(res) and the filtered one(hits)
write.csv(deseq_result,'de_result.all.csv')
write.csv(filtered,'de_result.filtered.csv')


#Save the normalized counts 
normalized_counts <- counts(dds,normalized=TRUE)
head(normalized_counts)
write.csv(normalized_counts,'normalized_counts.csv')




#VISUALIZATION 




#Dispersion plot
plotDispEsts(dds) 


#Examine the counts of reads for a single gene across the groups.
#Lets plot the count of the top gene (gene which had the smallest p value from the results).
plotCounts(dds, gene=which.min(deseq_result$padj), intgroup = "Treatment")
plotCounts(dds, gene=which.min(deseq_result$padj), intgroup = "Sequencing")

sample_info




#PCA
#pca stands for principal component analysis.
#it is a dimensionality reduction technique and in gene expression analysis, 
#it can be used to explain the variance in gene expression datasets.
#to generate the pca plot we will first perform a variance stabilizing transformation. 
#we will use the vst function in deseq.
#after that we use the transformed values to plot the PCA




#variance stabilizing transformation
vsd <- vst(dds,blind=FALSE)

#use transformed values to generate a pca plot
plotPCA(vsd,intgroup=c("Sequencing","Treatment"))








#Heatmaps
#R Package:pheatmap


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





#Heatmap of log transformed normalized counts. We will use the top 10 genes.

#top 10 genes
top_hits=deseq_result[order(deseq_result$padj),][1:10,]
top_hits=row.names(top_hits)
top_hits
rld <- rlog(dds, blind=FALSE)



pheatmap(assay(rld)[top_hits,], cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=FALSE)
pheatmap(assay(rld)[top_hits,])

#add annotation
annot_info <- as.data.frame(colData(dds)[,c("Sequencing","Treatment")])

pheatmap(assay(rld)[top_hits,],cluster_rows=FALSE,show_rownames=TRUE,cluster_cols=FALSE,
         annotation_col=annot_info)


#Heatmap of Z scores. We will use the top 10 genes.
cal_z_score <- function(x) {(x-mean(x)) / sd(x)}

zscore_all <- t(apply(normalized_counts, 1, cal_z_score))
zscore_subset <- zscore_all[top_hits,]
pheatmap(zscore_subset)



#MA Plot

plotMA(dds,ylim=c(-2,2))

#remove the noise associated with low count genes
resLFC <- lfcShrink(dds,coef="Treatment_treated_vs_untreated",type="apeglm")

plotMA(resLFC,ylim=c(-2,2))




#Volcano Plot


#change resLFC to a dataframe
resLFC <- as.data.frame(resLFC)

#label the genes

resLFC$diffexpressed <- "NO"
resLFC$diffexpressed[resLFC$log2FoldChange>0.1 & resLFC$padj<0.05] <- "UP"
resLFC$diffexpressed[resLFC$log2FoldChange<0.1 & resLFC$padj<0.05] <- "DOWN"

resLFC$delabel<-NA

#generate the volcano plot
ggplot(data=resLFC,aes(x=log2FoldChange,y=-log10(pvalue),col=diffexpressed,label=delabel))+
  geom_point()+
  theme_minimal()+
  geom_text_repel()+
  scale_color_manual(values=c('blue','black','red'))+
  theme(text=element_text(size=20))





