


#Differential expression

#load the libraries




library(ggplot2)
library(ballgown)
library(genefilter)
library(RSkittleBrewer)
library(devtools)
library(dplyr)
library(ggrepel)
library(pheatmap)
library(gplots)
library(GenomicRanges)
library(viridis)
library(cowplot)


#lets load the sample information
pheno_data <- read.csv("resources/geuvadis_phenodata.csv")

#let's show information for first 6 samples
head(pheno_data)

#Load the expression data using ballgown
bg_chrX <- ballgown(dataDir="abundance",samplePattern="ERR",pData=pheno_data)


#Lets filter out transcripts with low variance
#This is done to remove some genes that have few counts. Filtering improves the statistical power of differential expression analysis. 
#We use variance filter to remove transcripts with low variance( 1 or less)

bg_chrX_filt<- subset(bg_chrX,"rowVars(texpr(bg_chrX))>1",genomesubset=TRUE)


#Let's test on transcripts
de_transcripts <- stattest(bg_chrX_filt,feature="transcript",covariate="sex",adjustvars=c("population"),getFC=TRUE,meas="FPKM")
# the results_transcripts does not contain identifiers. we will therefore add this information

#add identifiers
de_transcripts = data.frame(geneNames=ballgown::geneNames(bg_chrX_filt), geneIDs=ballgown::geneIDs(bg_chrX_filt), de_transcripts)

# Let's test on genes
de_genes <- stattest(bg_chrX_filt,feature="gene",covariate="sex",adjustvars = c("population"),getFC=TRUE, meas="FPKM")


#lets get the gene names
bg_filt_table=texpr(bg_chrX_filt,'all')
gene_names=unique(bg_filt_table[,9:10])
features=de_genes$id
mapped_gene_names=vector()
for (i in features) 
{  query=gene_names%>%filter(gene_id==i & gene_name != '.') ; n_hit=dim(query)[1]; if (n_hit==1) {mapped_gene_names=append(mapped_gene_names,query$gene_name[[1]]) } else
{mapped_gene_names=append(mapped_gene_names,'.') }    
}
#add the mapped gene names to the de genes table
de_genes$gene_name <- mapped_gene_names
de_genes <- de_genes[, c('feature','gene_name','id','fc','pval','qval')]



de_genes[,"log2fc"] <- log2(de_genes[,"fc"])
de_transcripts[,"log2fc"] <- log2(de_transcripts[,"fc"])


#Let's arrange the results from the smallest P value to the largest
de_transcripts = arrange(de_transcripts,pval)
de_genes = arrange(de_genes,pval)


#Let's subset transcripts that are detected as differentially expressed at qval <0.05
subset_transcripts <- subset(de_transcripts,de_transcripts$qval<0.05)

#do same for the genes
subset_genes <- subset(de_genes,de_genes$qval<0.05)

dir.create('plots')

print('generating plots')


#gene expression for a isoforms of gene XIST
png('plots/XIST.png')
myplot=plotTranscripts(ballgown::geneIDs(bg_chrX)[ballgown::geneNames(bg_chrX) == "XIST"], bg_chrX, main=c('Gene XIST in sample ERR188234'), sample=c('ERR188234'))
print(myplot)
dev.off()
#DONE






#volcano plot


#https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html

#mycolors <- c("red", "green", "black")
#names(mycolors) <- c("DOWN", "UP", "NO")
#p3 <- p2 + scale_colour_manual(values = mycolors)

#de_genes$diffexpressed[de_genes$log2fc < -0.6 & de_genes$qval < 0.05] <- "DOWN"
de_genes$diffexpressed <- "NO"
de_genes$diffexpressed[de_genes$log2fc > 1 & de_genes$pval < 0.05] <- "UP"
de_genes$diffexpressed[de_genes$log2fc < -1 & de_genes$pval < 0.05] <- "DOWN"
de_genes$delabel <- NA
de_genes$delabel[de_genes$diffexpressed != "NO"] <- de_genes$id[de_genes$diffexpressed != "NO"]

options(ggrepel.max.overlaps = Inf)

png('plots/volcano.png',width = 1800, height = 1000) #,width = 1800, height = 1000
volcano=ggplot(data=de_genes, aes(x=log2fc, y=-log10(pval), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.8, 0.8), col="red") +
  theme(text=element_text(size=20))
  
  #geom_hline(yintercept=-log10(0.05), col="red")
print(volcano)
dev.off()
#DONE


#MAPLOT
#https://davetang.org/muse/2017/10/25/getting-started-hisat-stringtie-ballgown/

png('plots/maplot.png',width = 1800, height = 1000)
de_transcripts$mean <- rowMeans(texpr(bg_chrX_filt))
maplot=ggplot(de_transcripts, aes(log2(mean), log2(fc), colour = qval<0.05)) +
  scale_color_manual(values=c("#999999", "#FF0000")) +
  geom_point() +
  theme(legend.text=element_text(size=20),legend.title=element_text(size=20)) +
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20)) +
  geom_hline(yintercept=0)
  



print(maplot)
dev.off()
#DONE


#heatmap
#Visualise the differentially expressed genes as a heatmap

#https://pydupont.github.io/BallGownTuto/Protocol.html



#get fpkm values 

png('plots/heatmap_clustered.png')

fpkm = gexpr(bg_chrX_filt)
fpkm = log2(fpkm+1)
g_ids=subset_genes$id
hits= which (de_genes$id %in% g_ids)
hit_frame=fpkm[hits,]
row.names(hit_frame) <- g_ids
heatmap_image=pheatmap(hit_frame)
print(heatmap_image)
dev.off()

png('plots/heatmap_unclustered.png')
heatmap_image=pheatmap(hit_frame,cluster_rows = F,cluster_cols=F)
print(heatmap_image)
dev.off()


print("plots have been generated and saved in 'plots directory' ")



#annotation_col <- data.frame(Group = pheno_data$sex)
#annotation_col[,1] <- gsub(pattern = "male", replacement = "Male", x = annotation_col[,1])
#annotation_col[,1] <- gsub(pattern = "female", replacement = "Female", x = annotation_col[,1])
#my_color=list(
#sample = c(male = "#5977ff", female = "#f74747"))





  
