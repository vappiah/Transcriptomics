
#Differential expression

#load the libraries

library(devtools)
library(dplyr)
library(ggplot2)
library(gplots)
library(GenomicRanges)
library(viridis)
library(ballgown)
library(RSikttleBrewer)
library(genefilter)

#lets load the sample information
pheno_data <- read.csv("resources/guevadis_phenodata.csv")

#let's show information for first 6 samples
head(pheno_data)

#Load the expression data using ballgown
bg_chrX <- ballgown(dataDir="ballgown",samplePattern="ERR",pData=pheno_data)

#Note:We specify three parameters:dataDir to specify where the estimates are stored, samplePattern to indicate a  name pattern that is present in the sample names and pData which specifies the phenotypic information (we loaded it as pheno_data).

#We can check some information about the loaded abundances
#check object type and attributes
class(bg_chrX)
attr("package")

#Let's do some queries on the ballgown object
bg_chrX

#Let's get what methods are available for the ballgown object.
methods(class="ballgown")

#We can display expression levels for genes, transcripts,exons and introns using the functions gexpr,texpr, eexpr(), and iexpr() respectively

#Let's get the gene and transcript expression levels
gexpr(bg_chrX)
texpr(bg_chrX

#Lets get some few rows
head(gexpr(bg_chrX),2)
head(texpr(bg_chrX),2)

#Lets filter out transcripts with low variance
#This is done to remove some genes that have few counts. Filtering improves the statistical power of differential expression analysis. 
#We use variance filter to remove transcripts with low variance( 1 or less)

bg_chrX_filt<- subset(bg_chrX,"rowVars(texpr(bg_chrX))>1",genomesubset=TRUE)

#How many transcripts were filtered out?

#Let's find out how many transcripts are now available
bg_chrX_filt

#We can now perform the differential expression analysis using the stattest() function. The parameters we use are, feature, covariate, adjustvars, getFC and meas. Because we are testing for genes and transcripts that are differentially expressed between males and females, sex becomes our covariate of interest. Transcripts might also be differentially expressed among differential populations but since this is not our group of interest it becomes a cofounder. We correct for that using the adjustvars parameter. We also want to get the confounder-adjusted fold change values so we set getFC to TRUE. Ballgown gives the expression values in different measures but we will use FPKM as the meas. Testing for differential expression at exons and introns can also be done by changing the feature parameter

#Let's test on transcripts
results_transcripts <- stattest(bg_chrX_filt,feature="transcript",covariate="sex",adjustvas=c("population"),getFC=TRUE,meas='FPKM")

# Let's test on genes
results_genes <- stattest(bg_chrX_filt,feature="gene",covariate="sex",adjustvars = c("population"),getFC=TRUE, meas="FPKM")


#lets check how many scripts are significantly expressed. we use adjusted p value (qval)
table(results_transcripts$qval < 0.05)

#we see that they are 13

# the results_transcripts does not contain identifiers. we will therefore add this information
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_chrX_filt), geneIDs=ballgown::geneIDs(bg_chrX_filt), results_transcripts)

#lets check if the identifiers have been added
head(results_transcripts)

#Now let's check which transcripts are detected and differentially expressed at qval < 0.05?
results_transcripts %>% filter(qval < 0.05)

#Let's arrange the results from the smallest P value to the largest
results_transcripts = arrange(results_transcripts,pval)
results_genes = arrange(results_genes,pval)

head(results_transcripts)
head(results_genes)

#We can now save the results to a file (csv recommended)
write.csv(results_transcripts, “fc_transcript_results.csv”, row.names=FALSE)
write.csv(results_genes, “fc_gene_results.csv”, row.names=FALSE)

#Let's subset transcripts that are detected as differentially expressed at qval <0.05
subset_transcripts <- subset(results_transcripts,results_transcripts$qval<0.05)

#do same for the genes
subset_genes <- subset(results_genes,results_genes$qval<0.05)

#Save the subsets to a file
write.csv(subset_transcripts, “fc_transcript_subset.csv”, row.names=FALSE)
write.csv(subset_genes, “fc_gene_subset.csv”, row.names=FALSE)















