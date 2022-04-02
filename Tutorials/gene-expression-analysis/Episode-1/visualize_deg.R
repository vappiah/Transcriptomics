


#download data
wget ftp://ftp.ccb.jhu.edu/pub/RNAseq_protocol/chrX_data.tar.gz

#downloaded data is compressed. So we have to decompress/extract the data from it
tar xvfz chrX_data.tar.gz

#lets check the contents of the extracted data
ls

#check all content
tree chrX_data

#prepare the data for analysis by executing the prepare_data.sh script
./prepare_data.sh chrX_data

#a new directory called resources will be created.
#check the content
ls resources

#lets check the content of the phenotype file. Normally this file will have to be prepared by you.
cat resources/geuvadis_phenodata.csv

#from the output there are 12 samples
ls resources/samples


#Lets start the analysis
#Step 1: Mapping
#The purpose of mapping is to infer which transcripts are expressed by identifying locations where a short read best matches the reference
#Mapping is performed using hisat2. The first step will be to index the reference genome and this will be followed by the mapping.

#lets make a directory to store the indexed genome. we also move into the directory

mkdir index

#splice-site and exon information.
extract_splice_sites.py resources/chrX.gtf >index/chrX.ss
extract_exons.py resources/chrX.gtf >index/chrX.exon

ls index

#display some lines of the generated files
head -n 5 chrX.ss
head -n 5 chrX.exon

#build the index using 8 threads. but first move into the index directory
cd index
hisat2-build -p 8 --ss chrX.ss --exon chrX.exon ../resources/chrX.fa chrX_tran

#cd to previous directory
cd ../

#Lets map the samples to the reference genome. Because we have multiple samples, it will be good
#to have a script that will automate the process for us. The map.sh script will do the job for us

#we will use 8 threads
#
./map.sh
 

#sort the generated sam files and convert to bam files.  This is done using samtools
#sorting arrange the data by chromosomes,contigs/scaffolds
#this enables the efficient access of the data
#the sorted data is saved as bamfile because bamfiles occupy less storage space
#you can execute the script (sort.sh)  instead of running the sort commands one by one.
./sort.sh   or bash sort.sh

#Lets remove the sam files
rm mapped/*.sam


#We now proceed to perform assembly. The assembly process involves the reconstruction of the mapped reads to enable accurate quantification of genes and their respective isoforms.
#assembly is performed using

# the assembly commands have been aggregated into a script called assembly.sh
./assembly.sh


#lets check the output directory

ls assembly


#lets check the content of one of the files
head assembly/ERR188044.gtf


#the next step is to merge all the files generated from the assembly process. The files contain transcripts and their respective quantities for each sample.

#place all file names in a single text file.

ls assembly/*.gtf > mergelist.txt

#lets confirm the operation was successful
cat mergelist.txt

#we now merge the transcripts using stringtie
stringtie --merge -p 8 -G resources/chrX.gtf -o stringtie_merged.gtf mergelist.txt

#lets read the content of the generated file (stringtie_merged.gtf)
head stringtie_merged.gtf
head -n 5 stringtie_merged.gtf


#lets find how many transcripts were identified
cat stringtie_merged.gtf |grep -v "^#" |awk '$3=="transcript" {print}' |wc -l

#lets compare the list of genes and transcripts to the known transcripts which are present in the reference annotation file.
gffcompare -r resources/chrX.gtf -G -o merged stringtie_merged.gtf

#open the output file
cat merged.stats


#Lets now estimate the abundance of the transcripts
./estimate_abundance.sh

