#!/bin/bash
#WARNING: This script is for education purposes. You modify this script at your own risk

#we are interested in the raw-data(samples) and reference genome(fasta and gff) and phenotype information. 
#So we will put them in a directory called resources.All other files will be removed

resourcedir=resources
mkdir $resourcedir

mv chrX_data/genes/chrX.gtf $resourcedir
mv chrX_data/genome/chrX.fa $resourcedir
mv chrX_data/geuvadis_phenodata.csv $resourcedir
mv -v chrX_data/samples $resourcedir


lets remove the unwanted files
rm -fr chrX_data

echo "data preparating is completed"
echo "all data can be found in the directory : \"resources\" " 
