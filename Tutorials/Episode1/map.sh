#!/bin/bash

fastqdir=resources/samples
mapdir=mapped
mkdir $mapdir

hisat2 -p 8 --dta -x index/chrX_tran -1 $fastqdir/ERR188044_chrX_1.fastq.gz -2 $fastqdir/ERR188044_chrX_2.fastq.gz -S $mapdir/ERR188044.sam
hisat2 -p 8 --dta -x index/chrX_tran -1 $fastqdir/ERR188104_chrX_1.fastq.gz -2 $fastqdir/ERR188104_chrX_2.fastq.gz -S $mapdir/ERR188104.sam
hisat2 -p 8 --dta -x index/chrX_tran -1 $fastqdir/ERR188234_chrX_1.fastq.gz -2 $fastqdir/ERR188234_chrX_2.fastq.gz -S $mapdir/ERR188234.sam
hisat2 -p 8 --dta -x index/chrX_tran -1 $fastqdir/ERR188245_chrX_1.fastq.gz -2 $fastqdir/ERR188245_chrX_2.fastq.gz -S $mapdir/ERR188245.sam
hisat2 -p 8 --dta -x index/chrX_tran -1 $fastqdir/ERR188257_chrX_1.fastq.gz -2 $fastqdir/ERR188257_chrX_2.fastq.gz -S $mapdir/ERR188257.sam
hisat2 -p 8 --dta -x index/chrX_tran -1 $fastqdir/ERR188273_chrX_1.fastq.gz -2 $fastqdir/ERR188273_chrX_2.fastq.gz -S $mapdir/ERR188273.sam
hisat2 -p 8 --dta -x index/chrX_tran -1 $fastqdir/ERR188337_chrX_1.fastq.gz -2 $fastqdir/ERR188337_chrX_2.fastq.gz -S $mapdir/ERR188337.sam
hisat2 -p 8 --dta -x index/chrX_tran -1 $fastqdir/ERR188383_chrX_1.fastq.gz -2 $fastqdir/ERR188383_chrX_2.fastq.gz -S $mapdir/ERR188383.sam
hisat2 -p 8 --dta -x index/chrX_tran -1 $fastqdir/ERR188401_chrX_1.fastq.gz -2 $fastqdir/ERR188401_chrX_2.fastq.gz -S $mapdir/ERR188401.sam
hisat2 -p 8 --dta -x index/chrX_tran -1 $fastqdir/ERR188428_chrX_1.fastq.gz -2 $fastqdir/ERR188428_chrX_2.fastq.gz -S $mapdir/ERR188428.sam
hisat2 -p 8 --dta -x index/chrX_tran -1 $fastqdir/ERR188454_chrX_1.fastq.gz -2 $fastqdir/ERR188454_chrX_2.fastq.gz -S $mapdir/ERR188454.sam
hisat2 -p 8 --dta -x index/chrX_tran -1 $fastqdir/ERR204916_chrX_1.fastq.gz -2 $fastqdir/ERR204916_chrX_2.fastq.gz -S $mapdir/ERR204916.sam
 

