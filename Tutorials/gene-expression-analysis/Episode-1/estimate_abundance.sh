#!/bin/bash

abundancedir=abundance
mapdir=mapped

echo "Estimating abundance"

stringtie -e -B -p 8 -G stringtie_merged.gtf -o $abundancedir/ERR188044/ERR188044_chrX.gtf $mapdir/ERR188044.bam

stringtie -e -B -p 8 -G stringtie_merged.gtf -o $abundancedir/ERR188104/ERR188104_chrX.gtf $mapdir/ERR188104.bam

stringtie -e -B -p 8 -G stringtie_merged.gtf -o $abundancedir/ERR188234/ERR188234_chrX.gtf $mapdir/ERR188234.bam

stringtie -e -B -p 8 -G stringtie_merged.gtf -o $abundancedir/ERR188245/ERR188245_chrX.gtf $mapdir/ERR188245.bam

stringtie -e -B -p 8 -G stringtie_merged.gtf -o $abundancedir/ERR188257/ERR188257_chrX.gtf $mapdir/ERR188257.bam

stringtie -e -B -p 8 -G stringtie_merged.gtf -o $abundancedir/ERR188273/ERR188273_chrX.gtf $mapdir/ERR188273.bam

stringtie -e -B -p 8 -G stringtie_merged.gtf -o $abundancedir/ERR188337/ERR188337_chrX.gtf $mapdir/ERR188337.bam

stringtie -e -B -p 8 -G stringtie_merged.gtf -o $abundancedir/ERR188383/ERR188383_chrX.gtf $mapdir/ERR188383.bam

stringtie -e -B -p 8 -G stringtie_merged.gtf -o $abundancedir/ERR188401/ERR188401_chrX.gtf $mapdir/ERR188401.bam

stringtie -e -B -p 8 -G stringtie_merged.gtf -o $abundancedir/ERR188428/ERR188428_chrX.gtf $mapdir/ERR188428.bam

stringtie -e -B -p 8 -G stringtie_merged.gtf -o $abundancedir/ERR188454/ERR188454_chrX.gtf $mapdir/ERR188454.bam

stringtie -e -B -p 8 -G stringtie_merged.gtf -o $abundancedir/ERR204916/ERR204916_chrX.gtf $mapdir/ERR204916.bam

echo "Job completed"
