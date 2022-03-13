#!/bin/bash

mapdir=mapped

samtools sort -@ 8 -o $mapdir/ERR188044.bam $mapdir/ERR188044.sam
samtools sort -@ 8 -o $mapdir/ERR188104.bam $mapdir/ERR188104.sam
samtools sort -@ 8 -o $mapdir/ERR188234.bam $mapdir/ERR188234.sam
samtools sort -@ 8 -o $mapdir/ERR188245.bam $mapdir/ERR188245.sam
samtools sort -@ 8 -o $mapdir/ERR188257.bam $mapdir/ERR188257.sam
samtools sort -@ 8 -o $mapdir/ERR188273.bam $mapdir/ERR188273.sam
samtools sort -@ 8 -o $mapdir/ERR188337.bam $mapdir/ERR188337.sam
samtools sort -@ 8 -o $mapdir/ERR188383.bam $mapdir/ERR188383.sam
samtools sort -@ 8 -o $mapdir/ERR188401.bam $mapdir/ERR188401.sam
samtools sort -@ 8 -o $mapdir/ERR188428.bam $mapdir/ERR188428.sam
samtools sort -@ 8 -o $mapdir/ERR188454.bam $mapdir/ERR188454.sam
samtools sort -@ 8 -o $mapdir/ERR204916.bam $mapdir/ERR204916.sam
 

