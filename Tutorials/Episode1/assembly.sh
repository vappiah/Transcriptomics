#!/bin/bash

gtf=resources/chrX.gtf
assembly=assembly
mapdir=mapped
mkdir $assembly


stringtie $mapdir/ERR188044.bam -l ERR188044 -p 8 -G $gtf -o $assembly/ERR188044.gtf
stringtie $mapdir/ERR188104.bam -l ERR188104 -p 8 -G $gtf -o $assembly/ERR188104.gtf
stringtie $mapdir/ERR188234.bam -l ERR188234 -p 8 -G $gtf -o $assembly/ERR188234.gtf
stringtie $mapdir/ERR188245.bam -l ERR188245 -p 8 -G $gtf -o $assembly/ERR188245.gtf
stringtie $mapdir/ERR188257.bam -l ERR188257 -p 8 -G $gtf -o $assembly/ERR188257.gtf
stringtie $mapdir/ERR188273.bam -l ERR188273 -p 8 -G $gtf -o $assembly/ERR188273.gtf
stringtie $mapdir/ERR188337.bam -l ERR188337 -p 8 -G $gtf -o $assembly/ERR188337.gtf
stringtie $mapdir/ERR188383.bam -l ERR188383 -p 8 -G $gtf -o $assembly/ERR188383.gtf
stringtie $mapdir/ERR188401.bam -l ERR188401 -p 8 -G $gtf -o $assembly/ERR188401.gtf
stringtie $mapdir/ERR188428.bam -l ERR188428 -p 8 -G $gtf -o $assembly/ERR188428.gtf
stringtie $mapdir/ERR188454.bam -l ERR188454 -p 8 -G $gtf -o $assembly/ERR188454.gtf
stringtie $mapdir/ERR204916.bam -l ERR204916 -p 8 -G $gtf -o $assembly/ERR204916.gtf
