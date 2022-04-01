#!/bin/bash

gtf=resources/chrX.gtf
assembly=assembly
mapdir=mapped
mkdir $assembly

echo "assemblying ERR188044"
stringtie $mapdir/ERR188044.bam -l ERR188044 -p 8 -G $gtf -o $assembly/ERR188044.gtf

echo "assemblying ERR188104"
stringtie $mapdir/ERR188104.bam -l ERR188104 -p 8 -G $gtf -o $assembly/ERR188104.gtf

echo "assemblying ERR188234"
stringtie $mapdir/ERR188234.bam -l ERR188234 -p 8 -G $gtf -o $assembly/ERR188234.gtf

echo "assemblying ERR188245"
stringtie $mapdir/ERR188245.bam -l ERR188245 -p 8 -G $gtf -o $assembly/ERR188245.gtf

echo "assemblying ERR188257"
stringtie $mapdir/ERR188257.bam -l ERR188257 -p 8 -G $gtf -o $assembly/ERR188257.gtf

echo "assemblying ERR188273"
stringtie $mapdir/ERR188273.bam -l ERR188273 -p 8 -G $gtf -o $assembly/ERR188273.gtf

echo "assemblying ERR188337"
stringtie $mapdir/ERR188337.bam -l ERR188337 -p 8 -G $gtf -o $assembly/ERR188337.gtf

echo "assemblying ERR188383"
stringtie $mapdir/ERR188383.bam -l ERR188383 -p 8 -G $gtf -o $assembly/ERR188383.gtf

echo "assemblying ERR188401"
stringtie $mapdir/ERR188401.bam -l ERR188401 -p 8 -G $gtf -o $assembly/ERR188401.gtf

echo "assemblying ERR188428"
stringtie $mapdir/ERR188428.bam -l ERR188428 -p 8 -G $gtf -o $assembly/ERR188428.gtf

echo "assemblying ERR188454"
stringtie $mapdir/ERR188454.bam -l ERR188454 -p 8 -G $gtf -o $assembly/ERR188454.gtf

echo "assemblying ERR204916"
stringtie $mapdir/ERR204916.bam -l ERR204916 -p 8 -G $gtf -o $assembly/ERR204916.gtf
