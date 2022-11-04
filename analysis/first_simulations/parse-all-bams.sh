#!/bin/bash

basedirectory="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/simulated-reads"
resultsdirectory=$basedirectory"/parsed-bams"

for directoryname in "$basedirectory"/*; do 
    if [ -d "$directoryname" ]; then 
        species=${directoryname##*_}
        bamfile=$directoryname"/"$species"_to_human-bam2bam_subsampled.bam"
        tsvfile=$basedirectory"/parsed-bams/"$species"_subsampled.tsv"

        if [ "$species" != "mismapping/data/simulated-reads/parsed-bams" ]; then
            python3 bam-parser.py $bamfile $tsvfile

        fi
    fi
done 
