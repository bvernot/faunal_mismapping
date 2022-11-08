#!/bin/bash

basedirectory="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/simulated-reads"
resultsdirectory=$basedirectory"/parsed-bams"


## we also did ..
# for file in *_*/*_to_human-bam2bam*subsampled.bam ; do samtools sort $file > ${file/.bam}-sorted.bam ; done

for directoryname in "$basedirectory"/*; do 
    if [ -d "$directoryname" ]; then 
        species=${directoryname##*_}
        bamfile=$directoryname"/"$species"_to_human-bam2bam_subsampled-sorted.bam"
        tsvfile=$basedirectory"/parsed-bams/"$species"_subsampled.tsv"

        if [ "$species" != "mismapping/data/simulated-reads/parsed-bams" ]; then
            echo $species
            python3 bam-parser.py $bamfile $tsvfile

        fi
    fi
done 
