#!/bin/bash

basedirectory="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/simulated-reads"

echo "indexing"
for directoryname in "$basedirectory"/*; do 
    if [ -d "$directoryname" ]; then 
        echo $directoryname
        species=${directoryname##*_}
        bamfilename=$directoryname"/"$species"_from_fa.bam"

        samtools index $bamfilename &
    fi
done
wait

echo "subsampling"
for directoryname in "$basedirectory"/*; do 
    if [ -d "$directoryname" ]; then 
        echo $directoryname
        species=${directoryname##*_}
        bamfilename=$directoryname"/"$species"_from_fa.bam"
        outbam=$directoryname"/"$species"_subsampled.bam"

        samtools view -b -s 0.1 $bamfilename > $outbam &
    fi
done
wait