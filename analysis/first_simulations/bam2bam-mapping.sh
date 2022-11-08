#!/bin/bash

basedirectory="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/simulated-reads"
sizedistribution=$basedirectory"/size_distribution.tsv"
genomefilename="/mnt/expressions/robin_warner/robin_genomes/human/human_genome.fna"

for directoryname in "$basedirectory"/*; do
    if [ -d "$directoryname" ]; then
        echo $directoryname
        species=${directoryname##*_}
        queryfilename=$directoryname"/"$species"_subsampled.bam"
        outfile=$directoryname"/"$species"_to_human-bam2bam_subsampled.bam"

        bwa bam2bam --only-aligned -g $genomefilename -n 0.01 -o 2 -l 16500 -t 4 -f $outfile $queryfilename &
    fi
done
wait
