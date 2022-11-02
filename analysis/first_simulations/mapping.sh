#!/bin/bash

basedirectory="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/simulated-reads"
sizedistribution=$basedirectory"/size_distribution.tsv"
genomefilename="/mnt/expressions/robin_warner/robin_genomes/human/human_genome.fna"

for directoryname in "$basedirectory"/*; do 
    if [ -d "$directoryname" ]; then 
        species=${directoryname##*_}
        queryfilename=$directoryname"/"$species"_d.fa.gz"
        outfile=$directoryname"/"$species"_to_human.sai"
        logfile=$directoryname"/mapping.log"

        bwa aln $genomefilename $queryfilename 1> $outfile 2> $logfile &
    fi
done 

wait

for directoryname in "$basedirectory"/*; do 
    if [ -d "$directoryname" ]; then 
        species=${directoryname##*_}
        queryfilename=$directoryname"/"$species"_d.fa.gz"
        outfile=$directoryname"/"$species"_to_human.sai"
        logfile=$directoryname"/conversion.log"
        bamoutfile=$directoryname"/"$species"_to_human.bam"

        bwa samse $genomefilename $outfile $queryfilename 2> $logfile | samtools view -Sb > $bamoutfile &
    fi
done 

