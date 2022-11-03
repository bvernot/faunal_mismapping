#!/bin/bash

basedirectory="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/simulated-reads"
sizedistribution=$basedirectory"/size_distribution.tsv"
genomefilename="/mnt/expressions/robin_warner/robin_genomes/human/human_genome.fna"

echo "mapping"
for directoryname in "$basedirectory"/*; do 
    if [ -d "$directoryname" ]; then 
        echo $directoryname
        species=${directoryname##*_}
        queryfilename=$directoryname"/"$species"_d.fa.gz"
        outfile=$directoryname"/"$species"_to_human.sai"
        logfile=$directoryname"/mapping.log"

        bwa aln $genomefilename -n 0.01 -o 2 -l 16500 $queryfilename 1> $outfile 2> $logfile &
    fi
done 

wait

echo "converting to bam"
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
