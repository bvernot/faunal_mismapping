#!/bin/bash

basedirectory="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/simulated-reads"
sizedistribution=$basedirectory"/size_distribution.tsv"

echo "mapping"
for directoryname in "$basedirectory"/*; do 
    if [ -d "$directoryname" ]; then 
        echo $directoryname
        species=${directoryname##*_}
        queryfilename=$directoryname"/subsampled_seqs.fna"
        outfile=$directoryname"/"$species"_to_self.sai"
        logfile=$directoryname"/mapping_to_self.log"
        genomefilename=$directoryname"/endo/"$species"_genome.fna"
        echo $genomefilename

        bwa aln $genomefilename -n 0.01 -o 2 -l 16500 $queryfilename 1> $outfile 2> $logfile &
    fi
done

wait

echo "converting to bam"
for directoryname in "$basedirectory"/*; do 
    if [ -d "$directoryname" ]; then 
        species=${directoryname##*_}
        queryfilename=$directoryname"/subsampled_seqs.fna"
        outfile=$directoryname"/"$species"_to_self.sai"
        logfile=$directoryname"/conversion_to_self.log"
        bamoutfile=$directoryname"/"$species"_to_self.bam"
        genomefilename=$directoryname"/endo/"$species"_genome.fna"

        bwa samse $genomefilename $outfile $queryfilename 2> $logfile | samtools view -Sb > $bamoutfile &  
    fi
done 
