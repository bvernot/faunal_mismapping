#!/bin/bash
mergedBamDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/6_merged-bams"
resultsDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/7_subsets-vcfs-plink"
kSites="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/sites/all_1240k_sites.bed"

for speciesDir in $mergedBamDir"/"*; do
    species=`echo $speciesDir | awk -F "/" '{print $NF}'`
    echo ">" $species

    mkdir -p $resultsDir"/"$species

    for inFile in $speciesDir"/"*; do
        number=`echo $inFile | awk -F "_" '{print $NF}' | awk -F "." '{print $1}'`

        echo "  # "$number

        # subset bam
        subsetBam=$resultsDir"/"$species"/"$species"-"$number"_1240k.bam"
        bedtools intersect -a $inFile -b $kSites > $subsetBam

        # index bam
        samtools index $subsetBam

        # create vcf
        sampleName=$species"-"$number"_1240k"
        /home/mateja_hajdinjak/src/bam-caller/bam-caller.py --bam $subsetBam --strategy random --minmq 25 --sample-name $sampleName --output $resultsDir"/"$species"/"$sampleName

        # create plink
        plink --vcf $resultsDir"/"$species"/"$sampleName".vcf" --recode --make-bed --out $resultsDir"/"$species"/"$sampleName

    done
done
