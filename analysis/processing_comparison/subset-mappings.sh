# Copy the third allele mappings into the directory

# cp /mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/4_mappings/1_third/2_dog/dog-to-human_THIRD.bam /mnt/expressions/benjamin_vernot/faunal_mismapping/data/processing_comparisons/
# cp /mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/4_mappings/1_third/0_human/human-to-human_THIRD.bam /mnt/expressions/benjamin_vernot/faunal_mismapping/data/processing_comparisons/

datadirectory="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/processing_comparisons"

for bams in  0_human-to-human_REF-to-human_THIRD 2_dog-to-human_REF-to-human_THIRD dog-to-human_THIRD human-to-human_THIRD; do
    inputbam=$datadirectory/$bams".bam"
    # ls $inputbam
    output=$datadirectory"/"$bams".subset.bam"
    # echo $output
    derivedsites="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/sites/hominin_derived_sites.bed"

    bedtools intersect -a $inputbam -b $derivedsites > $output
done

