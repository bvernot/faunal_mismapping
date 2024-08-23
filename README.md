# Faunal mismapping "Speedrun" Project

Details in the Google Doc file.

Project directory on the servers:  
/mnt/expressions/benjamin_vernot/faunal_mismapping

Keep data here (but don't check it in!):  
/mnt/expressions/benjamin_vernot/faunal_mismapping/data


## Resources

### Simulating ancient DNA

You would have to change the number of chromosomes (it takes the first N chromosomes in the file - so for hg19 it would be good to put in 24 (1-22,X,Y), I think.

```
time python3 /mnt/expressions/benjamin_vernot/soil_capture_2017/ancient_animal_dna/fred_fake_seqdata/bin/chunk_genome.py \
     --num_seq $nseqs \
     --minlen 35 \
     --maxlen 100 \
     --substitution_file /mnt/scratch/fred/kraken/full_refseq/simulated_dataset/nonudg_error_forFred.txt \
     --outfile $ofile \
     --chromosomes 1 \
     --length /mnt/expressions/benjamin_vernot/soil_capture_2017/ancient_animal_dna/fred_fake_seqdata/A9180_final_sorted.deam53x3.bam.MQ25.lens.tsv \
     $fa
```

### Mapping w/ ancient parameters

### Calculating faunal mismapping proportions

This is still a work in progress.

For coverage based analyses, example command line runs:

```
time python3 bin/mismapping.py \
                 --bam data/test-mismapping-estimates/subsamples.subs_v18_human_REF_to_THIRD.1240k_props/dog_0.00.hum_1.00.200k.bam \
                 --sites data/sites/twist.1240k.burden.anno.in_apes.txt \
                 --minmq 25 \
                 -use-cats --categories mam_div_cat --strategy coverage --report-sim-truth dog human     
```

```
time python3 bin/mismapping.py \
                 --bam data/test_bams/A34692.pendant.uniq.L35MQ25.s0.1.bam \
                 --sites data/sites/twist.1240k.burden.anno.in_apes.txt \
                 --minmq 25 \
                 -use-cats --categories mam_div_cat --strategy coverage   
```


For hominin diagnostic position based analyses:

(still have to modify this sites file so that it has the mam_div_cat column, with appropriate categories - so the coverage estimates will not be exactly correct).
 
```
time python3 bin/mismapping.py \
                 --bam data/test_bams/A34692.pendant.uniq.L35MQ25.s0.1.bam \
                 --sites data/sites/hominin_derived_sites.fixcats.with_third.txt \
                 --minmq 25 \
                 -use-cats --categories b_3.mod --strategy coverage likelihood   
```

### Doing basic popgen stuff with Poseiden
