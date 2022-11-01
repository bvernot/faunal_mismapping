# Faunal mismapping Speedrun Project

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

The script is currently in a really bad place, but I think you *can* use it if you really want to:

```
time python3 bin/mismapping.py \
     --bam $bam \
     --strategy random \
     --sites /mnt/expressions/benjamin_vernot/pendant_2022/sites.2.txt \
     --flush 1000000 > $outfile
```

This just generates a table of alleles at each site, and doesn't actually do the calculations we want to do - I plan to add that. Currently I then merge that data w/ the allele data found here:

```
/mnt/expressions/benjamin_vernot/soil_capture_2017/site_categories_for_capture/lineage_and_qc_sites_with_burden_alleles_cats.bed
```

I have some R code for that but ultimately I want the python script to do all of the work, so I won't put it here unless someone wants it.


### Doing basic popgen stuff with Poseiden
