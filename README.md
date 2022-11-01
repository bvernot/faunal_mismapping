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

### Doing basic popgen stuff with Poseiden
