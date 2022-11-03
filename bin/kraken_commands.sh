
rule generic_kraken:
    input: fa = "{base}.fa.gz"
    output: kraken = "{base}.kraken", # was previously temp
    threads: 1
    shell: """
    
     if [ $(hostname) != "bionc13" ] ; then 
         echo KRAKEN MUST BE RUN ON BIONC13. SLEEPING 30 SECONDS.
         sleep 30
         stop
     fi

    echo THIS FAILS IF THE KRAKEN DB IS NOT COPIED - fix it

     nproc=$(time ~frederic_romagne/kraken/install/kraken --threads {threads} --db /mnt/ramdisk/refseqReleaseKraken \
     --output {output.kraken} {input.fa} 2>&1 | grep 'processed' | cut -f1 -d' ')
     echo "count seqs {input.fa}"
     nseq=$(gunzip -c {input.fa} | grep -c -e '>' -e '@') || echo "no seqs found $nseq"

     if [ ! $nseq -eq $nproc ] ; then echo KRAKEN RUN FAILED; stop; fi
    """

# ruleorder: generic_kraken_summary > generic_kraken_translate
    
rule generic_kraken_summary:
    input: kraken = "{base}.kraken"
    output: phylo = "{base}.kraken_phylo"
    output: 
    threads: 1
    shell: """
    
    if [ $(hostname) != "bionc13" ] ; then 
         echo KRAKEN MUST BE RUN ON BIONC13. SLEEPING 30 SECONDS.
         sleep 30
         stop
    fi

    echo 'making report..'
    ## translate file has info on each read (but not in an easy-to-digest way)
    ## phylo file has a summary for each level in the taxonomy
    python3 /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/kraken_report.py --db /mnt/ramdisk/refseqReleaseKraken {input.kraken} > {output.phylo}

    """

rule generic_kraken_translate:
    input: kraken = "{base}.kraken"
    # output: translate = temp("{base}.translate")
    output: translate = "{base}.translate"
    threads: 1
    shell: """
    
    if [ $(hostname) != "bionc13" ] ; then 
         echo KRAKEN MUST BE RUN ON BIONC13. SLEEPING 30 SECONDS.
         sleep 30
         stop
    fi

    echo 'making report..'
    ## translate file has info on each read (but not in an easy-to-digest way)
    ## phylo file has a summary for each level in the taxonomy
    python3 /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/kraken_report.py --db /mnt/ramdisk/refseqReleaseKraken {input.kraken} --translate {output.translate} > /dev/null

    """

## is not set up so you can extract from e.g. swp_masked kraken analyses
rule generic_kraken_extract:
    input: bam = "{base}/{base_name}.bam",
           kraken = "{base}/kraken/{base_name}.bam{kcat}.kraken",
           translate = "{base}/kraken/{base_name}.bam{kcat}.translate"
    output: bam = "{base}/split_kraken/{base_name}.k_{kraken_group}{kcat}.bam"
    wildcard_constraints:
        kcat = "|\.swp_third|\.swp_original|\.swp_masked",
        kraken_group = "[^\.]+"
    threads: 1
    shell: """
    
    if [ $(hostname) != "bionc13" ] ; then 
         echo KRAKEN MUST BE RUN ON BIONC13. SLEEPING 30 SECONDS.
         sleep 30
         stop
    fi

    echo 'finding clade'
    clade=$(awk '$2 == "'{wildcards.kraken_group}'" {{print $1}}' /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/data/clade_taxa_map_alphanum.txt)
    echo clade $clade
    ofile=$(dirname {output.bam})/$(basename {output.bam} .bam)
    echo ofile $ofile
    
    ## save the header first (doesn't kraken_report just overwrite this?)
    # samtools view -H -b {input.bam} > {output.bam}

    echo 'splitting bam..'
    ## translate file has info on each read (but not in an easy-to-digest way)
    ## have to write to the correct file - this doesn't take an ofile, which is maddening...
    python3 /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/kraken_report.py --db /mnt/ramdisk/refseqReleaseKraken {input.kraken} --extractFile {input.bam} \
       --clades $clade --extract-out-base $ofile  > /dev/null
    # --clades $clade --extract-outdir {wildcards.base}/split_kraken/ --suffix {wildcards.kcat} > /dev/null

    """


    
rule generic_kraken_byread:
    input: translate = "{base}.translate"
    output: byread = "{base}.byread"
    threads: 1
    shell: """
    
     echo 'parsing translate file..'
     python3 /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/parse_kraken_translate_file.py --translate-file {input.translate} > {output.byread}

    """

