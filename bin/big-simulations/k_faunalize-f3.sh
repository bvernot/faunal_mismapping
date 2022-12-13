# #!/bin/bash
EIGENSTRAT_subsetDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/12_EIGENSTRAT-subset"
F3Dir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/13_F3Stats"

### Create a new .ind file replacing Population IDs for chosen popuations (which is needed for pairwise analysis)
cat $EIGENSTRAT_subsetDir/poseidon-with-faunal-inds-triallelic-removed-1240k-subset.ind  | awk '{
    if ($3=="Sardinian" || $3=="Han" || $3=="Papuan" || $3=="French" || $3=="Orcadian" || $3=="Italian_North")
    print $1,$2,$1
    else 
    print $1,$2,$3
    }' > $EIGENSTRAT_subsetDir/poseidon-with-faunal-inds-triallelic-removed-1240k-subset-ADJUSTED-POP-IDS.ind

### Make a new directory for F3 Statistics outputs
mkdir -p $F3Dir

### Define common parameters
ind_file=$EIGENSTRAT_subsetDir"/poseidon-with-faunal-inds-triallelic-removed-1240k-subset-ADJUSTED-POP-IDS.ind"
snp_file=$EIGENSTRAT_subsetDir"/poseidon-with-faunal-inds-triallelic-removed-1240k-subset.snp"
original_geno=$EIGENSTRAT_subsetDir"/poseidon-with-faunal-inds-triallelic-removed-1240k-subset.eigenstratgeno"

#############################################################################
##################  Analysis #1 - run on broad populations ##################
#############################################################################
analysis_name="1_SNP-subset-broad_analysis"
mkdir -p $F3Dir/$analysis_name
broad_analysis_all_inds=$(echo HGDP00666 HGDP00667 HGDP00668 HGDP00669 HGDP00670 HGDP00671 HGDP00672 HGDP00774 HGDP00775 HGDP00776 HGDP00777 HGDP00779 HGDP00780 HGDP00781 HGDP00540 HGDP00541 HGDP00543 HGDP00545 HGDP00546 HGDP00547 HGDP00548)

# ### Make query file
# broad_inds_queryfile=$F3Dir/$analysis_name/$analysis_name"_query_file.tsv"
# echo HGDP00666 bear-1 | tr ' ' '\t' > $broad_inds_queryfile
# echo HGDP00667 dog-1 | tr ' ' '\t' >> $broad_inds_queryfile
# echo HGDP00668 sheep-1 | tr ' ' '\t' >> $broad_inds_queryfile
# echo HGDP00774 bear-2 | tr ' ' '\t' >> $broad_inds_queryfile
# echo HGDP00775 dog-2 | tr ' ' '\t' >> $broad_inds_queryfile
# echo HGDP00776 sheep-2 | tr ' ' '\t' >> $broad_inds_queryfile
# echo HGDP00540 bear-3 | tr ' ' '\t' >> $broad_inds_queryfile
# echo HGDP00541 dog-3 | tr ' ' '\t' >> $broad_inds_queryfile
# echo HGDP00543 sheep-3 | tr ' ' '\t' >> $broad_inds_queryfile

# ### Run the faunalizer to create the new .geno files
# ## Not lenient version
# not_lenient_geno=$F3Dir/$analysis_name/"poseidon-with-faunal-inds-triallelic-removed-1240k-subset-FAUNALIZED_not_lenient.eigenstratgeno"
# # python \
# #     /mnt/expressions/robin_warner/3_faunal-mismapping/bin/big-simulations/0_faunalizer/multi-faunalizer_v4.py \
# #     $original_geno \
# #     $ind_file \
# #     $broad_inds_queryfile \
# #     $not_lenient_geno

# # ## Lenient version
# lenient_geno=$F3Dir/$analysis_name/"poseidon-with-faunal-inds-triallelic-removed-1240k-subset-FAUNALIZED_lenient.eigenstratgeno"
# # python \
# #     /mnt/expressions/robin_warner/3_faunal-mismapping/bin/big-simulations/0_faunalizer/multi-faunalizer_v4.py \
# #     $original_geno \
# #     $ind_file \
# #     $broad_inds_queryfile \
# #     $lenient_geno \
# #     -l

# array=( "control" "not_lenient" "lenient" )
# array2=( $original_geno $not_lenient_geno $lenient_geno)


# # ## Make the pop list files
# for sample in $broad_analysis_all_inds; do
#     rm $F3Dir/$analysis_name/$sample.X.Mbuti.F3.list
#     for comparison in $broad_analysis_all_inds; do
#         echo $sample $comparison Mbuti >> $F3Dir/$analysis_name/$sample.X.Mbuti.F3.list
#     done

#     for type in "${!array[@]}"; do
#         val1=${array[type]}
#         val2=${array2[type]}
  
#         # Make directories for separate analysis
#         mkdir -p $F3Dir/$analysis_name/$val1
#         mkdir -p $F3Dir/$analysis_name/$val1/parameter_files
#         mkdir -p $F3Dir/$analysis_name/$val1/log_files

#         # Make parameter files
#         echo -e genotypename: ' \t '    $val2 > $F3Dir/$analysis_name/$val1/parameter_files/par.$sample.X.Mbuti.$val1.F3     
#         echo -e snpname: ' \t ' $snp_file >> $F3Dir/$analysis_name/$val1/parameter_files/par.$sample.X.Mbuti.$val1.F3
#         echo -e indivname: ' \t '  $ind_file >> $F3Dir/$analysis_name/$val1/parameter_files/par.$sample.X.Mbuti.$val1.F3
#         echo -e popfilename: ' \t '    $F3Dir/$analysis_name/$sample.X.Mbuti.F3.list >> $F3Dir/$analysis_name/$val1/parameter_files/par.$sample.X.Mbuti.$val1.F3

#         # Run F3 stats
#         qp3Pop -p $F3Dir/$analysis_name/$val1/parameter_files/par.$sample.X.Mbuti.$val1.F3 > $F3Dir/$analysis_name/$val1/log_files/$sample.X.Mbuti.$val1.F3.log & 
#     done
#     wait
# done

##  Make a heatmap table from F3 outputs
for type in control not_lenient lenient; do    
    ## Tidy up / delete outputs from previous runs
    rm $F3Dir/$analysis_name/$type/log_files/*inter.vals
    rm $F3Dir/$analysis_name/$type/log_files/*inter.csv
    rm $F3Dir/$analysis_name/$type/*temp*

    for sample in $broad_analysis_all_inds; do

    ## Adjust the log files
        cat $F3Dir/$analysis_name/$type/log_files/$sample.X.Mbuti.$type.F3.log | grep result | awk '{print$3,$5}' > $F3Dir/$analysis_name/$type/log_files/$sample.X.Mbuti.$type.F3_ADJUSTED_LOG.csv
        cat $F3Dir/$analysis_name/$type/log_files/$sample.X.Mbuti.$type.F3.log | grep "no data" | awk '{print$4,$6}' >> $F3Dir/$analysis_name/$type/log_files/$sample.X.Mbuti.$type.F3_ADJUSTED_LOG.csv
        
        ## Re-format the adjusted log .csv into an appropriate table for plotting
        for all_inds in $broad_analysis_all_inds; do      
            cat $F3Dir/$analysis_name/$type/log_files/$sample.X.Mbuti.$type.F3_ADJUSTED_LOG.csv | grep -w $all_inds >> $F3Dir/$analysis_name/$type/log_files/${sample}.inter.csv
        done
        
        ## Move the data around to make a table in the correct heatmap format
        sort -k 1 $F3Dir/$analysis_name/$type/log_files/${sample}.inter.csv > $F3Dir/$analysis_name/$type/log_files/temp_${sample}.inter.csv
        mv $F3Dir/$analysis_name/$type/log_files/temp_${sample}.inter.csv $F3Dir/$analysis_name/$type/log_files/${sample}.inter.csv
        awk '{print$2}' $F3Dir/$analysis_name/$type/log_files/${sample}.inter.csv > $F3Dir/$analysis_name/$type/log_files/${sample}.inter.vals
        
    done    
    
    for all_inds in $broad_analysis_all_inds; do
        echo $all_inds
    done | sort > $F3Dir/$analysis_name/$type/temp.list      
        
    cat '\t' $F3Dir/$analysis_name/$type/temp.list | tr '\n' '\t' > $F3Dir/$analysis_name/$type/temp.2.list
    echo "Population:" > $F3Dir/$analysis_name/$type/Population.word
    paste $F3Dir/$analysis_name/$type/Population.word $F3Dir/$analysis_name/$type/temp.2.list > $F3Dir/$analysis_name/$type/temp.3.list
    
    for line in $(cat $F3Dir/$analysis_name/$type/temp.2.list); do
        cat $F3Dir/$analysis_name/$type/log_files/$line.inter.vals | tr '\n' '\t'; echo "\n" | cut -f1 -d '\' -
    done > $F3Dir/$analysis_name/$type/$analysis_name-$type-Full_Table_of_Values.csv
    
    paste $F3Dir/$analysis_name/$type/temp.3.list > $F3Dir/$analysis_name/$type/$analysis_name-$type-HEATMAP_TABLE.csv
    paste $F3Dir/$analysis_name/$type/temp.list $F3Dir/$analysis_name/$type/$analysis_name-$type-Full_Table_of_Values.csv >> $F3Dir/$analysis_name/$type/$analysis_name-$type-HEATMAP_TABLE.csv
        
done


# # ##############################################################################
# # ####################  Analysis #2 - Modern Europeans Only ####################
# # ##############################################################################
analysis_name="2_SNP-subset-modern-europe"
mkdir -p $F3Dir/$analysis_name
modern_europeans_all_inds=$(echo HGDP00666 HGDP00667 HGDP00668 HGDP00669 HGDP00670 HGDP00671 HGDP00672 HGDP00511 HGDP00512 HGDP00513 HGDP00514 HGDP00515 HGDP00516 HGDP00517 HGDP00794 HGDP00796 HGDP00797 HGDP00798 HGDP00799 HGDP00800 HGDP00802 HGDP01147 HGDP01151 HGDP01152 HGDP01153 HGDP01155 HGDP01156 HGDP01157)

# # modern_europeans_queryfile=$F3Dir/$analysis_name/$analysis_name"_query_file.tsv"
# # echo HGDP00666 bear-1 | tr ' ' '\t' > $modern_europeans_queryfile
# # echo HGDP00667 dog-1 | tr ' ' '\t' >> $modern_europeans_queryfile
# # echo HGDP00668 sheep-1 | tr ' ' '\t' >> $modern_europeans_queryfile
# # echo HGDP00511 bear-2 | tr ' ' '\t' >> $modern_europeans_queryfile
# # echo HGDP00512 dog-2 | tr ' ' '\t' >> $modern_europeans_queryfile
# # echo HGDP00513 sheep-2 | tr ' ' '\t' >> $modern_europeans_queryfile
# # echo HGDP00794 bear-3 | tr ' ' '\t' >> $modern_europeans_queryfile
# # echo HGDP00796 dog-3 | tr ' ' '\t' >> $modern_europeans_queryfile
# # echo HGDP00797 sheep-3 | tr ' ' '\t' >> $modern_europeans_queryfile

# # ### Run the faunalizer to create the new .geno files
# # ## Not lenient version
# # not_lenient_geno=$F3Dir/$analysis_name/"poseidon-with-faunal-inds-triallelic-removed-1240k-subset-FAUNALIZED_not_lenient.eigenstratgeno"
# # python \
# #     /mnt/expressions/robin_warner/3_faunal-mismapping/bin/big-simulations/0_faunalizer/multi-faunalizer_v4.py \
# #     $original_geno \
# #     $ind_file \
# #     $modern_europeans_queryfile \
# #     $not_lenient_geno

# # # ## Lenient version
# # lenient_geno=$F3Dir/$analysis_name/"poseidon-with-faunal-inds-triallelic-removed-1240k-subset-FAUNALIZED_lenient.eigenstratgeno"
# # python \
# #     /mnt/expressions/robin_warner/3_faunal-mismapping/bin/big-simulations/0_faunalizer/multi-faunalizer_v4.py \
# #     $original_geno \
# #     $ind_file \
# #     $modern_europeans_queryfile \
# #     $lenient_geno \
# #      -l

# # array=( "control" "not_lenient" "lenient" )
# # array2=( $original_geno $not_lenient_geno $lenient_geno)


# # # ## Make the pop list files
# # for sample in $modern_europeans_all_inds; do
# #     rm $F3Dir/$analysis_name/$sample.X.Mbuti.F3.list
# #     for comparison in $modern_europeans_all_inds; do
# #         echo $sample $comparison Mbuti >> $F3Dir/$analysis_name/$sample.X.Mbuti.F3.list
# #     done

# #     for type in "${!array[@]}"; do
# #         val1=${array[type]}
# #         val2=${array2[type]}
  
# #         # Make directories for separate analysis
# #         mkdir -p $F3Dir/$analysis_name/$val1
# #         mkdir -p $F3Dir/$analysis_name/$val1/parameter_files
# #         mkdir -p $F3Dir/$analysis_name/$val1/log_files

# #         # Make parameter files
# #         echo -e genotypename: ' \t '    $val2 > $F3Dir/$analysis_name/$val1/parameter_files/par.$sample.X.Mbuti.$val1.F3     
# #         echo -e snpname: ' \t ' $snp_file >> $F3Dir/$analysis_name/$val1/parameter_files/par.$sample.X.Mbuti.$val1.F3
# #         echo -e indivname: ' \t '  $ind_file >> $F3Dir/$analysis_name/$val1/parameter_files/par.$sample.X.Mbuti.$val1.F3
# #         echo -e popfilename: ' \t '    $F3Dir/$analysis_name/$sample.X.Mbuti.F3.list >> $F3Dir/$analysis_name/$val1/parameter_files/par.$sample.X.Mbuti.$val1.F3

# #         # Run F3 stats
# #         qp3Pop -p $F3Dir/$analysis_name/$val1/parameter_files/par.$sample.X.Mbuti.$val1.F3 > $F3Dir/$analysis_name/$val1/log_files/$sample.X.Mbuti.$val1.F3.log & 
# #     done
# #     wait
# # done

##  Make a heatmap table from F3 outputs
for type in control not_lenient lenient; do    
    ## Tidy up / delete outputs from previous runs
    rm $F3Dir/$analysis_name/$type/log_files/*inter.vals
    rm $F3Dir/$analysis_name/$type/log_files/*inter.csv
    rm $F3Dir/$analysis_name/$type/*temp*

    for sample in $modern_europeans_all_inds; do

    ## Adjust the log files
        cat $F3Dir/$analysis_name/$type/log_files/$sample.X.Mbuti.$type.F3.log | grep result | awk '{print$3,$5}' > $F3Dir/$analysis_name/$type/log_files/$sample.X.Mbuti.$type.F3_ADJUSTED_LOG.csv
        cat $F3Dir/$analysis_name/$type/log_files/$sample.X.Mbuti.$type.F3.log | grep "no data" | awk '{print$4,$6}' >> $F3Dir/$analysis_name/$type/log_files/$sample.X.Mbuti.$type.F3_ADJUSTED_LOG.csv
        
        ## Re-format the adjusted log .csv into an appropriate table for plotting
        for all_inds in $modern_europeans_all_inds; do      
            cat $F3Dir/$analysis_name/$type/log_files/$sample.X.Mbuti.$type.F3_ADJUSTED_LOG.csv | grep -w $all_inds >> $F3Dir/$analysis_name/$type/log_files/${sample}.inter.csv
        done
        
        ## Move the data around to make a table in the correct heatmap format
        sort -k 1 $F3Dir/$analysis_name/$type/log_files/${sample}.inter.csv > $F3Dir/$analysis_name/$type/log_files/temp_${sample}.inter.csv
        mv $F3Dir/$analysis_name/$type/log_files/temp_${sample}.inter.csv $F3Dir/$analysis_name/$type/log_files/${sample}.inter.csv
        awk '{print$2}' $F3Dir/$analysis_name/$type/log_files/${sample}.inter.csv > $F3Dir/$analysis_name/$type/log_files/${sample}.inter.vals
        
    done    
    
    for all_inds in $modern_europeans_all_inds; do
        echo $all_inds
    done | sort > $F3Dir/$analysis_name/$type/temp.list      
        
    cat '\t' $F3Dir/$analysis_name/$type/temp.list | tr '\n' '\t' > $F3Dir/$analysis_name/$type/temp.2.list
    echo "Population:" > $F3Dir/$analysis_name/$type/Population.word
    paste $F3Dir/$analysis_name/$type/Population.word $F3Dir/$analysis_name/$type/temp.2.list > $F3Dir/$analysis_name/$type/temp.3.list
    
    for line in $(cat $F3Dir/$analysis_name/$type/temp.2.list); do
        cat $F3Dir/$analysis_name/$type/log_files/$line.inter.vals | tr '\n' '\t'; echo "\n" | cut -f1 -d '\' -
    done > $F3Dir/$analysis_name/$type/$analysis_name-$type-Full_Table_of_Values.csv
    
    paste $F3Dir/$analysis_name/$type/temp.3.list > $F3Dir/$analysis_name/$type/$analysis_name-$type-HEATMAP_TABLE.csv
    paste $F3Dir/$analysis_name/$type/temp.list $F3Dir/$analysis_name/$type/$analysis_name-$type-Full_Table_of_Values.csv >> $F3Dir/$analysis_name/$type/$analysis_name-$type-HEATMAP_TABLE.csv
        
done
