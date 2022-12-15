# #!/bin/bash
plinkfilesDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/7_subsets-vcfs-plink"
full_plink_file="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/9_merged-poseidon/poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only"
plink_subset_file="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/11_plink-subset-sites"
EIGENSTRAT_subsetDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/12_EIGENSTRAT-subset"


# ## Make a list of the simulated fauna to extract
# rm $plink_subset_file/simulated-fauna.list
# for speciesDir in $plinkfilesDir"/"*; do
#     species=`echo $speciesDir | awk -F "/" '{print $NF}'`
#     echo $species'-1' $species'-1' >> $plink_subset_file/simulated-fauna.list
#     echo $species'-2' $species'-2' >> $plink_subset_file/simulated-fauna.list
#     echo $species'-3' $species'-3' >> $plink_subset_file/simulated-fauna.list
#     done

# ## Make a dataset with just these simulated fauna
# echo $full_plink_file
# plink --bfile $full_plink_file --keep $plink_subset_file/simulated-fauna.list --make-bed --out $plink_subset_file/poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only-FAUNA-EXTRACTED

# ## Find the sites that are informative in the three fauna
# plink --bfile $plink_subset_file/poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only-FAUNA-EXTRACTED --missing --out $plink_subset_file/poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only-FAUNA-EXTRACTED
# cat $plink_subset_file/poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only-FAUNA-EXTRACTED.lmiss | awk '{
#     if($5!='1')
#     print $1,$2 
#     }'> $plink_subset_file/list-of-faunal-sites.txt

# ## Remove those sites from the full list
# plink --bfile $full_plink_file --exclude $plink_subset_file/list-of-informative-sites.txt --make-bed --out $plink_subset_file/poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only-NO-FAUNAL-SITES

############################# Select 50% of sites #############################
# cat $plink_subset_file/poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only-NO-FAUNAL-SITES.bim | shuf | head -514754 | sort > $plink_subset_file/non-faunal-sites-514754.txt

# ## Add sites together and make a dataset from it
# cat $plink_subset_file/list-of-faunal-sites.txt > $plink_subset_file/subset-database-sites.txt
# cat $plink_subset_file/non-faunal-sites-514754.txt >> $plink_subset_file/subset-database-sites.txt

# plink --bfile $full_plink_file --extract $plink_subset_file/subset-database-sites.txt --make-bed --out $plink_subset_file/poseidon-with-faunal-inds-triallelic-removed-1240k-subset
# plink --bfile $plink_subset_file/poseidon-with-faunal-inds-triallelic-removed-1240k-subset --recode --out $plink_subset_file/poseidon-with-faunal-inds-triallelic-removed-1240k-subset

# # Convert the subset database to EIGENSTRAT
# mkdir -p $EIGENSTRAT_subsetDir
# for i in poseidon-with-faunal-inds-triallelic-removed-1240k-subset; do 
#     echo -e genotypename: ' \t ' $plink_subset_file/$i.ped;
#     echo -e snpname: ' \t ' $plink_subset_file/$i.map ;
#     echo -e indivname: ' \t ' $plink_subset_file/$i.pedind;
#     echo -e outputformat: ' \t ' EIGENSTRAT;
#     echo -e genotypeoutname: ' \t '  $EIGENSTRAT_subsetDir/$i.eigenstratgeno;
#     echo -e snpoutname: ' \t '  $EIGENSTRAT_subsetDir/$i.snp;
#     echo -e indivoutname: ' \t '  $EIGENSTRAT_subsetDir/$i.ind;
#     echo -e familynames: ' \t '  NO;
# done > $EIGENSTRAT_subsetDir/par.poseidon-with-faunal-inds-triallelic-removed-1240k-subset.convertf

# for i in poseidon-with-faunal-inds-triallelic-removed-1240k-subset; do
#     more $plink_subset_file/${i}.fam | awk '{print$1,$2,$3,$4,$5, 1}' > $plink_subset_file/${i}.pedind
#     convertf -p $EIGENSTRAT_subsetDir/par.${i}.convertf> $EIGENSTRAT_subsetDir/${i}.convertf.log ;
#     paste $EIGENSTRAT_subsetDir/${i}.ind $plink_subset_file/${i}.pedind | awk '{print $1,$2,$4}' >$EIGENSTRAT_subsetDir/temp; 
#     mv $EIGENSTRAT_subsetDir/temp $EIGENSTRAT_subsetDir/${i}.ind;
#     awk '{print$3}' $EIGENSTRAT_subsetDir/${i}.ind > $EIGENSTRAT_subsetDir/${i}.poplist
# done

############################# Select 60% of sites #############################
# cat $plink_subset_file/poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only-NO-FAUNAL-SITES.bim | shuf | head -636546 | sort > $plink_subset_file/non-faunal-sites-636546.txt

# ## Add sites together and make a dataset from it
# cat $plink_subset_file/list-of-faunal-sites.txt > $plink_subset_file/subset-database-60-per-cent-sites.txt
# cat $plink_subset_file/non-faunal-sites-636546.txt >> $plink_subset_file/subset-database-60-per-cent-sites.txt

# plink --bfile $full_plink_file --extract $plink_subset_file/subset-database-60-per-cent-sites.txt --make-bed --out $plink_subset_file/poseidon-with-faunal-inds-triallelic-removed-1240k-subset-60-per-cent
# plink --bfile $plink_subset_file/poseidon-with-faunal-inds-triallelic-removed-1240k-subset-60-per-cent --recode --out $plink_subset_file/poseidon-with-faunal-inds-triallelic-removed-1240k-subset-60-per-cent

# # Convert the subset database to EIGENSTRAT
mkdir -p $EIGENSTRAT_subsetDir
for i in poseidon-with-faunal-inds-triallelic-removed-1240k-subset-60-per-cent; do 
    echo -e genotypename: ' \t ' $plink_subset_file/$i.ped;
    echo -e snpname: ' \t ' $plink_subset_file/$i.map ;
    echo -e indivname: ' \t ' $plink_subset_file/$i.pedind;
    echo -e outputformat: ' \t ' EIGENSTRAT;
    echo -e genotypeoutname: ' \t '  $EIGENSTRAT_subsetDir/$i.eigenstratgeno;
    echo -e snpoutname: ' \t '  $EIGENSTRAT_subsetDir/$i.snp;
    echo -e indivoutname: ' \t '  $EIGENSTRAT_subsetDir/$i.ind;
    echo -e familynames: ' \t '  NO;
done > $EIGENSTRAT_subsetDir/par.poseidon-with-faunal-inds-triallelic-removed-1240k-subset-60-per-cent.convertf

for i in poseidon-with-faunal-inds-triallelic-removed-1240k-subset-60-per-cent; do
    more $plink_subset_file/${i}.fam | awk '{print$1,$2,$3,$4,$5, 1}' > $plink_subset_file/${i}.pedind
    convertf -p $EIGENSTRAT_subsetDir/par.${i}.convertf> $EIGENSTRAT_subsetDir/${i}.convertf.log ;
    paste $EIGENSTRAT_subsetDir/${i}.ind $plink_subset_file/${i}.pedind | awk '{print $1,$2,$4}' >$EIGENSTRAT_subsetDir/temp; 
    mv $EIGENSTRAT_subsetDir/temp $EIGENSTRAT_subsetDir/${i}.ind;
    awk '{print$3}' $EIGENSTRAT_subsetDir/${i}.ind > $EIGENSTRAT_subsetDir/${i}.poplist
done