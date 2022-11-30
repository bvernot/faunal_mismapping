#!/bin/bash
plinkfilesDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/7_subsets-vcfs-plink"
poseidondatasetDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/8_poseidon-dataset-plink-format"
mergedPoseidonDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/9_merged-poseidon"
EIGENSTRATDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/10_EIGENSTRAT-file"

## Making the Poseidon dataset of modern individuals with references (NOTE - trident needs to be run in a conda environment)
mkdir -p $poseidondatasetDir
#trident forge -d /home/niall_cooke/poseidon_repo -o $poseidondatasetDir -n faunalised_inds_from_poseidon -f '*2012_PattersonGenetics*,*Reference_Genomes*,*Archaic_Humans*' --outFormat PLINK

## Remove sites from Poseidon in which A1 is the same as A2
# cat $poseidondatasetDir/faunalised_inds_from_poseidon.bim | awk '{
#     if ($5==$6)
#     print $0
#     }' > $poseidondatasetDir/sites-with-A1-eq-A2.list

# plink --bfile $poseidondatasetDir/faunalised_inds_from_poseidon --exclude $poseidondatasetDir/sites-with-A1-eq-A2.list --make-bed --out $poseidondatasetDir/faunalised_inds_from_poseidon_noA1_eq_A2

## Reformat site position IDs in .bim file for Poseidon
# cat $poseidondatasetDir/faunalised_inds_from_poseidon_noA1_eq_A2.bim | awk '{print$1,$1":"$4,$3,$4,$5,$6}' | tr ' ' '\t'> $poseidondatasetDir/temp
# mv $poseidondatasetDir/faunalised_inds_from_poseidon_noA1_eq_A2.bim $poseidondatasetDir/faunalised_inds_from_poseidon_noA1_eq_A2.bim_org
# mv $poseidondatasetDir/temp $poseidondatasetDir/faunalised_inds_from_poseidon_noA1_eq_A2.bim

## Reformat the site positions in each of the faunal .bim files
# for speciesDir in $plinkfilesDir"/"*; do
#     species=`echo $speciesDir | awk -F "/" '{print $NF}'`
#     echo ">" $species
#     for i in 1 2 3; do
#         cat $plinkfilesDir/$species/$species-${i}_1240k.bim | awk '{print$1,$1":"$4,$3,$4,$5,$6}' | tr ' ' '\t'> $plinkfilesDir/$species/$species-${i}_TEMP
#         mv $plinkfilesDir/$species/$species-${i}_1240k.bim $plinkfilesDir/$species/$species-${i}_1240k.bim_org
#         mv $plinkfilesDir/$species/$species-${i}_TEMP $plinkfilesDir/$species/$species-${i}_1240k.bim
#     done
# done

## Make the .list file for merging
# rm $poseidondatasetDir/faunal-merge.list
# for speciesDir in $plinkfilesDir"/"*; do
#     species=`echo $speciesDir | awk -F "/" '{print $NF}'`
#     #echo ">" $species
#     echo $plinkfilesDir/$species/$species-1_1240k.bed $plinkfilesDir/$species/$species-1_1240k.bim $plinkfilesDir/$species/$species-1_1240k.fam >> $poseidondatasetDir/faunal-merge.list
#     echo $plinkfilesDir/$species/$species-2_1240k.bed $plinkfilesDir/$species/$species-2_1240k.bim $plinkfilesDir/$species/$species-2_1240k.fam >> $poseidondatasetDir/faunal-merge.list
#     echo $plinkfilesDir/$species/$species-3_1240k.bed $plinkfilesDir/$species/$species-3_1240k.bim $plinkfilesDir/$species/$species-3_1240k.fam >> $poseidondatasetDir/faunal-merge.list
# done

## Merging the Poseidon dataset to the nine faunal inds
mkdir -p $mergedPoseidonDir
#plink --bfile $poseidondatasetDir/faunalised_inds_from_poseidon_noA1_eq_A2 --merge-list $poseidondatasetDir/faunal-merge.list --make-bed --out $mergedPoseidonDir/poseidon-with-faunal-inds

## Remove triallelic sites from each individual fauna; check missingness before and after

#  for speciesDir in $plinkfilesDir"/"*; do
#      species=`echo $speciesDir | awk -F "/" '{print $NF}'`
#      echo ">" $species
#      for i in 1 2 3; do
#          plink --bfile $plinkfilesDir/$species/$species-${i}_1240k --missing --out $plinkfilesDir/$species/$species-${i}_1240k
#          #plink --bfile $plinkfilesDir/$species/$species-${i}_1240k --exclude /mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/9_merged-poseidon/poseidon-with-faunal-inds-merge.missnp --make-bed --out $plinkfilesDir/$species/$species-${i}_1240k.no_triallelic
#          plink --bfile $plinkfilesDir/$species/$species-${i}_1240k.no_triallelic --missing --out $plinkfilesDir/$species/$species-${i}_1240k.no_triallelic
#      done
#  done

## Make the .list file for merging
#  rm $poseidondatasetDir/faunal-merge-post-triallelic.list
#  for speciesDir in $plinkfilesDir"/"*; do
#      species=`echo $speciesDir | awk -F "/" '{print $NF}'`
#      #echo ">" $species
#      echo $plinkfilesDir/$species/$species-1_1240k.no_triallelic.bed $plinkfilesDir/$species/$species-1_1240k.no_triallelic.bim $plinkfilesDir/$species/$species-1_1240k.no_triallelic.fam >> $poseidondatasetDir/faunal-merge-post-triallelic.list
#      echo $plinkfilesDir/$species/$species-2_1240k.no_triallelic.bed $plinkfilesDir/$species/$species-2_1240k.no_triallelic.bim $plinkfilesDir/$species/$species-2_1240k.no_triallelic.fam >> $poseidondatasetDir/faunal-merge-post-triallelic.list
#      echo $plinkfilesDir/$species/$species-3_1240k.no_triallelic.bed $plinkfilesDir/$species/$species-3_1240k.no_triallelic.bim $plinkfilesDir/$species/$species-3_1240k.no_triallelic.fam >> $poseidondatasetDir/faunal-merge-post-triallelic.list
#  done

#  ## Merging the Poseidon dataset to the nine faunal inds with triallelic sites removed
#  plink --bfile $poseidondatasetDir/faunalised_inds_from_poseidon_noA1_eq_A2 --merge-list $poseidondatasetDir/faunal-merge-post-triallelic.list --make-bed --out $mergedPoseidonDir/poseidon-with-faunal-inds-triallelic-removed

 ## Filter for 1240k sites only and check missingness
#plink --bfile $mergedPoseidonDir/poseidon-with-faunal-inds-triallelic-removed --extract $poseidondatasetDir/faunalised_inds_from_poseidon_noA1_eq_A2.bim --make-bed --out $mergedPoseidonDir/poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only
#plink --bfile $mergedPoseidonDir/poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only --missing --out $mergedPoseidonDir/poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only

## Re-label the .fam file for the faunal inds
# cat $mergedPoseidonDir/poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only.fam | awk '{
#     if ($2=="1240k")
#     print $1,$1,$3,$4,$5,$6
#     else
#     print $1,$2,$3,$4,$5,$6
#     }' > $mergedPoseidonDir/temp

# mv $mergedPoseidonDir/poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only.fam $mergedPoseidonDir/poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only.fam_org
# mv $mergedPoseidonDir/temp $mergedPoseidonDir/poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only.fam

# plink --bfile $mergedPoseidonDir/poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only --recode --out $mergedPoseidonDir/poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only
# plink --file $mergedPoseidonDir/poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only --recode --out $mergedPoseidonDir/poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only

# ## Convert to EIGENSTRAT format - make the parameter file and re-format the data
mkdir -p $EIGENSTRATDir
for i in poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only; do 
    echo -e genotypename: ' \t ' $mergedPoseidonDir/$i.ped;
    echo -e snpname: ' \t ' $mergedPoseidonDir/$i.map ;
    echo -e indivname: ' \t ' $mergedPoseidonDir/$i.pedind;
    echo -e outputformat: ' \t ' EIGENSTRAT;
    echo -e genotypeoutname: ' \t '  $EIGENSTRATDir/$i.eigenstratgeno;
    echo -e snpoutname: ' \t '  $EIGENSTRATDir/$i.snp;
    echo -e indivoutname: ' \t '  $EIGENSTRATDir/$i.ind;
    echo -e familynames: ' \t '  NO;
done > $EIGENSTRATDir/par.poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only.convertf

 for i in poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only; do
    more $mergedPoseidonDir/${i}.fam | awk '{print$1,$2,$3,$4,$5, 1}' > $mergedPoseidonDir/${i}.pedind
    convertf -p $EIGENSTRATDir/par.${i}.convertf> $EIGENSTRATDir/${i}.convertf.log ;
    paste $EIGENSTRATDir/${i}.ind $mergedPoseidonDir/${i}.pedind | awk '{print $1,$2,$4}' >$EIGENSTRATDir/temp; 
    mv $EIGENSTRATDir/temp $EIGENSTRATDir/${i}.ind;
    awk '{print$3}' $EIGENSTRATDir/${i}.ind > $EIGENSTRATDir/${i}.poplist
done