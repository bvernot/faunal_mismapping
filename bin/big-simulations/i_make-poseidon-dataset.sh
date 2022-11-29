#!/bin/bash
plinkfilesDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/7_subsets-vcfs-plink"
poseidondatasetDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/8_poseidon-dataset-plink-format"
mergedPoseidonDir="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/9_merged-poseidon"


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
plink --bfile $poseidondatasetDir/faunalised_inds_from_poseidon_noA1_eq_A2 --merge-list $poseidondatasetDir/faunal-merge.list --make-bed --out $mergedPoseidonDir/poseidon-with-faunal-inds

## Remove triallelic sites from each individual fauna

for speciesDir in $plinkfilesDir"/"*; do
    species=`echo $speciesDir | awk -F "/" '{print $NF}'`
    echo ">" $species
    for i in 1 2 3; do
        plink --bfile $plinkfilesDir/$species/$species-${i}_1240k --remove /mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/9_merged-poseidon/poseidon-with-faunal-inds-merge.missnp --make-bed --out $plinkfilesDir/$species/$species-${i}_1240k.no_triallelic
    done
done
