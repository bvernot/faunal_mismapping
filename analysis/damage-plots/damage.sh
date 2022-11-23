basedirectory="/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/4_mappings/0_ref"
referencegenome="/mnt/solexa/Genomes/hg19_evan/whole_genome.fa"


for species in 0_human 2_dog; do 
    if $species; then
    echo $species 
    directoryname="$basedirectory/$species"
    inputfilename=$directoryname/*"-to-human_REF.bam"
    
    mapDamage -i $inputfilename -r $referencegenome
done

exit