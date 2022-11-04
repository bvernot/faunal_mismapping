
## FORMAT: ./ancestral_derived_frequencies_from_plink.sh <plink_database> <population> <ancestral_derived_file>
## e.g. ./ancestral_derived_frequencies_from_plink.sh /mnt/expressions/niall_cooke/SpeedRun/SpeedRun_Dataset_1 Sardinian Ancestral_vs_derived.txt

# How ancestral and derived sites were initially determined using plink --freq and the following: 

#echo "Chromosome Position Ancestral Derived Category" | tr ' ' '\t' > Ancestral_vs_derived.txt
#cat Chimp.frq.strat | awk '{
#if ($8=='0')
#print $1,$2,"MISS","MISS","MISS"
#if ($8=='2' && $6==0)
#print $1,$2,$5,$4,"KEEP"
#if ($8=='2' && $6==1)
#print $1,$2,$4,$5,"SWAP"
#}' | tr ' ' '\t' >> Ancestral_vs_derived.txt



## Select population from the database, make a list file + print the info
cat $1.fam | grep -w $2 | grep -v "ignore" > $2.list
cat $2.list | awk '{print$1,$2}'
cat $2.list | awk '{print$1}' | sort | uniq -c

## Generate frequency stats
plink --bfile $1 --freq --within $2.list -out $2

## Compare the frequencies to the Chimp information file
echo "Chromosome Position Chimp? Ancestral Derived Freq(Derived) MAF" | tr ' ' '\t' > $2_ancestral_derived_freq.txt
paste $3 $2.frq.strat | awk '{ 
if ($5=="MISS")
print$1,$2,"N",$10,$9,"-",$11
if ($5=="KEEP")
print$1,$2,"Y",$10,$9,$11,"-"
if ($5=="SWAP")
print$1,$2,"Y",$9,$10,(1-$11),"-"
}' | tr ' ' '\t' >> $2_ancestral_derived_freq.txt

echo "Output written to " $2_ancestral_derived_freq.txt
