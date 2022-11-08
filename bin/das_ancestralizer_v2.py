###########################################################################################################################################################################################################################################
##	The function of this script is to revert a set proportion of an individual humans's called SNPs to their ancestral state.
##	The input is a ".geno" file (EIGENSTRAT format) containing both the individual of interest and a representitive of the ancestral state (i.e. Chimp). 
##	Individual and ancestral input numbers refer to their position in the order of the corresponding ".ind" file in EIGENSTRAT format.
##	Value for proportion must be between 0 and 1.
##	Note: in .geno format: "0" - 0 copies of the reference allele, "1" - one copy of the reference allele, "2" - two copies of the reference allele, "9" - missing data
##	To run:		python das_ancestrilizer_v2.py  --input <geno file> --individual <individual_position_in_ind> --ancestral <ancestral_position_in_ind> --prop <proportion_of_genome_to_ancestralize> > output.geno
###########################################################################################################################################################################################################################################

import argparse
import random

## Set up the parser
parser = argparse.ArgumentParser(prog='ancestralize Eigenstrat', description='revert a set proportion of an individual humans's called SNPs to their ancestral state')
parser._optionals.title = "Required Arguments"
optional_args=parser.add_argument_group('Optional Arguments')
parser.add_argument('--input', "-in", help='input file')
parser.add_argument('--individual', "-ind", type=int, help='input file')
parser.add_argument('--ancestral', "-anc", type=int, help='input file')
parser.add_argument('--prop', "-prop", default = 0.1, type=float, help='input file')
args=parser.parse_args()

## Adjust the position of individuals in the ind file for counting from 0 in Python
individual_position=int(args.individual-1)
ancestral_position=int(args.ancestral-1)

## Read in the file to 1) get the number of informative sites in an individual and 2) get number of sites in the individual that are "ancestralizable" (i.e. present in the individual while neither missing (9) or heterozygous (1) in the ancestral)
file = open(args.input)
non_missing_sites_in_ind_count = 0
ancestralizable_sites = 0
for read in file:
    geno = list(read.replace("\n", ""))
    if geno[individual_position] != "9":
        non_missing_sites_in_ind_count += 1
    if (geno[individual_position] != "9" and geno[ancestral_position] == "0") or (geno[individual_position] != "9" and geno[ancestral_position] == "2"):
        ancestralizable_sites +=1

## Adjust the specifed proportion to take into account site missing or heterozygous in the ancestral
adjusted_proportion= float(args.prop) / ((float(ancestralizable_sites) / float(non_missing_sites_in_ind_count)))

## Read in the file again to output the un-ancestralizable sites unchanged, and randomly change the adjusted set proportion of ancestralizable sites.
file = open(args.input)
i=0
n=0
for read in file: 
    geno = list(read.replace("\n", ""))
    if geno[individual_position] == "9":
	print("".join(geno))
	continue
    if (geno[ancestral_position] == "1") or (geno[ancestral_position] == "9"):
	print("".join(geno))
	continue
    n=n+1
    if random.random() <= (adjusted_proportion):
        geno[individual_position]=geno[ancestral_position]
        print("".join(geno))
        i=i+1 
    else:
        print("".join(geno))
