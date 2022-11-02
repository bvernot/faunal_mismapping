import argparse
from random import sample

#/mnt/expressions/niall_cooke/SpeedRun/SpeedRun_Dataset_1.geno

parser = argparse.ArgumentParser(prog='ancestralize Eigenstrat', description='ancestralizing sites')
parser._optionals.title = "Required Arguments"
optional_args=parser.add_argument_group('Optional Arguments')
parser.add_argument('--input', "-in", help='input file')
parser.add_argument('--individual', "-ind", type=int, help='input file')
parser.add_argument('--ancestral', "-anc", type=int, help='input file')
parser.add_argument('--prop', "-prop", default = 10, type=int, help='input file')
args=parser.parse_args()

file = open(args.input)

prop_list = [i for i in range(0, args.prop)]
i=0
n=0
for read in file: 
    geno = list(read.replace("\n", ""))
    if geno[args.individual] == "9":       
        print("".join(geno))
        #print(geno)
        continue
    if (geno[args.ancestral] == "1") or (geno[args.ancestral] == "9"):
        print("".join(geno))
        continue
    if geno[args.ancestral] == "0" :
        butterfly_kevin_loci = "0"
        #print(butterfly_kevin_loci)
    elif geno[args.ancestral] == "2" :
        butterfly_kevin_loci = "2" 
    #print(sample(prop_list, 1)[0])
    n=n+1
    if sample(prop_list, 1)[0] == 1:
        geno[args.individual]=geno[args.ancestral]
        print("".join(geno))
        i=i+1 
    else:
        print("".join(geno))

#print(i)
#print(n)