import pandas as pd
import argparse


## SET UP PARSER ##
parser = argparse.ArgumentParser(description='Fast "faunalizer" for multiple inputs.')
parser.add_argument("eigenstrat", help="the input .eigenstratget file")
parser.add_argument("individuals", help="the input .ind file containing all individuals from the eigenstrat")
parser.add_argument("query", help="the input .tsv file containing target and source individual names. Target will be changed by source. Example: target1'\\t'source1'\\n'target2'\\t'source2'")
parser.add_argument("output", help="the output .eigenstratget file")
parser.add_argument("-l", "--lenient", required=False, action="store_true", help="target SNP-positions with no information ('9' in eigenstrat) will still be faunalized")
parser.add_argument("-q", "--quiet", required=False, default=False, help="suppresses output")
args = parser.parse_args()


## LOAD INDIVIDUAL DATA ##
individual_df = pd.read_csv(args.individuals, delimiter=" ", names=["name", "sex", "group"])
individual_df["index"] = individual_df.index


## LOAD QUERY INFORMATION ##
query_df = pd.read_csv(args.query, delimiter="\t", names=["target", "source"])
target_df = pd.merge(left=query_df[["target"]], left_on="target", right=individual_df, right_on="name", how="left")
source_df = pd.merge(left=query_df[["source"]], left_on="source", right=individual_df, right_on="name", how="left")


## CHECK IF INDIVIDUALS EXIST ##
all_individuals = pd.concat([query_df["target"], query_df["source"]], axis=0)
for individual, is_valid in zip(all_individuals, all_individuals.isin(individual_df["name"])):
    if not is_valid:
        raise ValueError(f"Individual '{individual}' does not exist!")


## PRINT QUERY INFORMATION ##
if not args.quiet:
    print(f"Faunalizing {query_df.shape[0]} individuals")
    for (_, target_row), (_, source_row) in zip(target_df.iterrows(), source_df.iterrows()):
        print(f"> {target_row['name']} [{target_row['index']}, {target_row['sex']}, {target_row['group']}], faunalized by {source_row['name']} [{source_row['index']}, {source_row['sex']}, {source_row['group']}]")


## READ EIGENSTRAT ##
eigenstrat_df = pd.read_csv(args.eigenstrat, names=["full_str"])


## EXTRACT SNPS WHERE WE HAVE DATA ##
for ind, ((_, target_row), (_, source_row)) in enumerate(zip(target_df.iterrows(), source_df.iterrows())):
    eigenstrat_df[ind] = (eigenstrat_df["full_str"].str.slice(target_row["index"], target_row["index"]+1) != eigenstrat_df["full_str"].str.slice(source_row["index"], source_row["index"]+1)) & \
                         (eigenstrat_df["full_str"].str.slice(target_row["index"], target_row["index"]+1) != "9")

    if not args.lenient:
        eigenstrat_df[ind] = eigenstrat_df[ind] & (eigenstrat_df["full_str"].str.slice(source_row["index"], source_row["index"]+1) != "9")


## SPLIT THE STRINGS INTO COLUMNS FOR EACH INDIVIDUAL ##
reduced_eigenstrat_df = eigenstrat_df[eigenstrat_df[query_df.index].any(axis=1)]
split_eigenstrat_df = reduced_eigenstrat_df["full_str"].str.split("", expand=True)


## PRINTING SPLITTING INFO #
n_eigenstrat = eigenstrat_df.shape[0]
max_namelen = max(target_df["name"].str.len().max(), source_df["name"].str.len().max(), 7)
n_number = len(f"{n_eigenstrat:_}")

if not args.quiet:
    print(f"Faunalizing {reduced_eigenstrat_df.shape[0]:_} of {n_eigenstrat:_} SNPs ({reduced_eigenstrat_df.shape[0]/n_eigenstrat:.2e}), excluding sites where target==query")
    for ind, ((_, target_row), (_, source_row)) in enumerate(zip(target_df.iterrows(), source_df.iterrows())):
        n_neq9_target = (eigenstrat_df['full_str'].str.slice(target_row['index'], target_row['index']+1) != '9').sum()
        n_neq9_source = (eigenstrat_df['full_str'].str.slice(source_row['index'], source_row['index']+1) != '9').sum()
        n_neq9_both = ((eigenstrat_df['full_str'].str.slice(target_row['index'], target_row['index']+1) != '9') & \
                    (eigenstrat_df['full_str'].str.slice(source_row['index'], source_row['index']+1) != '9')).sum()
        print(f"> n({target_row['name']:<{max_namelen}} != 9) = {n_neq9_target:>{n_number}_} ({n_neq9_target/n_eigenstrat:.2e})  \ttarget")
        print(f"  n({source_row['name']:<{max_namelen}} != 9) = {n_neq9_source:>{n_number}_} ({n_neq9_source/n_eigenstrat:.2e})  \tsource")
        print(f"  overlap{' '*(max_namelen+1)} = {n_neq9_both:>{n_number}_} ({n_neq9_both/n_eigenstrat:.2e})")


## CHANGE THE NUCLEOTIDES ##
for ind, ((_, target_row), (_, source_row)) in enumerate(zip(target_df.iterrows(), source_df.iterrows())):
    split_eigenstrat_df.loc[reduced_eigenstrat_df[ind] == True, target_row["index"]+1] = split_eigenstrat_df[source_row["index"]+1]


## RECOMBINE THE SPLIT COLUMNS AND MERGE BACK ##
split_eigenstrat_df["unsplit"] = split_eigenstrat_df.values.sum(axis=1)
eigenstrat_df.loc[split_eigenstrat_df.index, "full_str"] = split_eigenstrat_df["unsplit"]


## SAVE DATA TO DISK ##
eigenstrat_df["full_str"].to_csv(args.output, header=False, index=False)
