import pandas as pd

eigenstrat_file = "/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/10_EIGENSTRAT-file/poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only.eigenstratgeno"
individual_file = "/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/10_EIGENSTRAT-file/poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only.ind"
query_file = "/mnt/expressions/robin_warner/3_faunal-mismapping/bin/big-simulations/faunalization_query.tsv"
output_file = "/mnt/expressions/robin_warner/3_faunal-mismapping/bin/big-simulations/0_faunalizer/result.eigenstratgeno"

## READ EIGENSTRAT ##
eigenstrat_df = pd.read_csv(eigenstrat_file, names=["full_str"])

## LOAD INDIVIDUAL DATA ##
individual_df = pd.read_csv(individual_file, delimiter=" ", names=["name", "sex", "group"])
individual_df["index"] = individual_df.index

## LOAD QUERY INFORMATION ##
query_df = pd.read_csv(query_file, delimiter="\t", names=["target", "source"])
target_df = pd.merge(left=query_df[["target"]], left_on="target", right=individual_df, right_on="name", how="left")
source_df = pd.merge(left=query_df[["source"]], left_on="source", right=individual_df, right_on="name", how="left")

## PRINT QUERY INFORMATION ##
print(f"Faunalizing {query_df.shape[0]} individuals")
for (_, target_row), (_, source_row) in zip(target_df.iterrows(), source_df.iterrows()):
    print(f"> {target_row['name']} [{target_row['index']}, {target_row['sex']}, {target_row['group']}], faunalized by {source_row['name']} [{source_row['index']}, {source_row['sex']}, {source_row['group']}]")

## EXTRACT SNPS WHERE WE HAVE DATA ##
for ind, ((_, target_row), (_, source_row)) in enumerate(zip(target_df.iterrows(), source_df.iterrows())):
    eigenstrat_df[ind] = (eigenstrat_df["full_str"].str.slice(target_row["index"], target_row["index"]+1) != eigenstrat_df["full_str"].str.slice(source_row["index"], source_row["index"]+1)) & \
                         (eigenstrat_df["full_str"].str.slice(target_row["index"], target_row["index"]+1) != "9") & \
                         (eigenstrat_df["full_str"].str.slice(source_row["index"], source_row["index"]+1) != "9")

## SPLIT THE INTO COLUMNS FOR EACH INDIVIDUAL ##
reduced_eigenstrat_df = eigenstrat_df[eigenstrat_df[query_df.index].any(axis=1)]
split_eigenstrat_df = reduced_eigenstrat_df["full_str"].str.split("", expand=True)

## CHANGE THE NUCLEOTIDES ##
for ind, ((_, target_row), (_, source_row)) in enumerate(zip(target_df.iterrows(), source_df.iterrows())):
    split_eigenstrat_df.loc[reduced_eigenstrat_df[ind] == True, target_row["index"]+1] = split_eigenstrat_df[source_row["index"]+1]

## RECOMBINE THE SPLIT COLUMNS AND MERGE BACK ##
split_eigenstrat_df["unsplit"] = split_eigenstrat_df.values.sum(axis=1)
eigenstrat_df.loc[split_eigenstrat_df.index, "full_str"] = split_eigenstrat_df["unsplit"]

## SAVE DATA TO DISK ##
eigenstrat_df["full_str"].to_csv(output_file, header=False, index=False)