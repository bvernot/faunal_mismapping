import pandas as pd
import numpy as np

eigenstrat_file = "/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/10_EIGENSTRAT-file/poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only.eigenstratgeno"
individual_file = "/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/10_EIGENSTRAT-file/poseidon-with-faunal-inds-triallelic-removed-1240k-sites-only.ind"
query_file = "/mnt/expressions/robin_warner/3_faunal-mismapping/bin/big-simulations/faunalization_query.tsv"

individual_df = pd.read_csv(individual_file, delimiter=" ", names=["name", "sex", "group"])
eigenstrat_df = pd.read_csv(eigenstrat_file, nrows=5, names=["full_str"])

# loading and printing query
query_df = pd.read_csv(query_file, delimiter="\t", names=["target", "source"])
target_df = pd.merge(left=query_df[["target"]], left_on="target", right=individual_df, right_on="name", how="left")
source_df = pd.merge(left=query_df[["source"]], left_on="source", right=individual_df, right_on="name", how="left")

print(f"Faunalizing {query_df.shape[0]} individuals")
for (_, target_row), (_, source_row) in zip(target_df.iterrows(), source_df.iterrows()):
    print(f"> {target_row['name']} [{target_row['sex']}, {target_row['group']}]")

# print("reading")
# eigenstrat_df = 

# print("converting")
# new_df = eigenstrat_df[0].str.split("", expand=True)

# print(new_df)

# arr = np.loadtxt(eigenstrat_file, dtype=str, delimiter="", max_rows=3)

# print(arr.shape)

# print(arr)