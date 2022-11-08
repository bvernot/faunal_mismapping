import pysam
import pandas as pd
import sys

_, bamfile, tsvfile = sys.argv

data = {
    "score": [],
    "len": [],
    "mismatch": []
}

bam = pysam.AlignmentFile(bamfile, "rb")
for index, read in enumerate(bam.fetch()):
    data["score"].append(read.mapping_quality)
    data["len"].append(read.infer_read_length())
    data["mismatch"].append(read.get_tag("NM"))

# bam.close()

df = pd.DataFrame.from_dict(data)
df.to_csv(tsvfile, sep="\t", index=False)
