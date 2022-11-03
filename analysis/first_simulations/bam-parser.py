import pysam
import pandas as pd

filename = "/mnt/expressions/benjamin_vernot/faunal_mismapping/data/simulated-reads/sorted_test.bam"

data = {
    "score": [],
    "len": []
}

bam = pysam.AlignmentFile(filename, "rb")
for index, read in enumerate(bam.fetch()):
    if index >= 10:
        break
    else:
        data["score"].append(read.mapping_quality)
        data["len"].append(read.infer_read_length())
        print(read.get_forward_qualities())


df = pd.DataFrame.from_dict(data)
print(df)



bam.close()