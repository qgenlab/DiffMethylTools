import pandas as pd
import sys
import os

path = sys.argv[1]
df = pd.read_csv(path, sep="\t", names=["chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes"])
output = pd.DataFrame()

output["chr"] = df["chrom"]
output["pos"] = df["chromStart"]
output["N"] = df["blockCount"]
output["X"] = (df["blockSizes"]/100) * df["blockCount"]

output["X"] = output["X"].astype(int)

f_output = path.rpartition('.')[0]

output.to_csv(f_output+"_dss_stranded.txt", sep="\t", index=False)
