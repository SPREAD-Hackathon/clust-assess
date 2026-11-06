import pandas as pd
import shutil
import os
import bz2
import re

# Define Paths
path = "find_c0s4_vs_reduced"
type_path = "database"
out_dir = "subset_fastas"
test_path = "testset_community0subcommunity4"
fasta_path = "database/sequences.fasta.bz2"

# Make output directories
os.makedirs(out_dir, exist_ok=True)

# Load and filter Accio results
accio_out = pd.read_csv(f"{path}/summary_report.csv")
pling_results = accio_out[["Sample", "pling_type"]]
clustered = pling_results.dropna()
known_clusters = clustered["pling_type"].loc[clustered["pling_type"] != "novel"].unique().tolist()

# Copy the test fastas
test_samples = clustered["Sample"].unique().tolist()
[shutil.copy2(f"{test_path}/{x}.fasta", f"subset_fastas/{x}.fasta") for x in test_samples]

# Extract the database fastas
database_types = pd.read_csv(f"{type_path}/typing.tsv", sep="\t")
database_samples = database_types["plasmid"].loc[database_types["type"].isin(known_clusters)].tolist()

def safe_name(name):
    """Make a header usable as a filename."""
    return re.sub(r"[^\w.\-]", "_", name)


def header_id(header_line):
    """Header line without '>' -> first whitespace-delimited token.
    Change this if your names in database_samples are the FULL header."""
    return header_line[1:].strip().split()[0]


found = set()
out_handle = None

with bz2.open(fasta_path, "rt") as fh:          # "rt" = text mode, streamed
    for line in fh:
        if line.startswith(">"):
            # close the previous output file
            if out_handle:
                out_handle.close()
                out_handle = None

            name = header_id(line)
            if name in database_samples:
                found.add(name)
                out_handle = open(os.path.join(out_dir, f"{safe_name(name)}.fasta"), "w")
                out_handle.write(line)
        elif out_handle:
            out_handle.write(line)

if out_handle:
    out_handle.close()