import pandas as pd
import numpy as np
from pathlib import Path

communities = "pling_plsdb/typing.tsv"
all_plasmid_distances = "pling_plsdb/all_plasmids_distances.tsv"
hub_plasmids = "pling_plsdb/hub_plasmids.csv"


# Load data
plasmid_communities = pd.read_csv(communities, sep="\t")
plasmid_communities["community"] = [f"community_{x.split("_")[1]}" for x in plasmid_communities["type"].values]

distance_info = pd.read_csv(all_plasmid_distances, sep="\t")
#containment_info = pd.read_csv("accio/test_runs/db_full430/pling/containment/all_pairs_containment_distance.tsv", sep="\t")

hub_plasmids = pd.read_csv("pling_plsdb/hub_plasmids.csv", sep=';')

# Merge the distances
pairwise_distances = distance_info.merge(plasmid_communities, how="left", left_on="plasmid_1", right_on="plasmid", suffixes=("", "_1"))
pairwise_distances = pairwise_distances.drop(columns="plasmid")
pairwise_distances = pairwise_distances.merge(plasmid_communities, how="left", left_on="plasmid_2", right_on="plasmid", suffixes=("", "_2"))
pairwise_distances = pairwise_distances.drop(columns="plasmid")

# Locate the hub plasmids
hub_set = set(hub_plasmids["hub_plasmids"].dropna())
pairwise_distances["hub_1"] = pairwise_distances["plasmid_1"].isin(hub_set).astype(str)
pairwise_distances["hub_2"] = pairwise_distances["plasmid_2"].isin(hub_set).astype(str)

# Remove the hubs
pairwise_distances["type"] = ["hub" if pd.isna(x) else x for x in pairwise_distances["type"].values]
pairwise_distances["community"] = ["hub" if pd.isna(x) else x for x in pairwise_distances["community"].values]
pairwise_distances["type_2"] = ["hub" if pd.isna(x) else x for x in pairwise_distances["type_2"].values]
pairwise_distances["community_2"] = ["hub" if pd.isna(x) else x for x in pairwise_distances["community_2"].values]
pairwise_communities = pairwise_distances.loc[pairwise_distances["type"] != "hub"]
pairwise_communities = pairwise_communities.loc[pairwise_distances["type_2"] != "hub"]

# Limit to same (sub)communitiy plasmids
df_pling = pairwise_communities[pairwise_communities['type'] == pairwise_communities['type_2']].copy()
df_community = pairwise_communities[pairwise_communities['community'] == pairwise_communities['community_2']].copy()

output = Path("database")

for group_col, group_col_2 in [
    ("community", "community_2"),
    ("type", "type_2")
]:
    if group_col == "":
        df = df_pling.copy()
    else:
        df = df_community.copy()

    # Only retain pairs belonging to the same group
    subset = df.loc[df[group_col] == df[group_col_2]]

    for group, data in subset.groupby(group_col):
        folder = output / f"{group_col}_custom"
        folder.mkdir(parents=True, exist_ok=True)

        values = np.sort(
            data["distance"].dropna().to_numpy(dtype=np.float64)
        )

        np.save(folder / f"{group}_dcj.npy", values)
