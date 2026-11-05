import pandas as pd
import numpy as np
from pathlib import Path


accio_test_info = "accio/test_runs/db_full430/plasmidDB_info.csv"
accio_test_distances = "accio/test_runs/db_full430/pling/all_plasmids_distances.tsv"
plsdb_version_info = "pling_plsdb/typing.tsv"
plsdb_distances = "pling_plsdb/all_plasmids_distances.tsv"


plasmid_info = pd.read_csv(plsdb_version_info, sep="\t")
plasmid_info["community"] = [f"community_{x.split("_")[1]}" for x in plasmid_info["type"].values]

print(plasmid_info)
#print(community_info)

distance_info = pd.read_csv(plsdb_distances, sep="\t")
#containment_info = pd.read_csv("accio/test_runs/db_full430/pling/containment/all_pairs_containment_distance.tsv", sep="\t")

print(distance_info)
#print(containment_info)

hub_plasmids = pd.read_csv("pling_plsdb/hub_plasmids.csv", sep=';')
print(hub_plasmids)

pairwise_distances = distance_info.merge(plasmid_info, how="left", left_on="plasmid_1", right_on="plasmid", suffixes=("", "_1"))
pairwise_distances = pairwise_distances.drop(columns="plasmid")

pairwise_distances = pairwise_distances.merge(plasmid_info, how="left", left_on="plasmid_2", right_on="plasmid", suffixes=("", "_2"))
pairwise_distances = pairwise_distances.drop(columns="plasmid")
#pairwise_distances = pairwise_distances.merge(containment_info, how="left", left_on=["plasmid_1", "plasmid_2"], right_on=["plasmid_1", "plasmid_2"], suffixes=("", "_contaiment"))


print(pairwise_distances["plasmid_1"].iloc[0], type(pairwise_distances["plasmid_1"].iloc[0]))
print(hub_plasmids["hub_plasmids"].iloc[0], type(hub_plasmids["hub_plasmids"].iloc[0]))

hub_set = set(hub_plasmids["hub_plasmids"].dropna())

pairwise_distances["hub_1"] = pairwise_distances["plasmid_1"].isin(hub_set).astype(str)
pairwise_distances["hub_2"] = pairwise_distances["plasmid_2"].isin(hub_set).astype(str)

print(pairwise_distances)
print(pairwise_distances.columns)

pairwise_clean = pairwise_distances.copy()
pairwise_clean = pairwise_clean.drop(columns=["hub_1", "hub_2", "check_1", "check_2"])
pairwise_clean["type"] = ["hub" if pd.isna(x) else x for x in pairwise_clean["type"].values]
pairwise_clean["community"] = ["hub" if pd.isna(x) else x for x in pairwise_clean["community"].values]
pairwise_clean["type_2"] = ["hub" if pd.isna(x) else x for x in pairwise_clean["type_2"].values]
pairwise_clean["community_2"] = ["hub" if pd.isna(x) else x for x in pairwise_clean["community_2"].values]
print(pairwise_clean)

pairwise_communities = pairwise_clean.loc[pairwise_clean["type"] != "hub"]
pairwise_communities = pairwise_communities.loc[pairwise_clean["type_2"] != "hub"]
print(pairwise_communities)

danish_accessions_path = "clust-assess/denmark_plsdb_accession.tsv"
danish_accessions = pd.read_csv(danish_accessions_path, sep="\t")
print(danish_accessions)

pairwise_danish = pairwise_communities.loc[pairwise_communities["plasmid_1"].isin(danish_accessions["NUCCORE_ACC"])]

df_pling = pairwise_communities[pairwise_communities['type'] == pairwise_communities['type_2']].copy()
df_community = pairwise_communities[pairwise_communities['community'] == pairwise_communities['community_2']].copy()


output = Path("accio")

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
        folder = output / group_col
        folder.mkdir(parents=True, exist_ok=True)

        values = np.sort(
            data["distance"].dropna().to_numpy(dtype=np.float64)
        )

        np.save(folder / f"{group}_dcj.npy", values)

        print(f"Saved {group_col}/{group}_dcj: {len(values):,} observations")