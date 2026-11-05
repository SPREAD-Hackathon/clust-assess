import pandas as pd
import numpy as np
from pathlib import Path

communities = "pling_plsdb/typing.tsv"
all_plasmid_distances = "outputfolder/pling/pling_dcj_dists_by_length.tsv"
all_plasmid_containment = "outputfolder/pling/containment/all_pairs_containment_distance.tsv"

# Load data
plasmid_communities = pd.read_csv(communities, sep="\t")
plasmid_communities["community"] = [f"community_{x.split("_")[1]}" for x in plasmid_communities["type"].values]

distance_info = pd.read_csv(all_plasmid_distances, sep="\t")
containment_info = pd.read_csv(all_plasmid_distances, sep="\t")


def save_database_files(df_in: pd.DataFrame, group_col1: str, group_col2: str, metric: str):
    df = df_in.copy()

    # Only retain pairs belonging to the same group
    subset = df.loc[df[group_col1] == df[group_col2]]

    for group, data in subset.groupby(group_col1):
        folder = f"database/{group_col}_custom"
        os.makedirs(folder, exist_ok=True)

        values = np.sort(data["distance"].dropna().to_numpy(dtype=np.float64))

        np.save(f"{folder}/{group}_{metric}.npy", values)


def generete_database(distance_file: pd.DataFrame, plasmid_communities: pd.DataFrame, metric: str):
    # Merge the distances
    pairwise = distance_file.merge(plasmid_communities, how="left", left_on="plasmid_1", right_on="plasmid", suffixes=("", "_1"))
    pairwise = pairwise.drop(columns="plasmid")
    pairwise = pairwise.merge(plasmid_communities, how="left", left_on="plasmid_2", right_on="plasmid", suffixes=("", "_2"))
    pairwise = pairwise.drop(columns="plasmid")

    # Locate the hub plasmids
    pairwise_communities = pairwise.dropna()

    # Limit to same (sub)communitiy plasmids
    df_pling = pairwise_communities[pairwise_communities['type'] == pairwise_communities['type_2']].copy()
    df_community = pairwise_communities[pairwise_communities['community'] == pairwise_communities['community_2']].copy()

    save_database_files(df_pling, "type", "type_2", metric)
    save_database_files(df_community, "community", "community_2", metric)


generete_database(distance_info, plasmid_communities, "distance")