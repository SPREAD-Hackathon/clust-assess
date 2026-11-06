import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
import numpy as np


pling_path = "accio_db_c5s0_subset_noclust/pling"
mash_path = "mash_distances/c5s0"


pling_dcj_df = pd.read_csv(f"{pling_path}/pling_dcj_dists_by_length.tsv", sep="\t").rename(columns={"distance": "dcj"})
pling_contain_df = pd.read_csv(f"{pling_path}/containment/all_pairs_containment_distance.tsv", sep="\t").rename(columns={"distance": "containment"})
pling_communities_df = pd.read_csv(f"{pling_path}/pling_dcj_length_subcommunities/objects/typing.tsv", sep="\t").rename(columns={"type": "type_1"})
mash_df = pd.read_csv(f"{mash_path}/plasmid_mash_distances.tsv", sep="\t")

pling_df = pling_contain_df.merge(pling_dcj_df, how="left", left_on=["plasmid_1", "plasmid_2"], right_on=["plasmid_1", "plasmid_2"])
pling_df = pling_df.merge(pling_communities_df, how="left", left_on=["plasmid_1"], right_on=["plasmid"])
pling_communities_df = pling_communities_df.rename(columns={"type_1": "type_2"})
pling_df = pling_df.merge(pling_communities_df, how="left", left_on=["plasmid_2"], right_on=["plasmid"])

all_distances = pling_df.merge(mash_df, how="left", left_on=["plasmid_1", "plasmid_2"], right_on=["plasmid A", "plasmid B"])
all_distances = all_distances.drop(columns=["Unnamed: 0", "plasmid_x", "plasmid_y", "plasmid A", "plasmid B"])
all_distances.loc[all_distances["type_1"] != all_distances["type_2"], ["dcj", "containment", "mash"]] = float("nan")

def convert_to_z(df, column_out, column_in):
    df[column_out] = (
        df[column_in] - df[column_in].mean()
    ) / df[column_in].std()
    return df

all_distances = convert_to_z(all_distances, "containment_z", "containment")
all_distances = convert_to_z(all_distances, "dcj_z", "dcj")
all_distances = convert_to_z(all_distances, "mash_z", "mash")

def single_linkage_clusters(df, distance_column, threshold=-2):
    # Get all unique plasmids
    plasmids = pd.unique(
        pd.concat([df["plasmid_1"], df["plasmid_2"]])
    )

    plasmid_to_idx = {p: i for i, p in enumerate(plasmids)}

    # Keep only significantly low Z-scores
    significant = df.loc[
        df[distance_column] <= threshold,
        ["plasmid_1", "plasmid_2"]
    ].dropna()

    # Construct an undirected adjacency matrix
    rows = significant["plasmid_1"].map(plasmid_to_idx).to_numpy()
    cols = significant["plasmid_2"].map(plasmid_to_idx).to_numpy()

    adjacency = coo_matrix(
        (
            np.ones(2 * len(rows)),
            (np.concatenate([rows, cols]), np.concatenate([cols, rows]))
        ),
        shape=(len(plasmids), len(plasmids))
    ).tocsr()

    # Connected components = single-linkage clusters
    n_clusters, labels = connected_components(
        adjacency,
        directed=False
    )

    return pd.DataFrame({
        "plasmid": plasmids,
        f"{distance_column}_cluster": labels + 1
    })


# Perform clustering independently for each distance measure
containment_clusters = single_linkage_clusters(
    all_distances, "containment_z", threshold=-2
)

dcj_clusters = single_linkage_clusters(
    all_distances, "dcj_z", threshold=-2
)

mash_clusters = single_linkage_clusters(
    all_distances, "mash_z", threshold=-2
)


# Perform clustering
containment_clusters = single_linkage_clusters(
    all_distances, "containment_z", threshold=-2
)

dcj_clusters = single_linkage_clusters(
    all_distances, "dcj_z", threshold=-2
)

mash_clusters = single_linkage_clusters(
    all_distances, "mash_z", threshold=-2
)

# Merge cluster assignments into all_distances
for cluster_df in [containment_clusters, dcj_clusters, mash_clusters]:
    all_distances = all_distances.merge(
        cluster_df,
        left_on="plasmid_1",
        right_on="plasmid",
        how="left"
    ).drop(columns="plasmid")

all_distances.to_csv("c5s0_out.csv")
