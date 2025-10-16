from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn.metrics import (
    precision_score,
    recall_score,
    f1_score,
    accuracy_score,
    balanced_accuracy_score,
    confusion_matrix,
    adjusted_rand_score,
    normalized_mutual_info_score,
)
from scipy.spatial import ConvexHull
import matplotlib.patches as patches
import networkx as nx
import matplotlib.patches as mpatches
import matplotlib.lines as mlines


# ======== Graph Creation & Clustering by Splitting ==========
def create_graph_based_on_reduction_method(
    method,
    n_samples,
    dist_matrix,
    sample_names,
    n_umap_neighbors,
    random_state,
    precomputed_graph,
    sample_index,
    n_neighbors,
):
    if method == "PCA":
        J = np.eye(n_samples) - np.ones((n_samples, n_samples)) / n_samples
        B = -0.5 * J @ (dist_matrix**2) @ J
        coords_2d = PCA(n_components=2).fit_transform(B)
    elif method == "UMAP":
        import umap

        reducer = umap.UMAP(
            n_components=2,
            n_neighbors=n_umap_neighbors,
            metric="precomputed",
            random_state=random_state,
        )
        coords_2d = reducer.fit_transform(dist_matrix)
        method = "UMAP"
    elif method == "MDS":
        coords_2d = MDS(
            n_components=2, dissimilarity="precomputed", random_state=random_state
        ).fit_transform(dist_matrix)

    G = nx.Graph()
    for sample in sample_names:
        G.add_node(sample)

    drawn_pairs = set()
    connected_samples = set()

    if precomputed_graph is not None:
        for u, v in precomputed_graph.edges():
            if u in sample_index and v in sample_index:
                d = dist_matrix[sample_index[u], sample_index[v]]
                G.add_edge(u, v, weight=d)
                drawn_pairs.add(tuple(sorted((u, v))))
                connected_samples.update([u, v])
    else:
        for i in range(n_samples):
            sorted_idx = np.argsort(dist_matrix[i])[1 : n_neighbors + 1]
            for j in sorted_idx:
                s1, s2 = sample_names[i], sample_names[j]
                dist = dist_matrix[i, j]
                pair = tuple(sorted((s1, s2)))
                G.add_edge(s1, s2, weight=dist)
                drawn_pairs.add(pair)
                connected_samples.update([s1, s2])
    return G, coords_2d


def split_big_component_edges_by_weight(G, max_size):
    """
    Removes edges with the highest weights from large components
    until all components are smaller than max_size.
    """
    G = G.copy()
    changed = True
    while changed:
        changed = False
        for comp in list(nx.connected_components(G)):
            if len(comp) <= max_size:
                continue
            sub = G.subgraph(comp).copy()
            # Get edge with maximum weight (i.e., weakest connection)
            edge_weights = [
                (e, sub[e[0]][e[1]].get("weight", 1.0)) for e in sub.edges()
            ]
            # Sort descending by weight
            edge_to_remove = max(edge_weights, key=lambda x: x[1])[0]
            G.remove_edge(*edge_to_remove)
            changed = True
            break  # Only remove one per iteration
    return G


# ========== Drawing Helpers ==========


def identify_clusters_singletons(
    G, sample_to_patient, samples_by_patient, drawn_pairs, sample_names
):
    # === FP nodes and singleton identification ===
    nodes_in_fp_cluster = {
        s
        for s1, s2 in drawn_pairs
        if sample_to_patient[s1] != sample_to_patient[s2]
        for s in (s1, s2)
    }
    singleton_nodes = {
        s for s in sample_names if len(samples_by_patient[sample_to_patient[s]]) == 1
    }
    singleton_graph = {s for s in singleton_nodes if G.degree(s) == 0}
    # === TP nodes identification ===
    nodes_in_tp_clusters = set()
    for component in nx.connected_components(G):
        component = list(component)
        if len(component) < 2:
            continue
        patients = {sample_to_patient[s] for s in component}
        if len(patients) == 1:
            pid = next(iter(patients))
            if len(samples_by_patient[pid]) >= 2:
                nodes_in_tp_clusters.update(component)

    # ===  Isolated nodes
    isolated_nodes = set()
    for s in sample_names:
        if (
            s in singleton_nodes
            or s in nodes_in_fp_cluster
            or s in nodes_in_tp_clusters
        ):
            continue
        pid = sample_to_patient[s]
        patient_mates = set(samples_by_patient[pid]) - {s}
        if not any(G.has_edge(s, mate) for mate in patient_mates if mate in G):
            isolated_nodes.add(s)

    return nodes_in_tp_clusters, nodes_in_fp_cluster, singleton_graph, isolated_nodes


def draw_convex_hull(
    samples_by_patient,
    connected_samples,
    ax,
    drawn_pairs,
    sample_index,
    coords_2d,
    tp_color,
):
    # === Draw convex hulls around TP clusters ===
    for pid, sample_list in samples_by_patient.items():
        filtered_samples = [
            s
            for s in sample_list
            if s in connected_samples
            and any(
                tuple(sorted((s, other))) in drawn_pairs
                for other in sample_list
                if s != other
            )
        ]

        if len(filtered_samples) < 2:
            continue

        if len(filtered_samples) == 2:
            s1, s2 = filtered_samples
            pair = tuple(sorted((s1, s2)))
            if pair in drawn_pairs:
                i, j = sample_index[s1], sample_index[s2]
                ax.plot(*zip(coords_2d[i], coords_2d[j]), color=tp_color, linewidth=3)
            continue

        coords = np.array([coords_2d[sample_index[s]] for s in filtered_samples])
        try:
            hull = ConvexHull(coords)
            vertices = coords[hull.vertices]
            patch = patches.Polygon(
                vertices,
                closed=True,
                facecolor=tp_color,
                edgecolor=(0.0, 0.4, 0.0, 0.5),
                linewidth=1,
            )
            ax.add_patch(patch)
        except Exception as e:
            print(f"Could not compute hull for patient {pid}: {e}")


# ========== Performance ==========
def transitive_performance(sample_names, drawn_pairs, sample_to_patient):
    # === Build Transitive Clusters (Connected Components) ===
    G = nx.Graph()
    G.add_nodes_from(sample_names)
    G.add_edges_from(drawn_pairs)

    # Cluster assignment: {sample: cluster_id}
    pred_clusters = {
        sample: i
        for i, comp in enumerate(nx.connected_components(G))
        for sample in comp
    }
    true_clusters = {sample: sample_to_patient[sample] for sample in sample_names}

    # === Pairwise Evaluation using Transitive Clusters ===
    y_true, y_pred = [], []
    false_negatives = []  # Store FN pairs here

    for s1, s2 in combinations(sample_names, 2):
        same_patient = true_clusters[s1] == true_clusters[s2]
        same_cluster = pred_clusters[s1] == pred_clusters[s2]
        y_true.append(same_patient)
        y_pred.append(same_cluster)

        # Collect FN pairs
        if same_patient and not same_cluster:
            false_negatives.append((s1, s2))

    cm = confusion_matrix(y_true, y_pred, labels=[True, False])
    TP = cm[0, 0]
    FN = cm[0, 1]
    FP = cm[1, 0]
    TN = cm[1, 1]

    precision = precision_score(y_true, y_pred)
    sensitivity = recall_score(y_true, y_pred)
    f1 = f1_score(y_true, y_pred)
    accuracy = accuracy_score(y_true, y_pred)
    bACC = balanced_accuracy_score(y_true, y_pred)

    print("=== Pairwise Clustering Performance (Transitive) ===")
    print(f"TP: {TP}, FP: {FP}, FN: {FN}, TN: {TN}")
    print(f"Precision: {precision:.3f}")
    print(f"Sensitivity:    {sensitivity:.3f}")
    print(f"F1 Score:  {f1:.3f}")
    print(f"Accuracy:  {accuracy:.3f}")
    print(f"Balanced Accuracy:  {bACC:.3f}")

    # === Clustering Metrics (ARI, NMI) ===
    y_true_cluster = [true_clusters[s] for s in sample_names]
    y_pred_cluster = [pred_clusters[s] for s in sample_names]

    ari = adjusted_rand_score(y_true_cluster, y_pred_cluster)
    nmi = normalized_mutual_info_score(y_true_cluster, y_pred_cluster)

    print("=== Overall Clustering Agreement ===")
    print(f"Adjusted Rand Index (ARI): {ari:.3f}")
    print(f"Normalized Mutual Info (NMI): {nmi:.3f}")

    # === Print False Negative Pairs ===
    if false_negatives:
        print(
            f"\nFalse Negative (FN) pairs (same patient but not transitively connected): {len(false_negatives)}"
        )
        for s1, s2 in false_negatives:
            print(f"  - {s1} <-> {s2}")
    else:
        print("\nNo False Negative (FN) pairs found.")

    return {
        "precision": precision,
        "sensitivity": sensitivity,
        "f1": f1,
        "accuracy": accuracy,
        "balanced accuracy": bACC,
        "ari": ari,
        "nmi": nmi,
    }
