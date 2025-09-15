import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from collections import defaultdict
import seaborn as sns
from itertools import combinations


def get_connected_samples(
    n_neighbour_in_belonging: bool,
    result: dict,
    sample_to_patient: dict,
    sample_index: dict,
    dist_matrix: np.ndarray,
    n_neighbors: int,
    sample_names: list,
    samples_by_patient: dict,
    subset_samples: list
) -> set:
    """
    Identify samples connected to a given subset either by belonging or distance criteria.

    Args:
        n_neighbour_in_belonging (bool): 
            If True, use belonging/not_belonging info from `result` to find connected samples.
            If False, use distance matrix and nearest neighbors.
        result (dict): 
            Output from distance evaluation containing "belonging" and "not_be" lists of sample pairs.
        sample_to_patient (dict): 
            Maps sample names to their patient IDs.
        sample_index (dict): 
            Maps sample names to their index in `dist_matrix` and `sample_names`.
        dist_matrix (np.ndarray): 
            Distance matrix between samples.
        n_neighbors (int): 
            Number of neighbors to consider per sample.
        sample_names (list): 
            List of all sample names corresponding to indices in `dist_matrix`.
        samples_by_patient (dict): 
            Maps patient IDs to lists of sample names belonging to that patient.
        subset_samples (list): 
            Subset of sample names for which connected samples are sought.

    Returns:
        set: 
            Set of sample names connected to the subset either by belonging or nearest neighbors.
    """
    subset_samples = set(subset_samples)
    connected_samples = set(subset_samples)

    if n_neighbour_in_belonging:
        for d in result["belonging"]:
            if d['sample1'] in subset_samples or d['sample2'] in subset_samples:
                connected_samples.update([d['sample1'], d['sample2']])

        for d in result["not_be"]:
            s1, s2 = d['sample1'], d['sample2']
            if (s1 in subset_samples or s2 in subset_samples) and sample_to_patient[s1] == sample_to_patient[s2]:
                connected_samples.update([s1, s2])

    else:
        one_hop = set()
        for name in subset_samples:
            if name in sample_index:
                i = sample_index[name]
                sorted_idx = np.argsort(dist_matrix[i])[1:n_neighbors+1]
                one_hop.update(sample_names[j] for j in sorted_idx)

        for j, other_name in enumerate(sample_names):
            if other_name in subset_samples:
                continue
            sorted_idx = np.argsort(dist_matrix[j])[1:n_neighbors+1]
            if any(sample_names[k] in subset_samples for k in sorted_idx):
                one_hop.add(other_name)

        connected_samples.update(one_hop)

    # Add violet partners: same-patient pairs where only one is in subset
    extra_violet_samples = set()
    for samples in samples_by_patient.values():
        for s1, s2 in combinations(samples, 2):
            if (s1 in subset_samples) != (s2 in subset_samples):
                extra_violet_samples.update([s1, s2])
    connected_samples.update(extra_violet_samples)

    return connected_samples
        


def plot_distances_neighbours_with_coloring(result,df, method = 'mds', n_neighbour_in_belonging= False, n_neighbors=1, subset_samples=None, draw_fn = True):
    """
    Visualizes pairwise sample distances in 2D space using PCA or MDS, with colored edges to indicate relationships 
    between samples based on patient ID and neighborhood relationships.

    Parameters:
    -----------
    result : dict
        Dictionary containing two lists:
        - 'belonging': list of dicts, each with keys 'sample1', 'sample2', and 'distance' (pairs considered belonging).
        - 'not_be': list of dicts with similar structure (pairs not considered belonging).
    
    df : pd.DataFrame
        DataFrame containing 'Sample_ID' and 'Patient_ID' columns, used to map samples to patients.

    pca : bool
        Whether to use PCA (True) or MDS (False) for dimensionality reduction and 2D projection.

    n_neighbour_in_belonging : bool
        If True, use the result sets 'belonging' and 'not_be' as basis for the coloring.
        If False, use the nearest neighbours as basis for the coloring.

    n_neighbors : int, optional (default=1)
        Number of nearest neighbors to consider when `n_neighbour_in_belonging` is False.

    subset_samples : list of str, optional
        List of sample IDs to focus the visualization on. If provided, the plot only includes these
        samples and their neighbors (depending on settings).

    Visual Elements & Coloring:
    ----------------
    - **Nodes**: Each point represents a sample, colored by its patient ID.
    - **Edges**:
        - Green: Same-patient pairs that are in 'belonging' or nearest neighbors.
        - Red: Cross-patient pairs that are in 'belonging' or nearest neighbors.
        - Violet: Same-patient pairs not in 'belonging' or not nearest neighbors.
    """
    n_neighbors = int(n_neighbors)  # Ensure Python integer type for slicing (for R)
    
    all_distances = result["belonging"] + result["not_be"]
    sample_to_patient = dict(zip(df['Sample_ID'], df['Patient_ID']))

    samples_by_patient = defaultdict(list)
    for sample, pid in sample_to_patient.items():
        samples_by_patient[pid].append(sample)

    sample_names = sorted(set(sample_to_patient.keys()))
    sample_index = {name: i for i, name in enumerate(sample_names)}
    n_samples = len(sample_names)

    # Build distance matrix
    dist_matrix = np.full((n_samples, n_samples), np.nan)
    for d in all_distances:
        i, j = sample_index[d['sample1']], sample_index[d['sample2']]
        dist_matrix[i, j] = dist_matrix[j, i] = float(d['distance'])
    np.fill_diagonal(dist_matrix, 0)
    dist_matrix = np.nan_to_num(dist_matrix, nan=1.0)

    # Determine which samples to display
    if subset_samples:
       connected_samples= get_connected_samples(n_neighbour_in_belonging, result, sample_to_patient, sample_index, dist_matrix, n_neighbors, sample_names,samples_by_patient, subset_samples)
    else:
        connected_samples = set(sample_names)
        subset_samples = set(sample_names)
        
    # Create 2D coordinates from distances
    if method =='pca':
        J = np.eye(n_samples) - np.ones((n_samples, n_samples)) / n_samples
        B = -0.5 * J @ (dist_matrix ** 2) @ J
        coords_2d = PCA(n_components=2).fit_transform(B)
    elif method == 'umap':
        import umap
        reducer = umap.UMAP(n_components=2, metric='precomputed', random_state=42)
        coords_2d = reducer.fit_transform(dist_matrix)
        method = 'UMAP'
    else:

        coords_2d = MDS(n_components=2, dissimilarity='precomputed', random_state=42).fit_transform(dist_matrix)

    # Plot points
    plt.figure(figsize=(20, 20))
    unique_patients = sorted(set(sample_to_patient.values()))
    palette = sns.color_palette("tab10", n_colors=len(unique_patients))
    patient_to_color = {pid: palette[i] for i, pid in enumerate(unique_patients)}

    for name in connected_samples:
        if name not in sample_index:
            continue
        i = sample_index[name]
        x, y = coords_2d[i]
        pid = sample_to_patient[name]
        color = patient_to_color[pid]
        label = f"Patient {pid}"
        if label not in plt.gca().get_legend_handles_labels()[1]:
            plt.scatter(x, y, color=color, label=label)
        else:
            plt.scatter(x, y, color=color)
        plt.text(x, y, name, fontsize=9)

    drawn_pairs = set()

    def should_draw(s1, s2):
        return s1 in connected_samples and s2 in connected_samples

    if n_neighbour_in_belonging:
        for d in result["belonging"]:
            s1, s2 = d['sample1'], d['sample2']
            if not should_draw(s1, s2): continue
            pair = tuple(sorted((s1, s2)))
            if pair in drawn_pairs: continue
            drawn_pairs.add(pair)
            i, j = sample_index[s1], sample_index[s2]
            color = 'green' if sample_to_patient[s1] == sample_to_patient[s2] else 'red'
            plt.plot(*zip(coords_2d[i], coords_2d[j]), color=color, linewidth=2, alpha=0.6)

        for d in result["not_be"]:
            s1, s2 = d['sample1'], d['sample2']
            if sample_to_patient[s1] == sample_to_patient[s2] and should_draw(s1, s2):
                pair = tuple(sorted((s1, s2)))
                if pair not in drawn_pairs:
                    drawn_pairs.add(pair)
                    i, j = sample_index[s1], sample_index[s2]
                    plt.plot(*zip(coords_2d[i], coords_2d[j]), color='violet', linewidth=2, alpha=0.6)

    else:
        neighbors = defaultdict(set)
        for i in range(n_samples):
            sorted_idx = np.argsort(dist_matrix[i])[1:n_neighbors+1]
            for j in sorted_idx:
                neighbors[sample_names[i]].add(sample_names[j])
                neighbors[sample_names[j]].add(sample_names[i])

        for s1, s2_set in neighbors.items():
            for s2 in s2_set:
                pair = tuple(sorted((s1, s2)))
                if not should_draw(s1, s2) or pair in drawn_pairs: continue
                drawn_pairs.add(pair)
                i, j = sample_index[s1], sample_index[s2]
                color = 'green' if sample_to_patient[s1] == sample_to_patient[s2] else 'red'
                plt.plot(*zip(coords_2d[i], coords_2d[j]), color=color, linewidth=2, alpha=0.6)

    # Draw remaining same-patient edges in violet if not already drawn
    if(draw_fn):
        for samples in samples_by_patient.values():
            for s1, s2 in combinations(samples, 2):
                pair = tuple(sorted((s1, s2)))
                if pair in drawn_pairs: continue
                if should_draw(s1, s2) and (s1 in subset_samples or s2 in subset_samples):
                    i, j = sample_index[s1], sample_index[s2]
                    plt.plot(*zip(coords_2d[i], coords_2d[j]), color='violet', linewidth=1.5, alpha=0.6)

    
    if n_neighbour_in_belonging:
        legend_green = 'pairs in belonging and same patient ID'
        legend_red = 'pairs in belonging but different patient ID'
        legend_violet = 'pairs not in belonging but same patient ID'
    else:
        legend_green = 'nearest neighbor with same patient ID'
        legend_red = 'nearest neighbor with different patient ID'
        legend_violet = 'not nearest neighbor but same patient ID'

    plt.title(f'{method} Projection\n'
              f'Green = {legend_green}\n'
              f'Red = {legend_red}\n'
              f'Violet = {legend_violet}')
    plt.axis('equal')
    plt.grid(True)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig("distance_plot.svg")
    plt.show()