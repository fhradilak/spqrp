from importlib.resources import files
from typing import Literal
import itertools
import pandas as pd
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import networkx as nx
import matplotlib.lines as mlines

from .helpers import percentile_cutoff, plot_distribution_with_highlights
from .helpers import get_distances, get_nearest_neighbours, get_evaluation_metrics, check_input_data_format
from .helpers import get_sample_relations_by_cutoff
from .clustering_helpers import *

   
def perform_distance_evaluation_on_ranked_proteins(
    df,
    top_importance_path=  files("spqrp").joinpath("data/ranked_classification_importance_cohort_a.csv"),
    top_importance_df = None,
    n=10,
    p=0.5,
    remove_list=None,
    metric = "correlation",
    fractional_p=0.5,
    quiet = False,
    number_display_neighbours = 4,
    name = ""):
    """
    Return .
    For given n and parameters, runs the distance metric and returns the results in a list.
    Visualizes the distance distribution, fp, fn samples and the cutoff.
    

    Parameters:
        df: DataFrame containing protein intensity data.
        top_importance_path: Path to ranked protein importance data.
        n: Number of proteins used for the distance metric.
        p: Cutoff for distance metric evaluation.
        remove_list: Proteins to not include into the distance metric.
        metric: Similarity metric to use.
        range: Range of values to iterate over.
        fractional_p: Fractional value for fractional metric.
        number_display_neighbours: Number of nearest neighbours to return
    """
    # 1: retrieve ranking
    if top_importance_df is not None:
        top_importance = top_importance_df
    elif top_importance_path is not None:
        top_importance = pd.read_csv(top_importance_path)
    else:
        raise ValueError("Provide either top_importance_path or top_importance_df")

    top_importance = top_importance.sort_values(by="Importance", ascending=False)[:n]
    if remove_list is not None:
        top_importance = top_importance[~top_importance["Protein"].isin(remove_list)]
    
    # check for correct input format
    check_input_data_format(df,top_importance)
    # filter required proteins
    df_dist = df[df["Protein"].isin(top_importance["Protein"])]
    
    # 2: calculate distances and nearest neighbours
    distance_matrix, df_pivot = get_distances(df_dist, metric = metric, fractional_p=fractional_p)
    sample_order = list(df_pivot.index)
    
    available = len(df_dist["Sample_ID"].unique()) - 1  # max possible neighbors
    k = min(number_display_neighbours, available)
    nearest_neighbours = get_nearest_neighbours(df_pivot,distance_matrix, k)
    
    # 3: calculate percentile cutoff
    percentiles = [p]
    # Flatten the upper triangle of the distance matrix (excluding diagonal)
    distances = distance_matrix[np.triu_indices(len(distance_matrix), k=1)]
    cutoff = percentile_cutoff(distances=distances, percentile=p)
    
    # 4: Classify sample pairings as belonging or not belonging
    sample_patient_mapping = dict(df_dist[["Sample_ID","Patient_ID"]].drop_duplicates().values)
    print("Real Number Proteins", len(df_dist["Protein"].unique()))
    results = get_sample_relations_by_cutoff(distance_matrix,cutoff=cutoff, sample_patient_mapping=sample_patient_mapping, sample_order=sample_order)
    belonging, not_belonging = results["belonging"], results["not_belonging"]
    
    # 5: Evaluate Performance
    eval_metrics = get_evaluation_metrics(belonging, not_belonging, quiet = quiet)
    plot_distribution_with_highlights(distances,  eval_metrics["False_Negative_Distances"], eval_metrics["False_Positive_Distances"], percentiles=percentiles, name = name + f" {n} proteins")
    
    return {
        "top_importance": top_importance,
        "nearest_neighbours": nearest_neighbours,
        "cutoff": cutoff,
        "belonging": belonging,
        "not_be": not_belonging,
        "eval_metrics": eval_metrics,
        "distance_matrix": distance_matrix
    }
    

def is_better_result(result, strategy, thresholds):
    """
    Determine if the current `result` is better than given `thresholds` based on the chosen `strategy`.

    Args:
        result (dict): 
            Dictionary containing performance metrics with keys:
            - "FP" (False Positives)
            - "FN" (False Negatives)
            - "Precision"
            - "Sensitivity"
            - optionally "F1"
        strategy (str): 
            Strategy to evaluate the result. Supported values:
            - "fp+fn": compare sum of FP+FN, then Precision
            - "fp": compare FP, then sum FP+FN
            - "fn": compare FN, then sum FP+FN
            - "F1": compare F1 score
            - "sensitivity": compare Sensitivity
            - "precision": compare Precision
        thresholds (dict): 
            Threshold values for the metrics, e.g.:
            - "fp+fn", "fp", "fn", "f1", "sensitivity", "precision"

    Returns:
        bool: True if the current result is better according to the strategy, False otherwise.
    """
    fp, fn, precision, f1, sensitivity = result["FP"], result["FN"], result["Precision"], result["Sensitivity"], result.get("F1", 0)

    if strategy == "fp+fn":
        curr_sum = fp + fn
        return (curr_sum < thresholds["fp+fn"]
                or (curr_sum == thresholds["fp+fn"] and precision > thresholds["f1"]))

    if strategy == "fp":
        return (fp < thresholds["fp"]
                or (fp == thresholds["fp"] and (fp + fn) < thresholds["fp+fn"]))

    if strategy == "fn":
        return (fn < thresholds["fn"]
                or (fn == thresholds["fn"] and (fp + fn) < thresholds["fp+fn"]))

    if strategy == "F1":
        return f1 > thresholds["f1"]
    
    if strategy == "sensitivity":
        return sensitivity > thresholds["sensitivity"]
    
    if strategy == "precision":
        return precision > thresholds["precision"]
    

    return False

def optimize_parameters(
    df,
    metric="correlation",
    log_file="optimization_log.txt",
    top_importance_path= files("spqrp").joinpath("data/ranked_classification_importance_cohort_a.csv"),
    top_importance_df = None,
    range = range(2,50,1),
    optimization_strategy: Literal["fp+fn", "fp", "fn", "F1", "precision", "sensitivity"] = "fp+fn",
    remove_list=[],
    quiet = False):
    """
    Optimize for a distance metric a specified performance metric(F1, precision, sensitivity, #fn, #fp, #fp+fn) for each n in a range.
    Finds the best parameter(s) e.g. cutoff for a given n and returns the best results in a list.

    Parameters:
        df: DataFrame containing input data.
        metric: Similarity metric to use.
        log_file: File to write optimization logs.
        top_importance_path: Path to ranked protein importance data.
        range: Range of values to iterate over.
        optimization_strategy: Optimization objective:
            - "fp+fn": Minimize false positives + false negatives 
            - "fp":   Minimize false positives
            - "fn":   Minimize false negatives
            - "F1":   Maximize F1 score
    """
    
    fractional = metric == ("fractional")
    best_results = []  
    n_values = [int(v) for v in range]
    fractional_p_values = np.arange(0.1, 1.0, 0.01) if fractional else [None]
    percentile_values = [25, 10,]
    descending = list(np.arange(5.0, 1.0, -0.1)) 
    ascending = list(np.arange(0.1, 1.0, 0.01))
    percentile_values += descending + ascending

    with open(log_file, "w") as log:
        log.write("Starting parameter optimization...\n\n")

         if top_importance_df is not None:
             classifier_top_proteins = top_importance_df
         elif top_importance_path is not None:
             classifier_top_proteins = pd.read_csv(top_importance_path)
         else:
             raise ValueError("Provide either top_importance_path or top_importance_df")
            
        for n in n_values:
            
            #set up initial values for optimization per n
            best_result_for_n = None
            best_params_for_n = None
            thresholds = {
                "fp+fn": float("inf"),
                "fp": float("inf"),
                "fn": float("inf"),
                "f1": 0,
                "precision": 0,
                "sensitivity": 0,
            }
            
            # 1: Get protein ranking and Sample-Patient Mapping
           
             top_importance = classifier_top_proteins
             top_importance = top_importance[~top_importance["Protein"].isin(remove_list)]
             
             ### Get only overlapping proteins 
             # important since we only have a sub overlapp with the proteom
             top_prot = top_importance[top_importance["Protein"].isin(df["Protein"])].sort_values(by="Importance", ascending=False).head(n)

            # 2: Sanity checking
            check_input_data_format(df,top_importance)
            
            df_dist = df[df["Protein"].isin(top_prot["Protein"])]
            sample_patient_mapping = df_dist.drop_duplicates(subset=["Sample_ID", "Patient_ID"]).set_index("Sample_ID")["Patient_ID"].to_dict()
            
            
            for fractional_p, percentile in itertools.product(fractional_p_values, percentile_values):
                # 3: Compute distances & neighbors & get sample order
                distance_matrix, df_pivot = get_distances(df_dist, metric=metric, fractional_p=fractional_p)
                sample_order = list(df_pivot.index)
                log.write(f"Testing n={n}, fractional_p={fractional_p}, percentile={percentile}\n")
                
                # 4: Determine cutoff
                # Flatten the upper triangle of the distance matrix (excluding diagonal)
                distances = distance_matrix[np.triu_indices(len(distance_matrix), k=1)]
                cutoff = percentile_cutoff(distances=distances, percentile=percentile)
                
                # 5: Identify belonging and non-belonging samples
                results = get_sample_relations_by_cutoff(distance_matrix,cutoff=cutoff, sample_patient_mapping=sample_patient_mapping, sample_order=sample_order)
                belonging, not_belonging = results["belonging"], results["not_belonging"]
            
                # 6: Evaluate results
                result = get_evaluation_metrics(belonging, not_belonging, quiet=quiet)
                # Log current results
                log.write(f"Precision: {result['Precision']}, Sensitivity: {result['Sensitivity']}, TP:{result['TP']}, FP:{result['FP']}, TN:{result['TN']}, FN:{result['FN']}\n\n")

                #7: compare current values to best so far
                if is_better_result(result, optimization_strategy, thresholds):
                    thresholds.update({
                        "fp+fn": result["FP"] + result["FN"],
                        "fp": result["FP"],
                        "fn": result["FN"],
                        "f1": result.get("F1", 0),
                        "precision": result["Precision"],
                        "sensitivity":result["Sensitivity"]
                    })

                    best_result_for_n = result
                    best_params_for_n = (n, fractional_p, percentile)
                
                if not quiet:
                    print(f"Testing n={n}, fractional_p={fractional_p}, percentile={percentile}")
                    print("Real Number Proteins", len(df_dist["Protein"].unique()))
                    print(distance_matrix.size)
                    print("Cutoff", cutoff)
                    print(result["FP"] + result["FN"] + result["TP"] + result["TN"])

            # 8: Store best result for this n
            if best_result_for_n:
                best_results.append({
                    "n": best_params_for_n[0],
                    "fractional_p": best_params_for_n[1],
                    "percentile": best_params_for_n[2],
                    "FP": best_result_for_n["FP"],
                    "FN": best_result_for_n["FN"],
                    "TP": best_result_for_n["TP"],
                    "TN": best_result_for_n["TN"],
                    "Precision": best_result_for_n["Precision"],
                    "Sensitivity": best_result_for_n["Sensitivity"],
                    "F1":best_result_for_n["F1"],
                    "Protein": top_importance["Protein"].iloc[n-1]
                    
                })

        # Convert results to a DataFrame
        results_df = pd.DataFrame(best_results)
    return results_df


# =========== CLUSTERING ====================

def cluster_samples_iteratively(
    result,
    df,
    method,
    random_state=42,
    n_neighbors= 1,
    max_component_size=None,
    n_umap_neighbors=15,
    precomputed_graph = None,
):
    # create distance matrics & mappings 
    n_neighbors = int(n_neighbors)
    dist_matrix = result["distance_matrix"]
    sample_names = sorted(df['Sample_ID'].unique())
    sample_index = {name: i for i, name in enumerate(sample_names)}
    n_samples = len(sample_names)

    G,coords_2d = create_graph_based_on_reduction_method(method,n_samples,dist_matrix,sample_names, n_umap_neighbors, random_state, precomputed_graph, sample_index,n_neighbors)
        
    if max_component_size is not None:
        G = split_big_component_edges_by_weight(G, max_component_size)
    return G, coords_2d

def plot_distances_neighbours_with_coloring_hue(
    df,
    G,
    coords_2d,
    method='UMAP',
    subset_samples=None,
    highlight_singletons=True,
    highlight_single_samples_missing_connections=True,
    figsize=(20,20),
    label_patient_only=False,
    label_offset_x = 0.01,
    label_offset_y = 0.01,
    label_font = 6,
    df_name = "DF_NAME",
    return_clusters = False,
):
    fp_color = "#DC267F"#"#4B0092"
    tp_color="#1AFF1A"
    error_candidate_color ="#FFBBF6"# "#FFB000"#CC6677"
    
    sample_to_patient = dict(zip(df['Sample_ID'], df['Patient_ID']))
    samples_by_patient = defaultdict(list)
    for sample, pid in sample_to_patient.items():
        samples_by_patient[pid].append(sample)
    drawn_pairs = {tuple(sorted(e)) for e in G.edges()}
    connected_samples = set(G.nodes())
    sample_names = sorted(df['Sample_ID'].unique())
    sample_index = {name: i for i, name in enumerate(sample_names)}
    # === Define subset to be visualized ===
    if subset_samples:
        connected_samples &= set(subset_samples)
        
    # === Create Plot ===
    fig, ax = plt.subplots(figsize=figsize, dpi=150)
    # 1 === Draw TP clusters ===
    draw_convex_hull(samples_by_patient, connected_samples, ax, drawn_pairs, sample_index, coords_2d, tp_color)
    
    # 2 === Draw red edges for FP - error candidates===
    fp_linestyle = ":"  # dashed pattern
    for s1, s2 in drawn_pairs:
        if any(s in connected_samples for s in (s1, s2)):
            i, j = sample_index[s1], sample_index[s2]
            if sample_to_patient[s1] != sample_to_patient[s2]:
                x1, y1 = coords_2d[i]
                x2, y2 = coords_2d[j]
                ax.plot([x1, x2], [y1, y2],
                        color=fp_color,
                        linewidth=2,
                        linestyle=fp_linestyle)
    
    # === Retrieve clustering type for nodes  ===
    nodes_in_tp_clusters, nodes_with_red_connections, singleton_nodes, isolated_nodes = identify_clusters_singletons(G,sample_to_patient,samples_by_patient, drawn_pairs, sample_names)
    
    # === Draw nodes ===
    for name in connected_samples:
        if name not in sample_index:
            continue
        i = sample_index[name]
        x, y = coords_2d[i]

        # === Determine node color ===
        if (name in singleton_nodes) or (name in nodes_in_tp_clusters):
            color = tp_color
        elif name in nodes_with_red_connections:
                color = fp_color
        elif name in isolated_nodes:
                color = error_candidate_color
        # === Draw node ===
        plt.scatter(x, y, color=color, s=20)
        
        # Label control
        if label_patient_only:
            label_text = sample_to_patient[name]
        else:
            label_text = name
        plt.text(x+label_offset_x, y+label_offset_y, label_text, fontsize=label_font)

        # === Blue rectangle for singleton nodes ===
        if (highlight_singletons and 
        name in singleton_nodes):
            size = 0.05
            square = patches.Rectangle(
                (x - size/2, y - size/2),
                size, size,
                edgecolor='blue',
                facecolor='none',
                linewidth=1
            )
            ax.add_patch(square)

        # === Pink circle for disconnected-but-should-be-connected ===
        if highlight_single_samples_missing_connections and name in isolated_nodes:
            circle = plt.Circle((x, y), 0.05, edgecolor=error_candidate_color, facecolor='none', linewidth=2)
            ax.add_patch(circle)

    # === Metrics ===    
    transitive_results = transitive_performance(sample_names, drawn_pairs, sample_to_patient)
    
    # === Legends ===
    handles = []
    fp_edge_handle = mlines.Line2D([], [], color=fp_color,
                                linestyle=fp_linestyle,
                                linewidth=2,
                                label='Error Candidate')
    tp_edge_handle = mpatches.Patch(color=tp_color, label='True Positive Cluster')
    handles.append(fp_edge_handle)
    handles.append(tp_edge_handle)
    
    
    if highlight_singletons:
        square_legend = mpatches.Rectangle(
            (0, 0), 1, 1,
            edgecolor='blue',
            facecolor='none',
            linewidth=2,
            label='True Positive Singular Sample'
        )
        handles.append(square_legend)
    if highlight_single_samples_missing_connections:
        uncertain_handle = mlines.Line2D([], [], 
        color=error_candidate_color,
        marker='o',
        linestyle='None',
        markersize=10,
        markerfacecolor='none',
        label='Uncertain Sample'
        )
        handles.append(uncertain_handle)

    plt.title(f'{method} projection of {df_name} sample Clustering\n', fontsize=16) #Green = TP cluster | Red = FP node | Violet = not connected but same patient
    plt.axis('equal')
    plt.grid(True)
    plt.legend(handles=handles, loc="best")
    if method.upper() == 'PCA':
        ax.set_xlabel("PCA 1", fontsize=14)
        ax.set_ylabel("PCA 2", fontsize=14)
    elif method.upper() == 'MDS':
        ax.set_xlabel("MDS Dimension 1", fontsize=14)
        ax.set_ylabel("MDS Dimension 2", fontsize=14)
    elif method.upper() =="UMAP":
        ax.set_xlabel(f"{method} Dim 1", fontsize=14)
        ax.set_ylabel(f"{method} Dim 2", fontsize=14)
    
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig("distance_plot.svg")
    plt.show()
    
    cluster_assignments = None
    if return_clusters:
        components = list(nx.connected_components(G))
        cluster_assignments = {}
        for cluster_id, component in enumerate(components, start=1):
            for sample in component:
                cluster_assignments[sample] = cluster_id

        # Print clusters in readable form
        print("\n=== Cluster Assignments ===")
        for cluster_id, component in enumerate(components, start=1):
            print(f"Cluster {cluster_id} ({len(component)} samples): {sorted(component)}")

        return G, cluster_assignments, transitive_results 
    else:
        return G



