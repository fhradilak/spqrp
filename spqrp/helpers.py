from sklearn.metrics import pairwise_distances
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple

SAMPLE = "Sample_ID"
PATIENT = "Patient_ID"

##################################
# Input Checking
##################################
required_columns_df = ["Sample_ID", "Patient_ID", "Intensity", "Protein"]
required_columns_ranking = ["Protein", "Importance"]


def check_required_columns(df, required_columns, type):
    required_columns = required_columns
    missing_columns = [col for col in required_columns if col not in df.columns]
    print(f"Checking if columns: {required_columns} are present in {type}")
    if missing_columns:
        print(
            f"For {type} the following required columns are missing from the DataFrame: {missing_columns}. You probably have to rename the columns of the dataset to fit the requirements."
        )
        raise ValueError(
            f"For {type} the following required columns are missing from the DataFrame: {missing_columns}. You probably have to rename the columns of the dataset to fit the requirements."
        )
    else:
        print(f"All required columns are present for {type}.")


def check_input_data_format(df, importance_ranking):
    check_required_columns(
        df=df, required_columns=required_columns_df, type="protein matrix"
    )
    check_required_columns(
        df=importance_ranking,
        required_columns=required_columns_ranking,
        type="importance ranking",
    )


def percentile_cutoff(distances, percentile=25):
    return np.percentile(distances, percentile)


##################################
# Distance Calculation
##################################


def get_distances(
    df_dist,
    metric="correlation",
    intensity="Intensity",
    kernelpca=False,
    fractional_p=0.5,
    index=SAMPLE,
):

    df_pivot = df_dist.pivot_table(
        index=index, columns="Protein", values=intensity, aggfunc="mean"
    )
    df_pivot.fillna(0, inplace=True)
    fractional = metric == ("fractional")

    if fractional:
        distance_matrix = pairwise_distances(
            df_pivot, metric="minkowski", p=fractional_p
        )
    else:
        distance_matrix = pairwise_distances(df_pivot, metric=metric)

    return distance_matrix, df_pivot


def get_nearest_neighbours(df_pivot, distance_matrix, k=4):
    # check k and reassign if unsafe
    k = min(k, len(df_pivot) - 1)

    neighbor_indices = np.argsort(distance_matrix, axis=1)[:, 1 : k + 1]
    neighbor_distances = np.take_along_axis(distance_matrix, neighbor_indices, axis=1)
    neighbor_names = np.vectorize(lambda i: df_pivot.index[i])(neighbor_indices)

    # --- Create combined neighbor + distance dataframe ---
    nearest_neighbors = pd.DataFrame(index=df_pivot.index)
    for i in range(k):
        nearest_neighbors[f"Neighbor_{i+1}"] = neighbor_names[:, i]
        nearest_neighbors[f"Distance_{i+1}"] = np.round(neighbor_distances[:, i], 2)

    return nearest_neighbors


##################################
# Plotting
##################################


def plot_distribution_of_pairwise_dist(distances, percentiles=[1, 2, 5, 10]):
    plt.figure(figsize=(8, 5))
    plt.hist(distances, bins=50, edgecolor="black", alpha=0.7)
    plt.xlabel("Pairwise Distance")
    plt.ylabel("Frequency")
    plt.title("Distribution of Pairwise Distances")
    for p in percentiles:
        plt.axvline(
            np.percentile(distances, p),
            color="#DC267F",
            linestyle="dashed",
            linewidth=2,
            label=f"{p}th Percentile",
        )
    plt.legend()
    plt.show()


def plot_distribution_with_highlights(
    distances, fn_distances, fp_distances, percentiles=[1, 2, 5, 10], name=""
):
    plt.figure(figsize=(8, 5))
    plt.hist(distances, bins=40, alpha=0.5, color="grey", label="All Distances")
    for fn in fn_distances:
        plt.axvline(
            fn,
            color="#648FFF",
            linestyle="-",
            alpha=0.2,
            label=(
                "False Negatives"
                if "False Negatives" not in plt.gca().get_legend_handles_labels()[1]
                else ""
            ),
        )
    for fp in fp_distances:
        plt.axvline(
            fp,
            color="#FFB000",
            linestyle="-",
            alpha=0.2,
            label=(
                "False Positives"
                if "False Positives" not in plt.gca().get_legend_handles_labels()[1]
                else ""
            ),
        )

    for p in percentiles:
        plt.axvline(
            np.percentile(distances, p),
            color="#DC267F",
            linestyle="-",
            linewidth=3,
            label=f"{p}th Percentile",
        )

    plt.xlabel("Distance", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.title(f"Distribution of Pairwise Distances {name}", fontsize=16)
    plt.legend()
    plt.show()


###################################
# Retrieve Pairwise Classification Threshold
##################################


def get_sample_relations_by_cutoff(
    distance_matrix: np.ndarray,
    cutoff: float,
    sample_patient_mapping: Dict[str, str],
    sample_order: List[str],
    quiet: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Returns two DataFrames:
      - belonging_df: sample pairs with distance <= cutoff
      - not_belonging_df: sample pairs with distance > cutoff
    """
    # Upper triangle indices (no self-pairs, no duplicates)
    triu_i, triu_j = np.triu_indices(len(distance_matrix), k=1)
    distances = distance_matrix[triu_i, triu_j]

    # Get sample IDs in vectorized form
    sample1_ids = np.array(sample_order)[triu_i]
    sample2_ids = np.array(sample_order)[triu_j]

    # Map to patient IDs
    patient1_ids = [sample_patient_mapping[s1] for s1 in sample1_ids]
    patient2_ids = [sample_patient_mapping[s2] for s2 in sample2_ids]

    # Build one big DataFrame
    df = pd.DataFrame(
        {
            "sample1": sample1_ids,
            "sample2": sample2_ids,
            "patient_id_1": patient1_ids,
            "patient_id_2": patient2_ids,
            "distance": distances,
        }
    )

    # Split into belonging / not belonging
    belonging_df = df[df["distance"] <= cutoff].reset_index(drop=True)
    not_belonging_df = df[df["distance"] > cutoff].reset_index(drop=True)

    if not quiet:
        print("Belonging pairs:")
        print(belonging_df)
        print("Not belonging pairs:")
        print(not_belonging_df)

    return dict(belonging=belonging_df, not_belonging=not_belonging_df)


###################################
# Evaluation Metrics
##################################


def get_evaluation_metrics(belonging, not_belonging, quiet=True):
    false_negatives = []
    false_negative_distances = []
    false_positives = []
    false_positive_distances = []
    true_negatives = []
    true_negative_distances = []
    true_positives = []
    true_positive_distances = []

    # True positives = same patient
    tp_mask = belonging["patient_id_1"] == belonging["patient_id_2"]

    true_positives = list(
        zip(belonging.loc[tp_mask, "sample1"], belonging.loc[tp_mask, "sample2"])
    )
    true_positive_distances = belonging.loc[tp_mask, "distance"].tolist()

    # False positives = different patients
    fp_mask = ~tp_mask
    false_positives = list(
        zip(belonging.loc[fp_mask, "sample1"], belonging.loc[fp_mask, "sample2"])
    )
    false_positive_distances = belonging.loc[fp_mask, "distance"].tolist()

    true_belonging = len(true_positives)
    false_belonging = len(false_positives)

    # True negatives = different patients
    tn_mask = not_belonging["patient_id_1"] != not_belonging["patient_id_2"]

    true_negatives = list(
        zip(
            not_belonging.loc[tn_mask, "sample1"], not_belonging.loc[tn_mask, "sample2"]
        )
    )
    true_negative_distances = not_belonging.loc[tn_mask, "distance"].tolist()

    # False negatives = same patients
    fn_mask = ~tn_mask
    false_negatives = list(
        zip(
            not_belonging.loc[fn_mask, "sample1"], not_belonging.loc[fn_mask, "sample2"]
        )
    )
    false_negative_distances = not_belonging.loc[fn_mask, "distance"].tolist()

    true_not_be = len(true_negatives)
    false_not_be = len(false_negatives)

    n_all = true_belonging + true_not_be + false_belonging + false_not_be

    accuracy = (true_belonging + true_not_be) / n_all if n_all > 0 else 0
    precision = (
        true_belonging / (true_belonging + false_belonging)
        if (true_belonging + false_belonging) > 0
        else 0
    )
    sensitivity = (
        true_belonging / (true_belonging + false_not_be)
        if (true_belonging + false_not_be) > 0
        else 0
    )
    balanced_accuracy = (
        0.5
        * (
            (true_belonging / (true_belonging + false_not_be))
            + (true_not_be / (true_not_be + false_belonging))
        )
        if (((true_belonging + false_not_be) > 0) & (true_belonging + false_not_be))
        else 0
    )
    f1 = (
        (2 * precision * sensitivity) / (precision + sensitivity)
        if (precision + sensitivity) > 0
        else 0
    )
    if not quiet:
        print("FP + FN:", false_belonging + false_not_be)
        print("TP:", true_belonging, "FP:", false_belonging)
        print("FN:", false_not_be, "TN:", true_not_be)
        print("Accuracy:", accuracy)
        print("balanced Accuracy:", balanced_accuracy)
        print("Precision:", precision)
        print("Sensitivity (Recall):", sensitivity)
        print("F1 Score:", f1)

    return dict(
        TP=true_belonging,
        FP=false_belonging,
        FN=false_not_be,
        TN=true_not_be,
        Accuracy=accuracy,
        Balanced_Accuracy=balanced_accuracy,
        Precision=precision,
        Sensitivity=sensitivity,
        F1=f1,
        False_Negative_Pairs=false_negatives,
        False_Negative_Distances=false_negative_distances,
        False_Positive_Pairs=false_positives,
        False_Positive_Distances=false_positive_distances,
        True_Negative_Pairs=true_negatives,
        True_Negative_Distances=true_negative_distances,
        True_Positive_Pairs=true_positives,
        True_Positive_Distances=true_positive_distances,
    )


###################################
# Optimization
##################################


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
    fp, fn, precision, f1, sensitivity = (
        result["FP"],
        result["FN"],
        result["Precision"],
        result["Sensitivity"],
        result.get("F1", 0),
    )

    if strategy == "fp+fn":
        curr_sum = fp + fn
        return curr_sum < thresholds["fp+fn"] or (
            curr_sum == thresholds["fp+fn"] and precision > thresholds["f1"]
        )

    if strategy == "fp":
        return fp < thresholds["fp"] or (
            fp == thresholds["fp"] and (fp + fn) < thresholds["fp+fn"]
        )

    if strategy == "fn":
        return fn < thresholds["fn"] or (
            fn == thresholds["fn"] and (fp + fn) < thresholds["fp+fn"]
        )

    if strategy == "F1":
        return f1 > thresholds["f1"]

    if strategy == "sensitivity":
        return sensitivity > thresholds["sensitivity"]

    if strategy == "precision":
        return precision > thresholds["precision"]

    return False
