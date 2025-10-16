from sklearn.preprocessing import LabelEncoder
import plotly.io as pio
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from imblearn.ensemble import BalancedRandomForestClassifier
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    roc_auc_score,
    roc_curve,
    precision_recall_curve,
    auc,
)
from sklearn.metrics.pairwise import nan_euclidean_distances
import statsmodels.api as sm
from spqrp.spqrp.filtering import by_isolation_forest, by_isolation_forest_plot


def pivot_df(
    df, index=["Patient_ID", "Sample_ID"], columns="Protein", values="Intensity"
):
    df_pivot = df.pivot(index=index, columns=columns, values=values).reset_index()

    df_pivot.columns.name = None  # Remove index name from columns
    df_pivot = df_pivot.rename_axis(None, axis=1)
    return df_pivot


pio.renderers.default = "browser"


def get_and_plot_outliers(train_df, test_df, title="", out_dir="."):

    forest_dict_train = by_isolation_forest(train_df, impute_median=True)
    outlier_train_df = forest_dict_train["outlier_list"]
    forest_dict_test = by_isolation_forest(test_df, impute_median=True)
    outlier_test_df = forest_dict_test["outlier_list"]

    ##Save
    out_dir = os.path.join(out_dir, "outlier_sample_ids")

    # Create folder if it does not exist
    os.makedirs(out_dir, exist_ok=True)
    with open(
        os.path.join(out_dir, "train_df_isolation_forest_median_imputed.txt"), "w"
    ) as fp:
        for item in outlier_train_df:
            fp.write(f"{item}\n")
        print("Train file written.")
    with open(
        os.path.join(out_dir, "test_df_isolation_forest_median_imputed.txt"), "w"
    ) as fp:
        for item in outlier_test_df:
            fp.write(f"{item}\n")
        print("Test file written.")
    ##Plot
    by_isolation_forest_plot(forest_dict_train["anomaly_df"], title=f"train_df_{title}")
    by_isolation_forest_plot(forest_dict_test["anomaly_df"], title=f"test_df_{title}")

    print("Number Outliers Train: ", len(outlier_train_df))
    print("Number Outliers Test: ", len(outlier_test_df))

    return outlier_train_df, outlier_test_df


def compute_feature_differences(sample1, sample2):
    """Compute the feature-wise differences between two samples."""
    return sample1 - sample2


def compute_nan_euclidean_distance(sample1, sample2):
    """Compute the Euclidean distance between two samples."""
    sample1 = np.array(sample1).reshape(1, -1)  # Convert to 2D
    sample2 = np.array(sample2).reshape(1, -1)  # Convert to 2D
    return nan_euclidean_distances(sample1, sample2)[0, 0]


def create_pairs(X_data, y_data, seen_pairs, pairwise, labels, pair_indices):
    for i in range(len(X_data)):
        for j in range(i + 1, len(X_data)):  # Ensure i < j to avoid duplicate pairs
            pair = tuple(
                sorted([i, j])
            )  # Sorting ensures that (i, j) and (j, i) are considered the same
            if pair not in seen_pairs:
                pairwise.append((X_data.iloc[i].values, X_data.iloc[j].values))
                labels.append(1 if y_data.iloc[i] == y_data.iloc[j] else 0)
                pair_indices.append((i, j))  # Store the indices of the pair
                seen_pairs.add(pair)


###############################
#
# Missclassified Samples
#
###############################


def get_and_print_missclassified_samples(
    df_pivot,
    misclassified_indices,
    pair_indices,
    seen_pairs,
    misclassified_samples,
    y_pairs,
    y_pred_adjusted,
    results_dict,
):
    for idx in misclassified_indices:
        # Retrieve the pair indices from pair_indices_test
        pair_indices_new = pair_indices[
            idx
        ]  # This gives us a tuple (i, j) with the original indices
        # Get the original samples using the pair indices
        sample_1_idx, sample_2_idx = pair_indices_new
        # Get the actual sample names using these indices from df_pivot
        name_1 = df_pivot.iloc[sample_1_idx]["Sample_ID"]
        name_2 = df_pivot.iloc[sample_2_idx]["Sample_ID"]

        # Create a tuple of the sample names
        sample_pair = tuple(
            sorted([name_1, name_2])
        )  # Sorting to ensure consistent order
        # If we haven't already seen this pair, add it to the list of misclassified pairs
        if sample_pair not in seen_pairs:
            seen_pairs.add(sample_pair)
            misclassified_samples.append(
                {
                    "Sample 1": name_1,  # Sample 1's name
                    "Sample 2": name_2,  # Sample 2's name
                    "True Label": y_pairs[idx],  # True label for the pair
                    "Predicted Label": y_pred_adjusted[
                        idx
                    ],  # Predicted label for the pair
                }
            )
    results_dict["misclassified_pairs"] = misclassified_samples

    # Print the names of the misclassified pairs
    print("Misclassified Pairs:")
    for pair in misclassified_samples:
        print(
            f"Sample 1: {pair['Sample 1']}, Sample 2: {pair['Sample 2']}, True:{pair['True Label']}, Predicted Label: {pair['Predicted Label']} "
        )
        print("********")
    return results_dict


###############################
#
# Model Helpers
#
###############################


def print_and_add_feature_importance(
    importances, feature_names, results_dict, test=False
):
    sorted_features = sorted(zip(importances, feature_names), reverse=True)
    key_importances = "feature_importances"
    if test:
        key_importances = "feature_importances_test"
    n = 10

    # Display top 3 important features (both original and engineered)
    print(f"Top {n} Important Features:")
    for importance, name in sorted_features[:n]:
        print(f"{name}: {importance:.4f}")
    for importance, name in sorted_features:  # Top 5 features (adjust if needed)
        results_dict[key_importances].append(
            {"feature": name, "importance": importance}
        )


def plot_evaluation(y_test, y_probs, y_prob):
    # ROC Curve
    fpr, tpr, thresholds = roc_curve(y_test, y_probs)
    roc_auc = auc(fpr, tpr)

    # Precision-Recall Curve
    precision, recall, pr_thresholds = precision_recall_curve(y_test, y_probs)

    # Plot ROC curve
    plt.figure(figsize=(10, 6))
    plt.plot(fpr, tpr, color="b", label=f"ROC curve (AUC = {roc_auc:.2f})")
    plt.plot([0, 1], [0, 1], color="gray", linestyle="--")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("Receiver Operating Characteristic (ROC) Curve")
    plt.legend(loc="lower right")
    plt.show()

    # Plot Precision-Recall curve
    plt.figure(figsize=(10, 6))
    plt.plot(recall, precision, color="b", label="Precision-Recall curve")
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title("Precision-Recall Curve")
    plt.legend(loc="lower left")
    plt.show()

    # Optionally, you can plot the decision threshold vs. accuracy curve (or F1 score, etc.)
    accuracies = []
    for threshold in thresholds:
        predictions = (y_probs >= threshold).astype(
            int
        )  # Apply threshold to get predicted class labels
        accuracy = np.mean(
            predictions == y_test
        )  # Calculate accuracy for that threshold
        accuracies.append(accuracy)

    # Plot accuracy vs threshold
    plt.figure(figsize=(10, 6))
    plt.plot(thresholds, accuracies, color="b", label="Accuracy vs Threshold")
    plt.xlabel("Decision Threshold")
    plt.ylabel("Accuracy")
    plt.title("Threshold vs Accuracy")
    plt.legend(loc="best")
    plt.show()

    _, threshold = get_threshold(y_test, y_prob)
    _, roc_threshold = get_threshold(y_test, y_prob, method="ROC")
    _, J_threshold = get_threshold(y_test, y_prob, method="J")
    _, Min_FP_threshold = get_threshold(y_test, y_prob, method="MinFP")

    # Plot probability distributions for each class
    plt.figure(figsize=(10, 6))
    sns.histplot(
        y_prob[y_test == 0],
        color="blue",
        bins=30,
        alpha=0.6,
        label="Class 0 (Not Belonging)",
    )
    sns.histplot(
        y_prob[y_test == 1],
        color="red",
        bins=30,
        alpha=0.6,
        label="Class 1 (Belonging)",
    )

    # Plot threshold line
    plt.axvline(
        x=threshold, color="black", linestyle="--", label=f"Threshold = {threshold:.2f}"
    )
    plt.axvline(
        x=roc_threshold,
        color="green",
        linestyle="--",
        label=f"Threshold = {roc_threshold:.2f}",
    )
    plt.axvline(
        x=J_threshold,
        color="blue",
        linestyle="--",
        label=f"Threshold = {J_threshold:.2f}",
    )
    plt.axvline(
        x=Min_FP_threshold,
        color="yellow",
        linestyle="--",
        label=f"Min_FP_Threshold = {Min_FP_threshold:.2f}",
    )

    # Labels and title
    plt.xlabel("Predicted Probability")
    plt.ylabel("Count")
    plt.title("Decision Threshold for Classification")
    plt.legend()
    plt.show()


def evaluate_model(y_test_pairs, y_pred_adjusted, y_prob, sample_decision_curve=True):
    auc_score = roc_auc_score(y_test_pairs, y_pred_adjusted)
    print(f"AUC-ROC Score: {auc_score}")
    print(f"Accuracy: {accuracy_score(y_test_pairs, y_pred_adjusted)}")
    print(
        f"Classification Report:\n{classification_report(y_test_pairs, y_pred_adjusted)}"
    )
    print("*************")
    false_negatives = np.where((y_pred_adjusted == 0) & (y_test_pairs == 1))[0]
    false_positives = np.where((y_pred_adjusted == 1) & (y_test_pairs == 0))[0]
    true_negatives = np.where((y_pred_adjusted == 0) & (y_test_pairs == 0))[0]
    true_positives = np.where((y_pred_adjusted == 1) & (y_test_pairs == 1))[0]
    print(f"False Negatives (FN) indices: {false_negatives}")
    print(f"False Positives (FP) indices: {false_positives}")

    print(f"Number: FP:{len(false_positives)}, TP:{len(true_positives)}")
    print(f"Number: FN:{len(false_negatives)}, TN:{len(true_negatives)}")
    if sample_decision_curve:
        plot_evaluation(y_test_pairs, y_pred_adjusted, y_prob)


def get_threshold(
    y_test_pairs, y_prob, method="ROC", labels_new=[], y_prob_new=[], max_fpr=0.01
):

    if method == "ROC":
        fpr, tpr, thresholds = roc_curve(y_test_pairs, y_prob)
        optimal_idx = np.argmax(tpr - fpr)
        optimal_threshold = thresholds[optimal_idx]
        print(optimal_threshold)

    elif method == "F1":
        precision, recall, thresholds_f1 = precision_recall_curve(
            labels_new, y_prob_new
        )
        f1_scores = 2 * (precision * recall) / (precision + recall)
        optimal_threshold = thresholds_f1[np.argmax(f1_scores)]
        print("Optimal threshold based on F1 score:", optimal_threshold)
    elif method == "J":
        fpr, tpr, thresholds = roc_curve(y_test_pairs, y_prob)
        # Calculate Youden's J statistic
        youden_j = tpr - fpr
        optimal_idx = np.argmax(youden_j)
        optimal_threshold = thresholds[optimal_idx]
        print(f"Optimal threshold (Youden’s J): {optimal_threshold:.3f}")
    elif method == "MinFP":
        fpr, tpr, thresholds = roc_curve(y_test_pairs, y_prob)
        valid_thresholds = thresholds[fpr <= max_fpr]

        if len(valid_thresholds) > 0:
            optimal_threshold = valid_thresholds[
                -1
            ]  # Take the highest threshold under max_fpr
        else:
            optimal_threshold = thresholds[
                0
            ]  # If all thresholds have high FPR, take the lowest

        print(f"Optimal threshold (FPR ≤ {max_fpr}): {optimal_threshold:.3f}")

    y_pred_adjusted = (y_prob > optimal_threshold).astype(int)
    return y_pred_adjusted, optimal_threshold


def plot_pairwise_probabilities_for_sample(
    sample_index, y_prob, pair_indices_test, df_pivot, optimal_threshold
):
    """
    Plot pairwise classification probabilities for a single sample against all others using precomputed y_prob.

    Parameters:
    - sample_index: index (row) in df_pivot for the target sample
    - y_prob: array of predicted probabilities for test pairs
    - pair_indices_test: list of index pairs (tuples) for each test pair
    - df_pivot: original DataFrame containing sample names (must include "Sample_ID" column)
    - compute_euclid: currently unused, included for compatibility
    """

    sample_name = df_pivot.iloc[sample_index]["Sample_ID"]

    pairwise_probs = []
    compared_sample_names = []

    for i, (idx1, idx2) in enumerate(pair_indices_test):
        if sample_index in [idx1, idx2]:
            # Identify the *other* sample in the pair
            other_idx = idx2 if sample_index == idx1 else idx1
            other_sample_name = df_pivot.iloc[other_idx]["Sample_ID"]

            pairwise_probs.append(y_prob[i])
            compared_sample_names.append(other_sample_name)

    # Plot
    plt.figure(figsize=(12, 6))
    plt.scatter(range(len(pairwise_probs)), pairwise_probs, c="blue", alpha=0.7)
    plt.axhline(
        y=optimal_threshold,
        color="red",
        linestyle="--",
        label=f"Optimal Threshold ({optimal_threshold})",
    )
    plt.title(f"Pairwise Probabilities for Sample: {sample_name}")
    plt.xlabel("Compared Sample Index")
    plt.ylabel("Predicted Probability (Class 1)")
    plt.xticks(
        ticks=range(len(compared_sample_names)),
        labels=compared_sample_names,
        rotation=90,
    )
    plt.grid(True)
    plt.tight_layout()
    plt.show()


###############################
#
# Training
#
###############################


def train_pairwise_balanced_rand_forest(
    X_train,
    y_train,
    X_test,
    y_test,
    df_pivot_test,
    k=0,
    compute_euclid=True,
    method="F1",
    importance_method="impurity",
    plots_per_sample=False,
    sample_decision_curve=False,
    absolute=False,
):
    pairwise_train = []
    pairwise_test = []
    labels_train = []
    labels_test = []

    results_dict = {
        "feature_importances": [],
        "feature_importances_test": [],
        "misclassified_pairs": [],
    }
    # print(y_train)
    # print(y_test)
    # Generate pairwise training samples
    seen_train_pairs = set()
    pair_indices_train = []
    create_pairs(
        X_data=X_train,
        y_data=y_train,
        seen_pairs=seen_train_pairs,
        pairwise=pairwise_train,
        labels=labels_train,
        pair_indices=pair_indices_train,
    )

    # Set to track seen pairs for testing
    seen_test_pairs = set()
    pair_indices_test = []

    create_pairs(
        X_data=X_test,
        y_data=y_test,
        seen_pairs=seen_test_pairs,
        pairwise=pairwise_test,
        labels=labels_test,
        pair_indices=pair_indices_test,
    )
    # print("Pairwise_train",pairwise_train)
    original_feature_names = X_train.columns.tolist()
    feature_names = []
    for feature in original_feature_names:
        feature_names.append(f"diff_{feature}")

    if compute_euclid:
        feature_names.append("euclidean_distance")

    print("Train_indices_test", pair_indices_train)
    print("Pair_indices_test", pair_indices_test)
    # Convert pairwise data into feature engineering (difference + distance)
    X_train_pairs = []
    X_test_pairs = []

    # For each pair, create feature differences and Euclidean distance
    for pair in pairwise_train:
        feature_diff = compute_feature_differences(pair[0], pair[1])
        if compute_euclid:
            euclid_dist = compute_nan_euclidean_distance(pair[0], pair[1])
            X_train_pairs.append(np.concatenate([feature_diff, [euclid_dist]]))
        else:
            X_train_pairs.append(feature_diff)

    for pair in pairwise_test:
        feature_diff = compute_feature_differences(pair[0], pair[1])
        if compute_euclid:
            euclid_dist = compute_nan_euclidean_distance(pair[0], pair[1])
            X_test_pairs.append(np.concatenate([feature_diff, [euclid_dist]]))
        else:
            X_test_pairs.append(feature_diff)

    # Convert list to DataFrame , to numpy arrays
    X_train_pairs = pd.DataFrame(X_train_pairs, columns=feature_names)
    X_test_pairs = pd.DataFrame(X_test_pairs, columns=feature_names)
    y_train_pairs = np.array(labels_train)
    y_test_pairs = np.array(labels_test)

    ######### Train model on the pairwise features
    clf = BalancedRandomForestClassifier(n_estimators=100, random_state=42)
    clf.fit(X_train_pairs, y_train_pairs)
    #########Test
    y_prob = clf.predict_proba(X_test_pairs)[:, 1]  # Get probability for class `1`

    ########Get optimal Threshold  ROC, F1
    y_pred_adjusted, optimal_treshold = get_threshold(
        y_test_pairs, y_prob, labels_new=labels_test, y_prob_new=y_prob, method=method
    )

    ########Evaluate the model
    k = k + 1
    print(f"Fold {k}")
    evaluate_model(y_test_pairs, y_pred_adjusted, y_prob, sample_decision_curve)
    if absolute:
        X_train_pairs = X_train_pairs.abs()
        X_test_pairs = X_test_pairs.abs()
    if plots_per_sample:
        for i in range(0, len(df_pivot_test)):
            plot_pairwise_probabilities_for_sample(
                sample_index=i,
                y_prob=y_prob,
                pair_indices_test=pair_indices_test,
                df_pivot=df_pivot_test,
                optimal_threshold=optimal_treshold,
            )

    #######Missclassified
    # Identify the misclassified pairwise instances
    misclassified_indices = np.where(y_pred_adjusted != y_test_pairs)[0]
    misclassified_samples = []
    seen_pairs = set()
    results_dict = get_and_print_missclassified_samples(
        df_pivot_test,
        misclassified_indices,
        pair_indices_test,
        seen_pairs,
        misclassified_samples,
        y_test_pairs,
        y_pred_adjusted,
        results_dict,
    )

    #######Importance
    if importance_method == "impurity":
        importances = clf.feature_importances_
        feature_names = X_train_pairs.columns
        print_and_add_feature_importance(importances, feature_names, results_dict)
    train_pairs_names = [
        (X_train.index[i], X_train.index[j]) for i, j in pair_indices_train
    ]
    test_pairs_names = [
        (X_test.index[i], X_test.index[j]) for i, j in pair_indices_test
    ]
    return {
        "results_dict": results_dict,
        "train_pairs_names": train_pairs_names,
        "test_pairs_names": test_pairs_names,
        "clf": clf,
        "X_train_pairs": X_train_pairs,
    }


def filter_by_occurrence(df, cutoff=0.7):
    """
    Keep proteins that are present in at least `cutoff`% of samples.

    Parameters
    ----------
    df : pd.DataFrame
        Long-format dataframe with columns ["Protein", "Sample_ID", "Intensity"].
    cutoff : int
        Percentage of samples in which a protein must appear (0–1).

    Returns
    -------
    pd.DataFrame
        Filtered dataframe with only proteins passing the cutoff.
    """
    if not (0 <= cutoff <= 1):
        raise ValueError("Cutoff must be between 0 and 1.")
    n_samples = df["Sample_ID"].nunique()
    counts = df.groupby("Protein")["Intensity"].apply(lambda x: x.notna().sum())
    fractions = counts / n_samples
    proteins_to_keep = fractions[fractions >= cutoff].index

    return df[df["Protein"].isin(proteins_to_keep)]


def log_transform(df):
    df["Intensity"] = np.log2(df["Intensity"])
    return df


def revert_log_transform(df):
    df["Intensity"] = 2 ** df["Intensity"]
    return df


def plate_correct_residuals_by_protein(
    group_data,
    individual="Patient_ID",
    sample="Sample_ID",
    figsize=(20, 12),
    impute=False,
    verbose=False,
):
    # Drop NA values from 'Intensity'
    if impute:
        group_data["Intensity"] = group_data.groupby(["Protein", individual])[
            "Intensity"
        ].transform(lambda x: x.fillna(x.median()))
    else:
        group_data = group_data.dropna(subset=["Intensity"])

    plate_col = next(
        (col for col in ["plate", "Plate"] if col in group_data.columns), None
    )

    if plate_col is None:
        print("⚠️ No 'plate' or 'Plate' column found. Skipping plate effect correction.")
        return group_data  # Return unchanged

    # --- Encode plate if categorical ---
    label_encoder = LabelEncoder()
    try:
        group_data[plate_col] = label_encoder.fit_transform(
            group_data[plate_col].astype(str)
        )
    except Exception as e:
        print(f"⚠️ Warning: Could not encode plate column ({e}). Skipping correction.")
        return group_data

    # Store original Intensity for visualization
    original_intensity = group_data["Intensity"].copy()

    # Fit the OLS model (Intensity ~ plate)
    X = sm.add_constant(group_data[plate_col])  # Adding constant for the intercept
    y = group_data["Intensity"]
    model = sm.OLS(y, X).fit()

    print(model.summary())

    # Check if the plate effect is statistically significant
    plate_p_value = model.pvalues[plate_col]

    if plate_p_value >= 0.05:
        print(
            "\nWARNING: The plate effect is not statistically significant (p-value >= 0.05)."
        )
        print("You may want to reconsider applying plate correction.")

    # Add residuals to the dataframe
    group_data["Intensity"] = model.resid

    if verbose:
        # Create a 2x2 grid for boxplots (before and after correction for Protein and Sample)
        fig, axes = plt.subplots(2, 2, figsize=figsize)

        # Before correction: Boxplot for original intensity per Protein
        sns.boxplot(ax=axes[0, 0], data=group_data, x="Protein", y=original_intensity)
        axes[0, 0].set_title("Before Plate Effect Correction (Protein)")
        axes[0, 0].set_xlabel("Protein")
        axes[0, 0].set_ylabel("Intensity")
        axes[0, 0].tick_params(axis="x", labelsize=7, rotation=90)

        # After correction: Boxplot for intensity after residual correction per Protein
        sns.boxplot(ax=axes[0, 1], data=group_data, x="Protein", y="Intensity")
        axes[0, 1].set_title("After Plate Effect Correction (Protein)")
        axes[0, 1].set_xlabel("Protein")
        axes[0, 1].set_ylabel("Intensity")
        axes[0, 1].tick_params(axis="x", labelsize=7, rotation=90)

        # Before correction: Boxplot for original intensity per Sample
        sns.boxplot(ax=axes[1, 0], data=group_data, x=sample, y=original_intensity)
        axes[1, 0].set_title("Before Plate Effect Correction (Sample)")
        axes[1, 0].set_xlabel("Sample")
        axes[1, 0].set_ylabel("Intensity")
        axes[1, 0].tick_params(axis="x", labelsize=7, rotation=90)

        # After correction: Boxplot for intensity after residual correction per Sample
        sns.boxplot(ax=axes[1, 1], data=group_data, x=sample, y="Intensity")
        axes[1, 1].set_title("After Plate Effect Correction (Sample)")
        axes[1, 1].set_xlabel("Sample")
        axes[1, 1].set_ylabel("Intensity")
        axes[1, 1].tick_params(axis="x", labelsize=7, rotation=90)

        plt.tight_layout()
        plt.show()

        #################################################Checking
        mean_difference = (group_data["Intensity"] - original_intensity).mean()
        print(
            f"Mean difference between original and corrected intensity: {mean_difference}"
        )

        # Store residuals separately for inspection
        residuals = model.resid
        group_data["Residuals"] = residuals

        sns.countplot(x=plate_col, data=group_data)
        plt.title("Distribution of Plates")
        plt.show()

        plt.figure(figsize=(20, 5))
        sns.boxplot(data=group_data, x="Protein", y="Intensity")
        sns.boxplot(
            data=group_data, x="Protein", y=original_intensity, color="red"
        )  # Original in red
        plt.title("Original vs Corrected Intensities")
        plt.xlabel("Protein")
        plt.ylabel("Intensity")
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.show()

    return group_data


def normalize_medianintensity(
    dataset, string_of_pool="", revert_log=False, sample="Sample_ID"
):

    if revert_log:
        revert_log_transform(dataset)

    # Calculate per-sample median intensity
    sample_median = dataset.groupby(sample)["Intensity"].median()
    median_of_summed_intensities = sample_median.median()

    dataset_before = dataset.copy()

    # Filter out specified pool (if provided)
    if string_of_pool:
        dataset = dataset[~dataset[sample].str.contains(string_of_pool)]

    # Normalize intensity
    if revert_log:
        dataset = (
            dataset.groupby(sample)
            .apply(
                lambda x: x.assign(
                    NormalizeFactor=median_of_summed_intensities
                    / x["Intensity"].median(),
                    Intensity=x["Intensity"]
                    * (median_of_summed_intensities / x["Intensity"].median()),
                )
            )
            .reset_index(drop=True)
        )
    else:
        # division in raw space is equivalent to subtraction in log space, we normalize using
        dataset = dataset.groupby(sample).apply(
            lambda x: x.assign(
                NormalizeFactor=x["Intensity"].median() - dataset["Intensity"].median(),
                Intensity=x["Intensity"]
                - (x["Intensity"].median() - dataset["Intensity"].median()),
            )
        )

    # Create boxplots before and after normalization
    plt.figure(figsize=(25, 5))

    # Boxplot before normalization
    plt.subplot(1, 2, 1)
    sns.boxplot(data=dataset_before, x="Sample_ID", y="Intensity")
    plt.title("Before Normalization")
    plt.xlabel("Sample_ID")
    plt.ylabel("log2(Intensity)")
    plt.xticks(rotation=90, fontsize=8)

    # Boxplot after normalization
    plt.subplot(1, 2, 2)
    sns.boxplot(data=dataset, x="Sample_ID", y="Intensity")
    plt.title("After Normalization")
    plt.xlabel("Sample_ID")
    plt.ylabel("log2(Intensity)")
    plt.xticks(rotation=90, fontsize=8)

    plt.tight_layout()
    plt.show()

    return dataset


def train_with_normalise(
    df,
    threshold=0.7,
    test_size=0.3,
    plate_corrected=True,
    individual="Patient_ID",
    sample="Sample_ID",
    compute_euclid=False,
    method="F1",
    outlier_removal=True,
    train_individuals=None,
    test_individuals=None,
    sample_decision_curve=False,
    importance_method="impurity",
    plot_per_sample=False,
    absolute=False,
):
    # Reindex to full grid and fill missing Intensity with NaN
    # This catches missing samples in protein groups influence on occurence filtering
    proteins = df["Protein"].unique()
    samples = df[sample].unique()
    full_index = pd.MultiIndex.from_product(
        [proteins, samples], names=["Protein", sample]
    )
    df_complete = df.set_index(["Protein", sample]).reindex(full_index).reset_index()
    # filter by threshold
    df_clean = filter_by_occurrence(df_complete, threshold)

    log_transform(df_clean)

    if train_individuals is None and test_individuals is None:
        # Get unique Individuals for Grouping
        unique_individuals = df_clean[individual].unique()
        train_individuals, test_individuals = train_test_split(
            unique_individuals, test_size=test_size, random_state=42
        )
        print("train_individuals", train_individuals)
        print("test_individuals", test_individuals)

    df_train = df_clean[df_clean[individual].isin(train_individuals)]
    df_test = df_clean[df_clean[individual].isin(test_individuals)]

    train_proteins = set(df_train["Protein"].unique())
    test_proteins = set(df_test["Protein"].unique())

    extra_in_test = test_proteins - train_proteins
    print("Proteins only in test set:", extra_in_test)

    df_train_norm = normalize_medianintensity(df_train)
    df_test_norm = normalize_medianintensity(df_test)
    df_test_normalized = df_test_norm
    df_train_normalized = df_train_norm

    if plate_corrected:
        df_test_normalized = plate_correct_residuals_by_protein(
            df_test_normalized, individual=individual
        )
        df_train_normalized = plate_correct_residuals_by_protein(
            df_train_normalized, individual=individual
        )

    df_train_pivot = pivot_df(df_train_normalized, index=[individual, sample])
    df_test_pivot = pivot_df(df_test_normalized, index=[individual, sample])
    ###Ensure both sets have same features + order
    df_test_pivot = df_test_pivot.reindex(columns=df_train_pivot.columns)

    ##### Outlier removal
    if outlier_removal:
        ##### Outlier identification
        train_outliers, test_outliers = get_and_plot_outliers(
            df_train_normalized, df_test_normalized
        )
        df_train_pivot_clean = df_train_pivot[
            ~df_train_pivot[sample].isin(train_outliers)
        ]
        df_test_pivot_clean = df_test_pivot[~df_test_pivot[sample].isin(test_outliers)]
    else:
        df_train_pivot_clean = df_train_pivot
        df_test_pivot_clean = df_test_pivot
    ##### Prepare X,y
    X_train = df_train_pivot_clean.drop([individual, sample], axis=1)
    y_train = df_train_pivot_clean[individual]

    X_test = df_test_pivot_clean.drop([individual, sample], axis=1)
    y_test = df_test_pivot_clean[individual]

    return train_pairwise_balanced_rand_forest(
        X_train=X_train,
        y_train=y_train,
        X_test=X_test,
        y_test=y_test,
        df_pivot_test=df_test_pivot_clean,
        compute_euclid=compute_euclid,
        method=method,
        importance_method=importance_method,
        plots_per_sample=plot_per_sample,
        sample_decision_curve=sample_decision_curve,
        absolute=absolute,
    )
