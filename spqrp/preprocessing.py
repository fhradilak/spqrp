from sklearn.preprocessing import LabelEncoder
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm


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


def plate_correct_residuals_by_protein(
    group_data,
    individual="Patient_ID",
    sample="Sample_ID",
    figsize=(20, 12),
    impute=False,
    verbose=False,
):
    """
    Corrects for plate effects by regressing Intensity on plate and replacing intensities with residuals.
    Safely handles missing plate columns.
    """

    # --- Handle missing Intensity values ---
    if impute:
        group_data["Intensity"] = group_data.groupby(["Protein", individual])[
            "Intensity"
        ].transform(lambda x: x.fillna(x.median()))
    else:
        group_data = group_data.dropna(subset=["Intensity"])

    # --- Identify plate column (if any) ---
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

    # --- Store original Intensity for visualization ---
    original_intensity = group_data["Intensity"].copy()

    # --- Fit OLS model: Intensity ~ plate ---
    X = sm.add_constant(group_data[plate_col])
    y = group_data["Intensity"]

    try:
        model = sm.OLS(y, X).fit()
    except Exception as e:
        print(f"⚠️ OLS model fitting failed: {e}")
        return group_data

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
