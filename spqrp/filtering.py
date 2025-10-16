import pandas as pd
from typing import Optional
import logging
from sklearn.ensemble import IsolationForest
from plotly.graph_objects import Figure
import plotly.express as px

SAMPLE = "Sample_ID"
PATIENT = "Patient_ID"

##################################
# Filtering
##################################


def remove_proteins_from_filtering(df, sample="Sample_ID"):
    pio.renderers.default = "browser"
    forest_dict = by_isolation_forest(df, impute_median=True)
    by_isolation_forest_plot(forest_dict["anomaly_df"])
    df_filtered = df[~df[sample].isin(forest_dict["outlier_list"])]
    return df_filtered


default_intensity_column = "Intensity"


def long_to_wide(intensity_df: pd.DataFrame, value_name: Optional[str] = None):
    """
    This function transforms the dataframe to a wide format that
    can be more easily handled by packages such as sklearn.
    Each sample gets one row with all observations as columns.

    :param intensity_df: the dataframe that should be transformed into
        long format
        :type intensity_df: pd.DataFrame

    :return: returns dataframe in wide format suitable for use by
        packages such as sklearn
    :rtype: pd.DataFrame
    """
    values_name = default_intensity_column
    return pd.pivot(intensity_df, index=SAMPLE, columns="Protein", values=values_name)


def by_isolation_forest(
    protein_df: pd.DataFrame,
    peptide_df: pd.DataFrame = None,
    n_estimators: int = 100,
    n_jobs: int = -1,
    impute_zero=False,
    impute_median=False,
) -> dict:
    """
    This function filters out outliers using a clustering
    isolation forest approach.

    :param protein_df: a dataframe in typical protzilla long format
        on which the outlier detection is performed
    :type protein_df: pandas DataFrame
    :param n_estimators: the number of estimators used by the algorithm,
        default: 100
    :type n_estimators: integer
    :param n_jobs: Number kernels used by algorithm, default:
        all kernels (-1)
    :type n_jobs: integer

    :return: returns a Dataframe containing all samples that are not outliers and a
        dict with list of outlier sample names
    :rtype: Tuple[pandas DataFrame, dict]
    """
    try:
        transformed_df = long_to_wide(protein_df)
        if impute_zero:
            transformed_df.fillna(0, inplace=True)
        elif impute_median:
            transformed_df.fillna(transformed_df.median(skipna=True), inplace=True)

        clf = IsolationForest(
            random_state=0,
            max_samples=(len(transformed_df) // 2),
            n_jobs=n_jobs,
            n_estimators=n_estimators,
        )

        df_isolation_forest_data = pd.DataFrame(index=transformed_df.index)
        df_isolation_forest_data["IF Outlier"] = clf.fit_predict(
            transformed_df.loc[:, transformed_df.columns != SAMPLE]
        )
        df_isolation_forest_data["Anomaly Score"] = clf.decision_function(
            transformed_df
        )
        df_isolation_forest_data["Outlier"] = (
            df_isolation_forest_data["IF Outlier"] == -1
        )
        outlier_list = df_isolation_forest_data[
            df_isolation_forest_data["Outlier"]
        ].index.tolist()

        protein_df = protein_df[~(protein_df[SAMPLE].isin(outlier_list))]
        peptide_df = (
            None
            if peptide_df is None
            else peptide_df[~(peptide_df[SAMPLE].isin(outlier_list))]
        )

        return dict(
            protein_df=protein_df,
            peptide_df=peptide_df,
            outlier_list=outlier_list,
            anomaly_df=df_isolation_forest_data[["Anomaly Score", "Outlier"]],
        )
    except ValueError as e:
        msg = "Outlier Detection by IsolationForest does not accept missing values \
            encoded as NaN. Consider preprocessing your data to remove NaN values."
        return dict(
            protein_df=protein_df,
            peptide_df=peptide_df,
            outlier_list=None,
            anomaly_df=None,
            messages=[dict(level=logging.ERROR, msg=msg, trace=str(e))],
        )


PLOT_COLOR_SEQUENCE = [
    "#4A536A",
    "#CE5A5A",
    "#87A8B9",
    "#8E3325",
    "#E2A46D",
]
"""List of colors to use in plots."""

PLOT_PRIMARY_COLOR = PLOT_COLOR_SEQUENCE[0]
"""First color in list. Conventionally used for visualizing outliers."""

import plotly.io as pio

pio.renderers.default = "browser"
PLOT_SECONDARY_COLOR = PLOT_COLOR_SEQUENCE[1]
"""Second color in list."""


def by_isolation_forest_plot(output_anomaly_df, title=""):
    fig = create_anomaly_score_bar_plot(output_anomaly_df, title=title)
    fig.show()
    return fig


def create_anomaly_score_bar_plot(
    anomaly_df: pd.DataFrame,
    colour_outlier: str = PLOT_SECONDARY_COLOR,
    colour_non_outlier: str = PLOT_PRIMARY_COLOR,
    title: str = "",
) -> Figure:
    """
    This function creates a graph visualising the outlier
    and non-outlier samples using the anomaly score.

    :param anomaly_df: pandas Dataframe that contains the anomaly score for each\
    sample, including outliers and on-outliers samples
    :param colour_outlier: hex code for colour depicting the outliers.
    Default: PROTZILLA_SECONDARY_COLOR
    :param colour_non_outlier: hex code for colour depicting the
    non-outliers. Default: PROTZILLA_PRIMARY_COLOR
    :return: returns a plotly Figure object
    """

    fig = px.bar(
        anomaly_df,
        title=title,
        x=anomaly_df.index,
        y="Anomaly Score",
        hover_name=anomaly_df.index,
        hover_data={
            "Anomaly Score": True,
            "Outlier": True,
        },
        color="Outlier",
        color_discrete_map={
            False: colour_non_outlier,
            True: colour_outlier,
        },
        labels={
            "Sample_ID": "Sample_ID ",
            "Anomaly Score": "Anomaly Score ",
            "Outlier": "Outlier ",
        },
    )
    fig.update_coloraxes(showscale=False)
    fig.update_xaxes(
        categoryorder="category ascending", visible=False, showticklabels=False
    )
    fig.update_yaxes(visible=True, showticklabels=True)
    return fig
