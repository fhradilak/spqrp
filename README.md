[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
# SPQRP 
Sample Provenance Quality Resolver in Proteomics is a tool that provides quality assessment for plasma-MS-proteome studies. Recent advancements in MS technology and lab methods opened the door for large-scale proteomics but also led to a growing concern regarding sample mix-ups. We built this package to help scientists evaluate whether their sample data is safe for further analysis.

The package now offers functions that can also be used in R code.
Distance Calculation & Threshold-based approach:
- ```result <- spqrp$perform_distance_evaluation_on_ranked_proteins(df = df)``` can be used to apply a distance metric on your data set and retrieve performance metrics based on a threshold.
- ```results <- spqrp$optimize_parameters(df = df)``` can be used to find the optimized parameters for perform_distance_evaluation_on_ranked_proteins for a range of n  proteins used.
Clustering-based approach & visualization:
- `cluster_samples_iteratively`: creates a graph in a dimensionality reduced space (e.g. with UMAP) based on the distances from perform_distance_evaluation_on_ranked_proteins and can cluster iteratively based on the n nearest neighbours per sample.
-  `plot_distances_neighbours_with_coloring_hue`: visualizes the graph from cluster_samples_iteratively based on the original data labeling.


## 🗃️ Input Data Format
1. Cohort data: protein intensities per sample
2. Protein importance ranking list (optional - can use default list)

#### Cohort data
The following columns have to be present in the input df.

| Column Name | Description                        |
|-------------|------------------------------------|
| `Sample_ID` | Unique identifier for the sample   |
| `Patient_ID`| Identifier for the patient         |
| `Protein`   | Name or identifier of the protein  |
| `Intensity` | Measured intensity (numeric value) |
#### Protein ranking
| Column Name | Description                        |
|-------------|------------------------------------|
| `Protein`   | Name or identifier of the protein  |
| `Importance`| Importance/ Ranking of the protein |



# Requirements to use SPQRP in R/ RStudio
## 1. Python Installation

> ⚠️ **Note:** We recommend installing Python from the [official Python website](https://www.python.org/downloads/).  
> Issues have been observed when using the **Microsoft Store version**, particularly with the `reticulate` package in R.  
> You may want to **uninstall the Windows Store version** before proceeding with the official installation.

- ✅ **Best source**: Install from the official Python website:  
  <img src="https://github.com/user-attachments/assets/e303f6be-96a1-4422-9add-6b97980acde9" width="500"/>

- ✅ **During installation**, **tick the "Add Python to PATH" box**:  
  <img src="https://github.com/user-attachments/assets/5ed8d756-6d5c-4380-81d5-4ab1a6100ea2" width="500"/>

- 🔍 **Find the path to your Python installation:**
  - Open a command line tool (e.g., search for `cmd` in the Windows search bar).
  - On Windows, type:

    ```
    where python
    ```

  - Example output:  
    <img src="https://github.com/user-attachments/assets/3b593c31-a0df-4518-a82c-0be6834cf844" width="300"/>

  - Use the **path before `python.exe`** for integration into R (e.g., when setting `use_python` or `RETICULATE_PYTHON`).


# How to run the Package in Rstudio
## 1. Pull the python package

```r
library(reticulate)
#use_python("C:\\Program Files\\Python312\\", required = TRUE) --> you might need to explicitly set the path to your python installation from the last step from the Python installation
py_require("git+https://github.com/fhradilak/spqrp.git")
virtualenv_create("myenv")
use_virtualenv("myenv",required=TRUE)
py_install("git+https://github.com/fhradilak/spqrp.git",envname="myenv",method="virtualenv")
spqrp <- import("spqrp")
```

( alternatively if you manage your own virtual environment
```r
library(reticulate)
py_require("git+https://github.com/fhradilak/spqrp.git")
spqrp <- import("spqrp")

```
)
## 2. Running Analysis (with Default Settings) - perform_distance_evaluation_on_ranked_proteins

- `df`: **Protein Intensity Input** *(mandatory)*  
  Input dataframe containing protein intensity values.

- `top_importance_path = "path/to/ranking"`  
  Path to precomputed protein importance rankings.

- `n = 10`  
  Number of proteins used for the distance metric.

- `p = 0.5`  
  Cutoff threshold for the distance metric.

- `remove_list = None`  
  List of proteins to exclude from the distance calculation.

- `metric = "correlation"`  
  Distance metric to use. See [available metrics](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.pairwise_distances.html). Tested options: `correlation`, `euclidean`, `fractional`.

- `fractional_p = 0.5`  
  Parameter for the fractional metric. For more info, see [this reference](https://dl.acm.org/doi/10.5555/645504.656414).


### Example with default values
```r
df <- read.csv("path\to\data.csv")

# ensure data is in the right format e.g. 
# df$Patient_ID <- df$P
# df$Sample_ID <- df$Sample

result <- spqrp$perform_distance_evaluation_on_ranked_proteins(df = df)
```
### Example with set values
```r
df <- read.csv("path\to\data.csv")

# ensure data is in the right format e.g. 
# df$Patient_ID <- df$P
# df$Sample_ID <- df$Sample

top_importance_path ="path\to\top\importance.csv"
result <- spqrp$perform_distance_evaluation_on_ranked_proteins(df_filtered = df, top_importance_path = top_importance_path, n=10, p=0.5, remove_list=["p123], metric = "fractional", fractional_p=0.5)
```

## 3. Running analysis with optimizing cutoff for the distance metric -optimize_parameters

### Default Settings

- `remove_list = []`  
  List of proteins to exclude from consideration.

- `metric = "correlation"`  
  Metric used for similarity/distance (e.g., correlation, euclidean, etc.).

- `path_classification_importance = "path/to/ranking"`  
  Path to the file containing importance scores or rankings.

- `top_median = False`  
  If `False`, uses **standard deviation ratio** for protein ranking instead of top median.

- `log_file = "optimization_log.txt"`  
  File to log the output of all optimization runs.

- `range = range(2, 50, 1)`  
  Range of `n` (number of proteins) to evaluate during optimization.

- `optimization_strategy = "fpfn"`  
  Metric to optimize for. Options:
  - `fpfn`
  - `F1`
  - `fn`
  - `fp`


### Example with default settings
```r
result <- spqrp$optimize_parameters(df = df)
```

## 4. Clustering Approach and Visualizing Results — `cluster_samples_iteratively` & `plot_distances_neighbours_with_coloring_hue`

### 📌 Function Purpose

`cluster_samples_iteratively` clusters the samples based on the precalculated distances from the result from `result <- spqrp$perform_distance_evaluation_on_ranked_proteins()`. `n_neighbours` is based on the knowledge about the datasets expected sample size. `max_component_size` is a parameter to adjust the maximum final cluster size, usally set to n_neighbours+1.
### ⚙️ Parameters
- **`result`** *(mandatory)*
  Output from `spqrp$perform_distance_evaluation_on_ranked_proteins()`
- **`method`** *(mandatory)*
  Dimensionality Reduction Method ("UMAP", "PCA", "MDS").
- **`random_state`** *(default: `42`)*
  Random state seed for reproducability for the dimensionality reduction methods.  
- **`n_neighbors`** *(default: `1`)*
  Number of samples that should have an edge drawn for each sample. The n connected samples are n the nearest neighbours.
- **`max_component_size`** *(default: `None`)*
  If None, no further restrictions are performed on the graph. If equals an Integer, the respective longest edges are iteratively removed for all connected components until only components/ clusters with maximum max_component_size remain.
- **`n_umap_neighbors`** *(default: `15`)*
  Parameter that "controls how UMAP balances local versus global structure in the data".([https://www.python.org/downloads/](https://umap-learn.readthedocs.io/en/latest/parameters.html))
- **`precomputed_graph`** *(default: `None`)*
  Option to input an already precomputed graph and to apply the iterative clustering approach.

### ⚙️ Returns
- **`G`** The networkx Graph.
- **`coords_2d`** The coordinates of each sample/node in the dimensionality reduction view.

---
The function `plot_distances_neighbours_with_coloring_hue()` visualizes sample-to-sample relationships from the clustering using either PCA or MDS dimensionality reduction. It highlights patient-specific groupings and nearest neighbors using color-coded nodes and styled edges in a plot.

### ⚙️ Parameters

- **`df`** *(mandatory)*  
  Input dataframe containing **protein intensity values**.

- **`G`** *(mandatory)*  
  NetworkX Graph object representing sample connections, typically the output of `cluster_samples_iteratively()`.

- **`coords_2d`** *(mandatory)*  
  Coordinates for each sample/node in 2D space, produced by PCA, MDS, or UMAP.

- **`method`** *(default: `'UMAP'`)*  
  String specifying which dimensionality reduction method to label axes for:
  - `'PCA'` – principal component analysis  
  - `'MDS'` – multidimensional scaling  
  - `'UMAP'` – uniform manifold approximation and projection  

- **`subset_samples`** *(default: `None`)*  
  List of specific sample IDs to visualize. Dimensionality reduction is still computed on the full dataset, but only these samples are drawn.

- **`highlight_singletons`** *(default: `False`)*  
  If `True`, marks singleton nodes (samples forming their own cluster & being the only sample for their patient ID) with a blue square.

- **`highlight_single_samples_missing_connections`** *(default: `False`)*  
  If `True`, marks samples that are disconnected but there is at least 1 other sample with the same patient ID in the input df, with a pink circle.

- **`figsize`** *(default: `(20,20)`)*  
  Tuple specifying the figure size in inches `(width, height)`.

- **`label_patient_only`** *(default: `False`)*  
  If `True`, labels nodes by patient ID only; otherwise labels nodes by sample ID.

- **`label_offset_x`** *(default: `0.01`)*  
  Horizontal offset for node labels to avoid overlapping the node marker.

- **`label_offset_y`** *(default: `0.01`)*  
  Vertical offset for node labels.

- **`label_font`** *(default: `6`)*  
  Font size for node labels.

- **`df_name`** *(default: `'DF_NAME'`)*  
  Name of the dataframe, used in plot title.

- **`return_clusters`** *(default: `False`)*  
  If `True`, the function returns a tuple `(G, cluster_assignments)` where `cluster_assignments` maps each sample to its cluster ID. If `False`, only the graph `G` is returned.
---

### 📤 Output / Plot

- **Nodes**: Each point represents a sample, **colored by cluster status
  - 🟢 **Green**: Samples from clusters containing only samples with the same patient ID.
  - 🟣 **Magenta**: Samples from clusters with samples with different patient IDs.
  - 🟦 **Blue**: Samples that are the single representatives of their patient ID in the data and are singletons in the graph.
  - 🔴 **Pink**: Samples that have at least one other sample with the same patient ID in the data but are singletons in the graph.
- **Edges**:
  - 🟢 **Green**: Same-patient pairs that are nearest neighbors.
  - 🟣 **Magenta**: Cross-patient pairs that are nearest neighbors.
---


### 🧪 Example: Basic Usage in R

```r
res = spqrp$cluster_samples_iteratively(
  result = result,
  df = df,
  method = "UMAP",
  random_state = 42L,
  n_neighbors = 1L,
  max_component_size = 2L,
  n_neighbour_in_belonging = FALSE,
  n_umap_neighbors = 15L
)
g <- res[[1]]   # your Graph object from spqrp
coords_2d = res[[2]]

spqrp$plot_distances_neighbours_with_coloring_hue(
  df = df,
  G = g,
  coords_2d = coords_2d,
  method = "UMAP",
  subset_samples = NULL,
  highlight_singletons = TRUE,
  highlight_single_samples_missing_connections = TRUE,
  figsize = c(20, 20),
  label_font = 4.5,
  df_name = "Cohort A on Cohort A"
)
```
It shows the result in RStudio like <img src="https://github.com/user-attachments/assets/8a45e309-b77a-474d-972d-e61ba52424ca" width="500"/>
You can preview it directly in RStudio or save the plot with 
```r
plt <- import("matplotlib.pyplot")
plt$savefig("spqrp_graph.png", dpi=300, bbox_inches="tight")
```

### Visualizing a subset
To visualize only a subset of samples one can run the function with handing over ```subset_samples=list("9_0","23_1")``` to
```plot_distances_neighbours_with_coloring_hue```

<img src="https://github.com/user-attachments/assets/b4745354-a834-4678-8a29-f1c4a3ed1162" width="300"/>

### Additional Information

#### Note on UMAP with small sample sizes
UMAP relies on constructing a k-nearest-neighbor graph to learn the sample manifold. If the number of available samples is small (especially when fewer than ~10–20), this graph becomes unstable or degenerate. In such cases, UMAP may fail, produce warnings, or yield embeddings that are not meaningful. For very small datasets, alternative methods such as PCA or MDS are typically more robust.
