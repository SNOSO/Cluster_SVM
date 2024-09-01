# README for Spatial Transcriptomics SVM Analysis

## Overview
This project contains R functions designed to load, process, and analyze spatial transcriptomics data using Support Vector Machines (SVM). The analysis focuses on the classification of spatial clusters within different datasets (DLPFC, CODEX, Xenium) and evaluating the performance of the SVM model through metrics like the Adjusted Rand Index (ARI) and Confusion Matrix.

## Prerequisites
Ensure that the following R packages are installed:
- `e1071`: For running SVM models.
- `mclust`: For calculating the Adjusted Rand Index (ARI).
- `caret`: For generating Confusion Matrices.

```r
install.packages(c("e1071", "mclust", "caret"))
```

## Functions

### 1. `load_and_process_data`
This function loads and processes spatial transcriptomics data based on the specified dataset type. The function supports three data types: DLPFC, CODEX, and Xenium.

#### Parameters:
- `data_type`: A string indicating the type of data to load (`"DLPFC"`, `"CODEX"`, or `"Xenium"`).
- `gene_list_path`: Path to the RDS file containing the list of genes to subset.
- `data_path`: Path to the RDS file containing the dataset to be loaded.

#### Returns:
A list containing:
- `matrix`: The processed expression matrix.
- `coordinates`: The spatial coordinates corresponding to the expression data.
- `labels`: The cell type or niche labels associated with the data.

#### Example Usage:
```r
data_list <- load_and_process_data(data_type = "DLPFC", 
                                   gene_list_path = "~/path/to/DLPFCGeneList.RDS", 
                                   data_path = "~/path/to/SpatialSeurats.RDS")
```

### 2. `RunSVM`
This function applies an SVM model to the spatial transcriptomics data to classify spatial clusters. It supports neighborhood averaging of the input data and allows for flexible SVM kernel selection.

#### Parameters:
- `matrix`: The expression matrix to be used for training and testing.
- `coordinates`: The spatial coordinates associated with the expression matrix.
- `labels`: The cell type or niche labels.
- `nsize`: The neighborhood size for averaging (default is 3).
- `kernel_type`: The type of SVM kernel to use (`"linear"`, `"polynomial"`, `"radial"`, or `"sigmoid"`).

#### Returns:
A list containing:
- `ARI`: The Adjusted Rand Index for the model's performance.
- `ConfusionMatrix`: The Confusion Matrix comparing the predicted labels to the true labels.

#### Example Usage:
```r
SVM <- RunSVM(matrix = data_list$matrix, 
              coordinates = data_list$coordinates, 
              labels = data_list$labels, 
              nsize = 3, 
              kernel_type = "radial")
```

## Example Workflow

1. **Load and Process Data:**
   Load the DLPFC data, subset it based on a gene list, and extract relevant coordinates and labels.
   ```r
   data_list <- load_and_process_data(data_type = "DLPFC", 
                                      gene_list_path = "~/path/to/DLPFCGeneList.RDS", 
                                      data_path = "~/path/to/SpatialSeurats.RDS")
   ```

2. **Run SVM:**
   Apply an SVM model with a radial kernel to classify spatial clusters within the processed data.
   ```r
   SVM <- RunSVM(matrix = data_list$matrix, 
                 coordinates = data_list$coordinates, 
                 labels = data_list$labels, 
                 nsize = 3, 
                 kernel_type = "radial")
   ```

3. **Evaluate Results:**
   Access the ARI and Confusion Matrix from the SVM model's output.
   ```r
   print(SVM$ARI)
   print(SVM$ConfusionMatrix)
   ```

## Notes
- Ensure the paths provided to the functions are accurate and that the required RDS files are available.
- The functions are flexible and can be adjusted to work with other spatial transcriptomics datasets by modifying the data loading and processing steps.
