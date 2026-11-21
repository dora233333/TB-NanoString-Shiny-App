# TB NanoString Shiny App

This repository contains a Shiny application for dimension reduction and k-means clustering of TB NanoString gene-expression data.

## Features

The application allows the user to:

- Select PCA, UMAP, or NMF as the dimension-reduction method
- Select the number of k-means clusters
- View observations colored by k-means cluster
- View observations shaped by TB status
- View the k-means cluster centers
- View a cross-tabulation of clusters and TB status
- Calculate the adjusted Rand index

## Files

- `app.R`: Complete Shiny application
- `TBnanostring.rds`: TB NanoString gene-expression dataset
- `app_pca.png`: Screenshot of the PCA app
- `app_umap.png`: Screenshot of the UMAP app
- `app_nmf.png`: Screenshot of the NMF app

## Required R packages

The application requires the following R packages:

```r
install.packages(
  c(
    "shiny",
    "ggplot2",
    "uwot",
    "NMF",
    "mclust"
  )
)
```

The `NMF` package may also require the Bioconductor package `Biobase`:

```r
install.packages("BiocManager")

BiocManager::install(
  "Biobase"
)
```

## Running the app

1. Download or clone this repository.
2. Keep `app.R` and `TBnanostring.rds` in the same folder.
3. Open `app.R` in RStudio.
4. Install the required packages if necessary.
5. Click **Run App**.

The application calculates PCA, UMAP, and rank-2 NMF when it starts. The user can then select a dimension-reduction method and choose the number of k-means clusters.

## Output

The application displays:

- A two-dimensional representation of the TB NanoString data
- Points colored by k-means cluster
- Points shaped by TB status
- K-means cluster centers
- A cross-tabulation of clusters and TB status
- The adjusted Rand index
