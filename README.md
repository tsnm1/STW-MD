# STW-MD
STW-MD: A Novel Spatio-Temporal Weighting and Multi-Step Decision Tree Method for Considering Spatial Heterogeneity in Brain Gene Expression Data

## [UMAP_OLD3.ipynb](https://github.com/tsnm1/STW-MD/blob/main/UMAP_OLD3.ipynb "UMAP_OLD3.ipynb")
This code processes data (AD dataset) from nineteen brain regions and developmental information of study subjects. 

Initially, it employs the limma model for differential gene expression selection. 

Following this, it adjusts the weights and assigns weights to the top 5% of differentially expressed genes based on the absolute value of the fold change (FC). After weighting, it conducts a second round of screening using the limma model. 

Finally, it utilizes a multi-step decision tree model, including umap, k-means, pca and decision tree,  to further process the differentially expressed genes selected in the second round of screening.

## [UMAP_NEW3.ipynb](https://github.com/tsnm1/STW-MD/blob/main/UMAP_NEW3.ipynb "UMAP_NEW3.ipynb")
This code processes data (brain development dataset) from sixteen brain regions and patient information data. 

Firstly, it conducts descriptive analysis on the data. 

Subsequently, it employs the limma model for differential gene expression selection. After adjusting the weights, it applies a second round of screening using the limma model. 

Finally, it employs a multi-step decision tree model, including umap, k-means, pca and decision tree, to further process the differentially expressed genes selected in the second round of screening.

## [R-code2.0.R](https://github.com/tsnm1/STW-MD/blob/main/R-code2.0.R "R-code2.0.R")

Part 1: This code segment preprocesses the brain development data. It starts by connecting the data, obtaining the source of target gene data, as well as information on brain regions, samples, etc. Then, it filters the gene expression data and separates them based on brain regions.

Part 2: Visualizing the results of gene enrichment analysis.

Part 3: Visualizing the number of upregulated and downregulated genes.







