# What is this?
This folder contains R files of the entire analysis, performed sequentially from loading raw 10X Genomics files to generating an integrated analysis profile.

### 1_environment_setup.R:
This file contains R script for installing the correct programs utilized by R to effectively perform the scRNAseq analysis. These software include:
1. devtools
2. pkgbuild
3. Rtools
4. pkgbuild
5. Matrix
6. ggridges
7. cowplot
8. ggrepel
9. R.utils
10. gridExtra
11. plotly
12. Seurat
13. Bioconductor/Biocmanager
14. Monocle


### 2_alk3n3_seurat_setup.R:
This file contains the R script required to setup a Seurat v3.0.0 object which will be later utilized during the downstream analysis. 

### 3_alk3n3_analysis.R:
A complete script of analysis performed on the entire alk3n3 dataset (alk3n3 or alk3_n=3 represents the n=3 donors analyzed in this code)

### 4_alk3n3NOMES_analysis.R:
A set of R code which allows for mesenchymal clusters to be subsetted out of the analysis. The resulting dataset, is what has been used throughout the analysis. 

### 5_alk3n3_integration_with_cellatlas.R:
A set of code which allows for the integration of the alk3n3 dataset into a islet pancreas single cell RNAseq gene expression atlas.

### 6_alk3n3_integration_visualization.R:
A set of code allowing for the visualization of the integrative analysis performed previously.

#### 7_alk3n3_integration_3Dvisualization.R:
A set of code used to allow the viewing of 3D tSNE visuals using the RShiny derived package plotly.
