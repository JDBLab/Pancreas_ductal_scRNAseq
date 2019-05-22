# Single cell resolution analysis of the human pancreatic ductal progenitor cell niche

## What is this?
This repository contains coding scripts utilized for the analysis performed in the "Single cell resolution analysis of the human pancreatic ductal progenitor cell niche" publication (Qadir/Alvarez-Cubela et. al, 2019). The purpose of providing the code here is to allow for transparency and robust data-analysis reproducibility. Most of the steps used for data analysis and visualization have been optimised for an average computing environment (for the year 2019). Some analyses however, require a high-performace computing environment (see computing environment). The methodology has already been described extensively in the manuscript. However, this analysis relies heavily on scRNAseq analysis algorithms developed by the Satija lab, namely Seurat (for a complete list of dependencies and code utilized see analysis & visualization programs).

## Downloading Data files
Data files utilized in this analysis have been deposited in the Gene Expression Omnibus (GEO), gene expression data repository at the NIH. Data are part of the GSEXXXXXX repository and can be found here: (). please note, since this manuscript is currently under review, GSE raw data files will not be present for public download. If you would like to download the data, please email the project leader for requests. 

### Data sub-structure
We povide raw FASTQ files generated of single-cell cDNA libraries sequenced by the Illumina sequencing platform, filtered/unfiltered post-alignment count files generated by the Cellranger v3.0.1 software. In addition we also provide a gene expression matrix containing data from the [GSE81076](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81076) [(Grun et. al, 2016: Cell Stem Cell)](https://www.sciencedirect.com/science/article/pii/S1934590916300947?via%3Dihub) and [GSE85241](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85241) [(Muraro et. al, 2016: Cell Systems)](https://www.sciencedirect.com/science/article/pii/S2405471216302927?via%3Dihub) datasets. These two datasets represent the human pancreas atlas for single-cell gene transcription. We use these datasets for our integrated single-cell analysis.

#### FASTQ files
These are sequencing reads generated by the Illumina sequencing platform. Files contain raw reads and sequencing efficiency information.
These are the input files for the Cellranger software.

#### Cellranger output (processed gene-counts of single cells/barcodes)
This contains data outputs of Cellranger v3.0.1, which was run using default settings. Code used to analyze data is a part of this repository. This data contains filtered/unfiltered count files for gene expression across barcodes/cells. 

#### Data analysis of ALK3n3 dataset
Preliminary data-analyses involving n=3 de-identified human exocrine-pancreata derived ALK3+ cells are included in this file. This includes data thresholding, normalization, subsetting, linear dimensionality reduction (PCA), non-linear multimodal dimensionality reduction (PCA/tSNE), clustering, and data visualization.

#### Data analysis of integrated-dataset
In order to understand where our cells map against other pancretic cells, we mapped our cells against a human pancreas dataset for single-cell gene transcription. In doing so, we were able to understand where our cells reside in the context of other neighbouring pancreatic cells. We hope to use our data to expand the human pancreas single-cell transcriptional profile. These data were downloaded from a pre-analyzed gene expression matrix created as part of Seurat tutorial, which can be found and [downloaded from here](https://satijalab.org/seurat/v3.0/pancreas_integration_label_transfer.html).

## Analysis and visualization programs
### Cellranger software from 10X Genomics:
1. [Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation)

### R and R's integrated developmental environment RStudio:
1. [R v3.5.3 (x64 bit)](https://cran.r-project.org/bin/windows/base/old/)
2. [RStudio v1.2.1335 (x64 bit)](https://www.rstudio.com/products/rstudio/download/)
3. [RTools v3.3.X](https://cran.r-project.org/bin/windows/Rtools/index.html)
4. [Tutorial for R](https://cran.r-project.org/doc/manuals/r-release/R-intro.html)
5. [Tutorial for RStudio](https://resources.rstudio.com/)

### scRNAseq analysis pipeline SEURAT developed by the Satija lab:
1. [Source code for Seurat v3.0.0](https://cran.r-project.org/web/packages/Seurat/index.html)
2. [Tutorials for Seurat](https://satijalab.org/seurat/)

### Pseudotemporal gene expression analysis using Monocle developed by the Trapnell Lab:
1. [Source code for Monocle v2.8.0](https://bioconductor.org/packages/release/bioc/html/monocle.html)
2. [Tutorial for Monocle](http://cole-trapnell-lab.github.io/monocle-release/docs/#constructing-single-cell-trajectories)

### 3D tSNE analysis and Gene expression plotting
1. [Source code for plotly](https://cran.r-project.org/web/packages/plotly/index.html)
2. [Tutorial for plotly](https://plot.ly/r/)

## Setting up the right environment
under construction.

## Computing environment
### Hardware and OS environment for running Cellranger
1. Processor: Intel Sandy Bridge E5-2670 (16cores x 16 threads)
2. RAM: 25GB
3. OS: CentOS 6.5

Hardware integrated into the [Pegasus Supercomputer](http://ccs.miami.edu/ac/service/pegasus/) at the University of Miami

### Hardware and OS environment for running Seurat and Monocle
1. Processor: Intel Core i7-6700 CPU (4cores x 8threads)
2. RAM: 32GB DDR3
3. OS: Windows 10 Enterprise (x64 bit)

### Hardware and OS environment for running Velocyto
Under construction
