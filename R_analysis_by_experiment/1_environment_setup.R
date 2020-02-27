# DESCRIPTION ####
# The following code will install all dependencies in R allowing for analysis to be performed on the alk3n3 dataset

# In order to utilize umap, which will work within a python v3.7 environment, install Anaconda3 in windows
# Install and run R v3.5.3 (x64 bit) and RStudio v1.2.1335 (x64 bit) for running this code (see readme on how to install)
# Seurat is a multimodal single Cell RNA seq analysis algorithm created by
# The Satija Lab. For more information please see: https://satijalab.org/seurat/
# Monocle is a psudotemporal differentiation trajectory analysis program crated by
# The Trapnell Lab. For more information please see: http://cole-trapnell-lab.github.io/monocle-release/

# INSTALLATION ####
# Please download and install Rtools 3.5 from http://cran.r-project.org/bin/windows/Rtools/
install.packages("devtools")
install.packages("pkgbuild")

library(devtools)
library(pkgbuild)
find_rtools() # This should come back as true

# Install additional packages
install.packages("Matrix")
install.packages("ggridges")
install.packages("cowplot")
install.packages('ggrepel')
install.packages("R.utils")
install.packages("gridExtra")
install.packages("Seurat")
install.packages("plotly")
install.packages("clustree")
install.packages('multtest')

# Install the scRNAseq analysis package Seurat
install.packages('Seurat') #This installs Seurat v3.0.0 as of 5/20/2019, if Seurat is updated direct installation to v3.0.0

# Install the scRNAseq pseudotemporal analysis package Monocle
# source("http://bioconductor.org/biocLite.R") #Monocle v2.8.0 is used
# biocLite()
# biocLite("monocle")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("monocle")

# Install multtest for Seurat
BiocManager::install("multtest")

# LOAD LIBRARIES ####
# Restart Rstudio or R
# Run the following code once you have Seurat installed

library(ggplot2)
library(cowplot)
library(Matrix)
library(ggridges)
library(ggrepel)
library(dplyr)
library(Seurat)
library(monocle)
library(plotly)
library(clustree)

# CONFIRM CORRECT INSTALL ####
# Confirm package version of Seurat and Monocle
packageVersion("Seurat")
packageVersion("monocle")


################## #
################## #
# CODE END ####
