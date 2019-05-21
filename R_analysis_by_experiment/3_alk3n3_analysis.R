# DESCRIPTION ####
# The following code ASSUMES all dependencies in R  have been installed and previous R code have been run successfully see file(s): 
# 1_environment_setup.R 
# 2_alk3n3_seurat_setup.R
# The purpose of this code is to analyze n=3 data, as a Seurat object.

# R v3.5.3 (x64 bit) and RStudio v1.2.1335 (x64 bit) were used for running this code (see readme on how to install)
# Seurat is a multimodal single Cell RNA seq analysis algorithm created by
# The Satija Lab. For more information please see: https://satijalab.org/seurat/

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
library(plot_ly)

# CONFIRM CORRECT INSTALL ####
# Confirm package version of Seurat and Monocle
packageVersion("Seurat")
packageVersion("monocle")

# LINEAR DIMENSIONALITY REDUCTION ####
# this Seurat object was created and normalized using code outlined in: 2_alk3n3_seurat_object_setup_05202019.R
# This Seurat object is called alk3n3.combined
# Without running the R files in sequence YOU WILL NOT get the data output. Running in sequence is essential.

# PCA analysis
alk3n3.combined <- RunPCA(object = alk3n3.combined, features = VariableFeatures(object = alk3n3.combined))

# Examine and visualize PCA results a few different ways
print(x = alk3n3.combined[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object = alk3n3.combined, dims = 1:2, reduction = "pca")
DimPlot(object = alk3n3.combined, reduction = "pca")
DimHeatmap(object = alk3n3.combined, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(object = alk3n3.combined, dims = 1:15, cells = 500, balanced = TRUE)

# Using JacStraw plots to visualize LDR measurements
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
alk3n3.combined <- JackStraw(object = alk3n3.combined, num.replicate = 100)
alk3n3.combined <- ScoreJackStraw(object = alk3n3.combined, dims = 1:20)
JackStrawPlot(object = alk3n3.combined, dims = 1:15)

# Elbow plot
ElbowPlot(object = alk3n3.combined)

# Clustering
alk3n3.combined <- FindNeighbors(object = alk3n3.combined, dims = 1:10)
alk3n3.combined <- FindClusters(object = alk3n3.combined, resolution = 0.4)

# Look at cluster IDs of the first 5 cells
head(x = Idents(object = alk3n3.combined), 5)

# NON-LINEAR DIMENSIONALITY REDUCTION ####
# Using tSNE for dimensionality eduction
alk3n3.combined <- RunTSNE(object = alk3n3.combined, dims = 1:10)

# Note: you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(object = alk3n3.combined, reduction = "tsne", pt.size = 2)

# DIFFERENTIALLY EXPRESSED-GENE ANALYSIS ####
# Find markers for every cluster compared to all remaining cells, report only the positive ones
# Here we defin a DE gene as a gene which has:
# Fold Change of >2
# Atleast 50% of cells express that gene
# p value < 0.001
alk3n3.combined.markers <- FindAllMarkers(object = alk3n3.combined, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.7)
alk3n3.combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(alk3n3.combined.markers, 'alk3n3.combined.markers.csv')

# You can look at the expression of a particular gene across an entire data set here
VlnPlot(object = alk3n3.combined, features = c("COL1A1"))

# Look at a heatmap of top 10 most differentially expressed genes
top10 <- alk3n3.combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(object = alk3n3.combined, features = top10$gene) + NoLegend()

################## #
################## #
# CODE END ####
