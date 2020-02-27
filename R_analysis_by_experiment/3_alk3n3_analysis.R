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
library(clustree)
library(plot_ly)

# CONFIRM CORRECT INSTALL ####
# Confirm package version of Seurat and Monocle
packageVersion("Seurat")
packageVersion("monocle")

# LINEAR DIMENSIONALITY REDUCTION ####
# this Seurat object was created and normalized using code outlined in: 2_alk3n3_seurat_object_setup_05202019.R
# This Seurat object is called alk3n3.integrated
# Without running the R files in sequence YOU WILL NOT get the data output. Running in sequence is essential.
# Look at your default assay
DefaultAssay(object = alk3n3.integrated)

# Change default assay to integrated, save information in the "integrated" assay
DefaultAssay(object = alk3n3.integrated) <- "integrated"

# PCA analysis data will be stored in the "reductions' slot
alk3n3.integrated <- RunPCA(object = alk3n3.integrated, features = VariableFeatures(object = alk3n3.integrated))

# Examine and visualize PCA results a few different ways
print(x = alk3n3.integrated[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object = alk3n3.integrated, dims = 1:2, reduction = "pca")
DimPlot(object = alk3n3.integrated, reduction = "pca")
DimHeatmap(object = alk3n3.integrated, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(object = alk3n3.integrated, dims = 1:15, cells = 500, balanced = TRUE)

# Using JacStraw plots to visualize LDR measurements
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
# alk3n3.integrated <- JackStraw(object = alk3n3.integrated, num.replicate = 100)
# alk3n3.integrated <- ScoreJackStraw(object = alk3n3.integrated, dims = 1:10)
# JackStrawPlot(object = alk3n3.integrated, dims = 1:10)

# Elbow plot
ElbowPlot(object = alk3n3.integrated)

# Clustering
alk3n3.integrated <- FindNeighbors(object = alk3n3.integrated, dims = 1:20)
alk3n3.integrated <- FindClusters(object = alk3n3.integrated, resolution = 0.4)

# NON-LINEAR DIMENSIONALITY REDUCTION ####
# RunUMAP
# reticulate::py_install(packages = 'umap-learn')

# We ran RunUMAP on the parameters: metric = cosine and umap.method = UWOT
# you are welcome to try UMAp-learn, the interpretation of the outcome doesnt
# adversally affect the cell clustering on superficial observation
alk3n3.integrated <- RunUMAP(alk3n3.integrated, dims = 1:20, metric = 'cosine', umap.method = 'uwot')

# Visualization
# Note: you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
# UMAP Visualization
p20 <- DimPlot(alk3n3.integrated, group.by = c("orig.ident", "integrated_snn_res.0.4"), combine = FALSE, pt.size = 1)
p20 <- lapply(X = p20, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3))))
CombinePlots(p20)

# DIFFERENTIALLY EXPRESSED-GENE ANALYSIS ####
# Find markers for every cluster compared to all remaining cells, report only the positive ones
# Here we define a DE gene as a gene which has:
# Fold Change of >1.5
# Atleast 10% of cells express that gene
# p value < 0.001
alk3n3.integrated.markers <- FindAllMarkers(object = alk3n3.integrated, only.pos = TRUE, logfc.threshold = 0.41, slot = 'data', test.use = 'wilcox')
alk3n3.integrated.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(alk3n3.integrated.markers, 'C:/Users/mxq52/Box/NPOD PAPER/Round 4/R_files/alk3n3.integrated.markers.csv')

# Look at your default assay
DefaultAssay(object = alk3n3.integrated)

# Change default assay to RNA, to visualize data on graph
# You can toggle between integrated, SCT and RNA to see different expression profiles/different normalizations
DefaultAssay(object = alk3n3.integrated) <- "RNA"

# You can look at the expression of a particular gene across an entire data set here
# As an example we need to remove all mesenchymal cells from the analysis
# We use COL1A1 and THY1 as mesenchymal identifiers
# Violin plot
VlnPlot(object = alk3n3.integrated, features = c("THY1"), group.by = "integrated_snn_res.0.4", ncol = 1)

# UMAP expression plot
FeaturePlot(object = alk3n3.integrated, 
            features = c("CPA2"),
            pt.size = 2,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            max.cutoff = 100,
            order = TRUE)


# Look at a heatmap of top 10 most differentially expressed genes
# Change default assay to RNA, save information in the "RNA" assay
DefaultAssay(object = alk3n3.integrated) <- "integrated"

# Create heatmap using doheatmap
top100.all <- alk3n3.integrated.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
DoHeatmap(object = alk3n3.integrated, features = top100.all$gene) + NoLegend()

# Save RDS file
# saveRDS(alk3n3.integrated, file = "C:/Users/mxq52/Box/NPOD PAPER/Round 4/R_files/alk3n3.integrated.code3.rds")

################## #
################## #
# CODE END ####
