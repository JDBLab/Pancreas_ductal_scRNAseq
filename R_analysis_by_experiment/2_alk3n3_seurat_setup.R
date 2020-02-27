# DESCRIPTION ####
# The following code ASSUMES all dependencies in R  have been installed See file(s): 
# 1_environment_setup.R)
# The purpose of this code is to setup+normalize+threshold data, as a Seurat object.

# R v3.5.3 (x64 bit) and RStudio v1.2.1335 (x64 bit) were used for running this code (see readme on how to install)
# Seurat is a multimodal single Cell RNA seq analysis algorithm created by
# The Satija Lab. For more information please see: https://satijalab.org/seurat/
# Code run-time 4m 25s till rds save, 5m 35s total.

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

# CONFIRM CORRECT INSTALL ####
# Confirm package version of Seurat and Monocle
packageVersion("Seurat")

# SETTING UP SEURAT OBJECT ####
# Data-integration
# Increase global calculation ####
options(future.globals.maxSize = 4000 * 1024^2)

# Loading data ####
HPD1.data <- Read10X(data.dir = "C:/Users/mqadir/Desktop/alk3n3 scRNAseq data/HPD1/filtered_feature_bc_matrix")
HPD1 <- CreateSeuratObject(counts = HPD1.data, 
                           project = "HPD1"
)
HPD2.data <- Read10X(data.dir = "C:/Users/mqadir/Desktop/alk3n3 scRNAseq data/HPD2/filtered_feature_bc_matrix")
HPD2 <- CreateSeuratObject(counts = HPD2.data,
                           project = "HPD2"
)
HPD3.data <- Read10X(data.dir = "C:/Users/mqadir/Desktop/alk3n3 scRNAseq data/HPD3/filtered_feature_bc_matrix")
HPD3 <- CreateSeuratObject(counts = HPD3.data,
                           project = "HPD3"
)

# Adding Cell ids ####
head(x = colnames(x = HPD1))
HPD1 <- RenameCells(object = HPD1, add.cell.id = "HPD1")
head(x = colnames(x = HPD1))
head(x = colnames(x = HPD2))
HPD2 <- RenameCells(object = HPD2, add.cell.id = "HPD2")
head(x = colnames(x = HPD2))
head(x = colnames(x = HPD3))
HPD3 <- RenameCells(object = HPD3, add.cell.id = "HPD3")
head(x = colnames(x = HPD3))

# THRESHOLDING ####
# The [[ operator can add columns to object metadata. Here we store MT percentages for each cell
HPD1[["percent.mt"]] <- PercentageFeatureSet(object = HPD1, pattern = "^MT-")
HPD2[["percent.mt"]] <- PercentageFeatureSet(object = HPD2, pattern = "^MT-")
HPD3[["percent.mt"]] <- PercentageFeatureSet(object = HPD3, pattern = "^MT-")

# Visualize QC metrics as a violin plot
p1 <- VlnPlot(object = HPD1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p2 <- VlnPlot(object = HPD2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p3 <- VlnPlot(object = HPD3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CombinePlots(plots = list(p1, p2, p3), ncol =1)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
p4 <- FeatureScatter(object = HPD1, feature1 = "nCount_RNA", feature2 = "percent.mt")
p5 <- FeatureScatter(object = HPD1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p6 <- FeatureScatter(object = HPD2, feature1 = "nCount_RNA", feature2 = "percent.mt")
p7 <- FeatureScatter(object = HPD2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p8 <- FeatureScatter(object = HPD3, feature1 = "nCount_RNA", feature2 = "percent.mt")
p9 <- FeatureScatter(object = HPD3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(p4, p5, p6, p7, p8, p9), ncol =2)

# MT based cell thresholding
HPD1 <- subset(x = HPD1, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
HPD2 <- subset(x = HPD2, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
HPD3 <- subset(x = HPD3, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)

# Data after thresholding
# Visualize QC metrics as a violin plot
p10 <- VlnPlot(object = HPD1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p11 <- VlnPlot(object = HPD2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p12 <- VlnPlot(object = HPD3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CombinePlots(plots = list(p10, p11, p12), ncol =1)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
p13 <- FeatureScatter(object = HPD1, feature1 = "nCount_RNA", feature2 = "percent.mt")
p14 <- FeatureScatter(object = HPD1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p15 <- FeatureScatter(object = HPD2, feature1 = "nCount_RNA", feature2 = "percent.mt")
p16 <- FeatureScatter(object = HPD2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p17 <- FeatureScatter(object = HPD3, feature1 = "nCount_RNA", feature2 = "percent.mt")
p18 <- FeatureScatter(object = HPD3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(p13, p14, p15, p16, p17, p18), ncol =2)

# Create a pancreas-list ####
alk3n3.list <- c(HPD1, HPD2, HPD3)

# run SCTransform on EACH OBJECT SEPERATELY ####
for (i in 1:length(alk3n3.list)) {
  alk3n3.list[[i]] <- SCTransform(alk3n3.list[[i]], variable.features.n = 3000, verbose = FALSE)
}

# Selecting features for downstream integration ####
alk3n3.features <- SelectIntegrationFeatures(object.list = alk3n3.list, nfeatures = 3000)
alk3n3.list <- PrepSCTIntegration(object.list = alk3n3.list, anchor.features = alk3n3.features, 
                                    verbose = FALSE)

# Anchor identification and integration ####
alk3n3.anchors <- FindIntegrationAnchors(object.list = alk3n3.list, normalization.method = "SCT", 
                                           anchor.features = alk3n3.features, verbose = FALSE)

alk3n3.integrated <- IntegrateData(anchorset = alk3n3.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)

# After integration the RNA and SCT slot no longer contain variable or SCTransformned data, so re-analyze enabling visualization
# Look at your default assay
DefaultAssay(object = alk3n3.integrated)

# Change default assay to RNA, save information in the "RNA" assay
DefaultAssay(object = alk3n3.integrated) <- "RNA"

# Find variable features and save them in the RNA assay slot
alk3n3.integrated <- FindVariableFeatures(object = alk3n3.integrated, selection.method = "vst", nfeatures = 3000, assay = "RNA")

# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = alk3n3.integrated), 10)

# Ploting variable features with and without labels
p19 <- VariableFeaturePlot(object = alk3n3.integrated)
p19
LabelPoints(plot = p19, points = top10, repel = TRUE)

# Check data
head(x = colnames(x = alk3n3.integrated))
tail(x = colnames(x = alk3n3.integrated))
unique(x = sapply(X = strsplit(x = colnames(x = alk3n3.integrated), split = "_"), FUN = "[", 1))
table(alk3n3.integrated$orig.ident)

# Final figures for supplementary figures in Qadir et al., 20??
# Visualize QC metrics as a violin plot
VlnPlot(object = alk3n3.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, cols = c("seagreen",
                                                                                                                 "deepskyblue4",
                                                                                                                 "firebrick2"))
FeatureScatter(object = alk3n3.integrated, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(object = alk3n3.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Save RDS file
# saveRDS(alk3n3.integrated, file = "C:/Users/mxq52/Box/NPOD PAPER/Round 4/R_files/alk3n3.integrated.code2.rds")

################## #
################## #
# CODE END ####
