# DESCRIPTION ####
# The following code ASSUMES all dependencies in R  have been installed See file(s): 
# 1_environment_setup.R)
# The purpose of this code is to setup+normalize+threshold data, as a Seurat object.

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

# CONFIRM CORRECT LOADING ####
# Confirm package version of Seurat and Monocle
packageVersion("Seurat")
packageVersion("monocle")

# SETTING UP SEURAT OBJECT ####
# Data-merge
HPD1.data <- Read10X(data.dir = "C:/Users/mxq52/Box/Fahd Shared JUAN/Dominguez-Bendala, Juan scRNAseq prj/10x Gen input files/HPD1/filtered_feature_bc_matrix")
HPD1 <- CreateSeuratObject(counts = HPD1.data, 
                             project = "HPD1"
)

HPD2.data <- Read10X(data.dir = "C:/Users/mxq52/Box/Fahd Shared JUAN/Dominguez-Bendala, Juan scRNAseq prj/10x Gen input files/HPD2/filtered_feature_bc_matrix")
HPD2 <- CreateSeuratObject(counts = HPD2.data,
                             project = "HPD2"
)

HPD3.data <- Read10X(data.dir = "C:/Users/mxq52/Box/Fahd Shared JUAN/Dominguez-Bendala, Juan scRNAseq prj/10x Gen input files/HPD3/filtered_feature_bc_matrix")
HPD3 <- CreateSeuratObject(counts = HPD3.data,
                             project = "HPD3"
)

alk3n3.combined <- merge(x = HPD1, y = c(HPD2, HPD3), add.cell.ids = c("HPD1", "HPD2", "HPD3"), project = "alk3n3")
alk3n3.combined

# Check data
head(x = colnames(x = alk3n3.combined))
tail(x = colnames(x = alk3n3.combined))
unique(x = sapply(X = strsplit(x = colnames(x = alk3n3.combined), split = "_"), FUN = "[", 1))
table(alk3n3.combined$orig.ident)

# NORMALIZATION AND THRESHOLDING ####
# The [[ operator can add columns to object metadata. Here we store MT percentages for each cell
alk3n3.combined[["percent.mt"]] <- PercentageFeatureSet(object = alk3n3.combined, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(object = alk3n3.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

FeatureScatter(object = alk3n3.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(object = alk3n3.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# MT based cell thresholding
alk3n3.combined <- subset(x = alk3n3.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 5)

# Normalization
alk3n3.combined <- NormalizeData(object = alk3n3.combined)

# Identifying highly variable genes
alk3n3.combined <- FindVariableFeatures(object = alk3n3.combined, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = alk3n3.combined), 10)

# Ploting variable features with and without labels
VariableFeaturePlot(object = alk3n3.combined)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# Data scaling against percentage of expressed mitochondrial genes
# This step is computationally intensive, its a good idea to save the global environment on your PC once it has run.
# OR save as a R project, so that you dont haveto re-run the analysis again and again
alk3n3.combined <- ScaleData(object = alk3n3.combined, vars.to.regress = 'percent.mt', 
                             features = rownames(alk3n3.combined), block.size = 33000, 
                             do.par = TRUE, num.cores = 8)
# Save RDS file
saveRDS(alk3n3.combined, file = "~/alk3n3.combined.rds")

################## #
################## #
# CODE END ####
