# DESCRIPTION ####
# The following code ASSUMES all dependencies in R  have been installed and previous R code have been run successfully see file(s): 
# 1_environment_setup.R 
# The purpose of this code is integrate the alk3n3 dataset with a human pancreas single cell atlas. This data has been utlized from 
# the following datasets: GSE81076 and GSE85241. A gene expression matrix of this dataset exist, as part of a Seurat v3.0.0 tutorial
# those files can be found here https://satijalab.org/seurat/v3.0/pancreas_integration_label_transfer.html
# These files are necessary for the analysis. We will attempt to upload them to GEO for asy access.

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
library(plotly)
library(clustree)

# CONFIRM CORRECT INSTALL ####
# Confirm package version of Seurat and Monocle
packageVersion("Seurat")
packageVersion("monocle")

# Set global environment parameter
options(future.globals.maxSize = 4000 * 1024^2)

# OBJECT SETUP AND NORMALIZATION ####
# Load the Human Pancreas Single-cell RNAseq Atlas.
pancreas.data <- readRDS(file = "C:/Users/mxq52/Box/NPOD PAPER/Rounds 1-3/NEW R files/pancreas_v3_files/pancreas_expression_matrix.rds")
metadata <- readRDS(file = "C:/Users/mxq52/Box/NPOD PAPER/Rounds 1-3/NEW R files/pancreas_v3_files/pancreas_metadata.rds")

# STEP 1: Add sample specific and experiment specific metadata ####
# Stage 1 is for the Muraro-Enge dataset
# Extract information from Muraro and Enge datasets
# Create Seurat obj: celseq and celseq2
pancreas <- CreateSeuratObject(counts = pancreas.data, meta.data = metadata)
pancreas.list <- SplitObject(object = pancreas, split.by = "tech")

# Extract and Organise Celseq dataset
celseq <- pancreas.list[["celseq"]]

# Sample specific Metadata addition
celseq$sample <- "celseq" 
celseq$exp <- "Muraro-Grun"

# Extract and Organise Celseq2 dataset
celseq2 <- pancreas.list[["celseq2"]]

# Sample specific Metadata addition
celseq2$sample <- "celseq2" 
celseq2$exp <- "Muraro-Grun"

# Stage 2 is for the Qadir dataset
# Create Seurat obj: HPD1
HPD1.data <- Read10X(data.dir = "C:/Users/mxq52/Desktop/cell_ranger/HPD1/filtered_feature_bc_matrix")
HPD1 <- CreateSeuratObject(counts = HPD1.data, 
                           project = "HPD1"
)

# Sample specific Metadata addition
HPD1$sample <- "HPD1"
HPD1$exp <- "Qadir"

# Create Seurat obj: HPD3
HPD2.data <- Read10X(data.dir = "C:/Users/mxq52/Desktop/cell_ranger/HPD2/filtered_feature_bc_matrix")
HPD2 <- CreateSeuratObject(counts = HPD2.data,
                           project = "HPD2"
)

# Sample specific Metadata addition
HPD2$sample <- "HPD2"
HPD2$exp <- "Qadir"

# Create Seurat obj: HPD3
HPD3.data <- Read10X(data.dir = "C:/Users/mxq52/Desktop/cell_ranger/HPD3/filtered_feature_bc_matrix")
HPD3 <- CreateSeuratObject(counts = HPD3.data,
                           project = "HPD3"
)

# Sample specific Metadata addition
HPD3$sample <- "HPD3"
HPD3$exp <- "Qadir"

#################### #
#################### #

# STEP 2: Thresholding ####
# The operator can add columns to object metadata. This is a great place to stash QC stats
HPD1[["percent.mt"]] <- PercentageFeatureSet(object = HPD1, pattern = "^MT-")
HPD2[["percent.mt"]] <- PercentageFeatureSet(object = HPD2, pattern = "^MT-")
HPD3[["percent.mt"]] <- PercentageFeatureSet(object = HPD3, pattern = "^MT-")

# RNA based cell thresholding
HPD1 <- subset(x = HPD1, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
HPD2 <- subset(x = HPD2, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
HPD3 <- subset(x = HPD3, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
celseq <- subset(x= celseq, subset = nFeature_RNA >200 & nFeature_RNA <10000)
celseq2 <- subset(x= celseq2, subset = nFeature_RNA >200 & nFeature_RNA <10000)

# Visualize QC metrics as a violin plot
VlnPlot(object = HPD1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HPD2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HPD3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = celseq, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(object = celseq2, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

#################### #
#################### #

# Step 3: Add cell IDs ####
# Add cell IDs
celseq <- RenameCells(celseq, add.cell.id = "celseq")
celseq2 <- RenameCells(celseq2, add.cell.id = "celseq2")
HPD1 <- RenameCells(HPD1, add.cell.id = "HPD1")
HPD2 <- RenameCells(HPD2, add.cell.id = "HPD2")
HPD3 <- RenameCells(HPD3, add.cell.id = "HPD3")

#################### #
#################### #

# Step 4: Merge Datasets
# Merge refrence-qadir datasets
pancreas.list <- list("HPD1" = HPD1, "HPD2" = HPD2, "HPD3" = HPD3, "celseq" = celseq, "celseq2" = celseq2)

# Step 5: Data normalization
# Normalize the dataset using SCTransform
for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- SCTransform(pancreas.list[[i]], verbose = TRUE)
}

#################### #
#################### #

# Step 6: Feature selection
# Select features for downstream integration
pancreas.features <- SelectIntegrationFeatures(object.list = pancreas.list, nfeatures = 3000)
pancreas.list <- PrepSCTIntegration(object.list = pancreas.list, anchor.features = pancreas.features, 
                                    verbose = TRUE)

# Organize a archetectural map of the pancreas derived from Muraro and Enge we choose datasets 4 and 5
pancreas.list
reference_dataset <- c(4, 5)

#################### #
#################### #

# Step 7: Anchor identification and data integration
# Identify anchors and integrate dataset
pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, normalization.method = "SCT", 
                                       anchor.features = pancreas.features, reference = reference_dataset)

pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT", 
                                     verbose = TRUE)

#################### #
#################### #

# Step 8: Linear dimensionality assessment
# Look at your default assay
DefaultAssay(object = pancreas.integrated)

# Change default assay to integrated, to view dimensionality
DefaultAssay(object = pancreas.integrated) <- "integrated"

# Dimensionality assessment
pancreas.integrated <- RunPCA(pancreas.integrated, features = VariableFeatures(object = pancreas.integrated))

# Examine data dimensionality
ElbowPlot(pancreas.integrated)

#################### #
#################### #

# Step 9: CLUSTER ANALYSIS ####
# Clustering, we run multiple permutations to allow clustree to analyze optimal clustering resolution.
pancreas.integrated <- FindNeighbors(object = pancreas.integrated, dims = 1:10)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.1)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.2)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.3)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.4)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.5)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.6)

# Cluster-tree analysis, looking appropriate non-anomalous clustering resolution 
clustree(pancreas.integrated, prefix = "integrated_snn_res.")

# Based of clustree assessment choose res = 0.4, stable yet biologically relevant
# Beyond 0.4 massive cluster destabilization occurs
pancreas.integrated <- FindClusters(pancreas.integrated, resolution = 0.4)

# Alternatively build a cluster tree
DefaultAssay(object = pancreas.integrated) <- "integrated"
pancreas.integrated=BuildClusterTree(pancreas.integrated, slot = "scale.data")
PlotClusterTree(pancreas.integrated)

# Intrestingly in cluster tree, the euclidean distance bewteen 11 and 0 is very small. So they must be very similar cells
# We now force clusters 11 and 0 together
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
new.cluster.ids <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "0", "12", "13")
pancreas.integrated@active.ident <- plyr::mapvalues(x = pancreas.integrated@active.ident, from = current.cluster.ids, to = new.cluster.ids)

# Re-run cluster tree
DefaultAssay(object = pancreas.integrated) <- "integrated"
pancreas.integrated=BuildClusterTree(pancreas.integrated, slot = "scale.data")
PlotClusterTree(pancreas.integrated)

#################### #
#################### #

# Step 10: non-linear dimensionality assessment ####
# Run UMAP calculations
pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:10)

# Change default assay to RNA for variable gene calculation, remember to save in RNA slot
# This is necessary to look at normalised gene expression from the RNA slot, looking at raw counts is incorrect.
# In our paper, we concentrate on SCT normalized data, or log normalized RNA NOT raw counts.
# In addition if you want, you may look at SCT or integrated slots, to look at gene expression
DefaultAssay(object = pancreas.integrated)
DefaultAssay(object = pancreas.integrated) <- "SCT"
DefaultAssay(object = pancreas.integrated) <- "RNA"
DefaultAssay(object = pancreas.integrated) <- "integrated"

# Find variable features and save them in the RNA assay slot
#Normalize and store RNA data:
pancreas.integrated <- NormalizeData(pancreas.integrated, assay = "RNA")

# Say hello to my little friend UMAP, Dont try to be too smart and give the group.by command, the plot works just fine as it is.
# UMAP by cluster
DimPlot(pancreas.integrated, reduction = "umap", 
        cols =  c("darkgreen",
                  "red",
                  "sienna3",
                  "mediumseagreen",
                  "turquoise4",
                  "black",
                  "royalblue1",
                  "yellow4",
                  "gray30",
                  "darkred",
                  "orange2",
                  "darkmagenta",
                  "deeppink2"
                  ),
        pt.size = 1)

# UMAP by sample ID
DimPlot(object = pancreas.integrated, reduction = "umap", group.by = "sample", pt.size = 1)

# UMAP by experimental ID
DimPlot(object = pancreas.integrated, reduction = "umap", group.by = "exp", pt.size = 1)

FeaturePlot(object = pancreas.integrated, 
            features = c("KRTAP2-3"),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            max.cutoff = 3,
            slot = 'data',
            order = TRUE)

# Visualization in Violin plots, we change data to RNA, as scalled data ranges from -ve to +ve values.
DefaultAssay(pancreas.integrated) <- "RNA"
VlnPlot(object = pancreas.integrated,
        split.by = "exp",
        features = c("BMPR1A"),  
        c("black",
          "red"
          ),
        pt.size = 0
        )

# Looking at acinar and ductal genes
feature.plot.integrated <- c("CPA2", "SYCN",
                             "OLFM4", "CEL",
                             "XIST", "WSB1",
                             "SPP1", "CRP",
                             "IGFBP3", "TFF1",
                             "AKAP12", "KRTAP2-3"
                             )

DotPlot(object = pancreas.integrated, 
        features = feature.plot.integrated, 
        dot.scale = 10, 
        col.min = 0,
        col.max = 2,
) + scale_colour_gradient2(low = "#0200ad", mid = "#ffe272", high = "#ff0000")

# DIFFERENTIAL GENE EXPRESSION ####
# find markers for every cluster compared to all remaining cells, report only the positive ones
# remember to run this analysis against default RNA data
DefaultAssay(pancreas.integrated) <- "integrated"
integration.combined.markers <- FindAllMarkers(object = pancreas.integrated, 
                                                  features = VariableFeatures(pancreas.integrated, assay = 'integrated'), 
                                                  only.pos = TRUE, 
                                                  min.pct = 0.1, 
                                                  logfc.threshold = 0.41, 
                                                  assay = 'SCT',
                                                  slot = c('data'))

integration.combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(integration.combined.markers, 'C:/Users/mxq52/Box/NPOD PAPER/Round 4/R_files/integration.combined.markers.csv')
getwd()

# Heatmap visualization of top 50 most DE genes
top50.pancreas.integrated.genes <- integration.combined.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC) 
DoHeatmap(object = pancreas.integrated, 
          features = top50.pancreas.integrated.genes$gene, 
          disp.min = -1,
          group.colors = 
            c("darkgreen",
              "red",
              "sienna3",
              "mediumseagreen",
              "turquoise4",
              "black",
              "royalblue1",
              "yellow4",
              "gray30",
              "darkred",
              "orange2",
              "darkmagenta",
              "deeppink2"),
          group.bar.height = 0.05,
          draw.lines = TRUE,
          lines.width = 25,
          disp.max = 1) + scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", 
                                                                           "#fbfcbd", 
                                                                           "#ff0000"))(256))
################## #
################## #
# CODE END ####
