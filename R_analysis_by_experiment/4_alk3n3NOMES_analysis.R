# DESCRIPTION ####
# The following code ASSUMES all dependencies in R  have been installed and previous R code have been run successfully see file(s): 
# 1_environment_setup.R 
# 2_alk3n3_seurat_setup.R
# 3_alk3n3_analysis.R
# The purpose of this code is to analyze n=3 data which has been subsetted to remove mesenchymal cells, as a Seurat object.

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
library(plotly)

# CONFIRM CORRECT INSTALL ####
# Confirm package version of Seurat and Monocle
packageVersion("Seurat")
packageVersion("monocle")
packageVersion("clustree")

# SUBSET DATA ####
# In order to remove mesenchymal populations from the analysis
alk3n3.integrated.nomes <- subset(alk3n3.integrated, idents = c("0", "1", "2", "3", "4", "7"))

# looking at data after subsetting
FeatureScatter(object = alk3n3.integrated.nomes, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(object = alk3n3.integrated.nomes, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(object = alk3n3.integrated.nomes, feature1 = "percent.mt", feature2 = "nFeature_RNA", cols = c("black",
                                                                                                              "black",
                                                                                                              "black",
                                                                                                              "black",
                                                                                                              "black",
                                                                                                              "black",
                                                                                                              "black"))

# CHECKING DATA ####
# Look at your default assay
DefaultAssay(object = alk3n3.integrated.nomes)

# Change default assay to RNA temporarily, so that we can run vst and normalize RNA counts data
DefaultAssay(object = alk3n3.integrated.nomes) <- "RNA"

# Find variable features and save them in the RNA assay slot
#Normalize and store RNA data:
alk3n3.integrated.nomes <- NormalizeData(alk3n3.integrated.nomes)

# Now over-write the SCT assay with new-analyzed data from this subsetted data
alk3n3.integrated.nomes <- SCTransform(alk3n3.integrated.nomes, assay = "RNA", new.assay.name = "SCT", verbose = TRUE, return.only.var.genes = TRUE)

# Look at your default assay (usually after SCTransform, teh defauls switches to SCT)
DefaultAssay(object = alk3n3.integrated.nomes)

# Change default assay to SCT, if not already
DefaultAssay(object = alk3n3.integrated.nomes) <- "SCT"

# Identify the 10 most highly variable genes, within the SCT slot
top10.nomes.vst <- head(x = VariableFeatures(object = alk3n3.integrated.nomes), 10)

# Ploting variable features with and without labels, note how you get geometric mean of expression and residual variance
p21 <- VariableFeaturePlot(object = alk3n3.integrated.nomes, assay = 'SCT', selection.method = c('sct'))
p21
LabelPoints(plot = p21, points = top10.nomes.vst, repel = TRUE)

# Viewing data
head(x = colnames(x = alk3n3.integrated.nomes))
tail(x = colnames(x = alk3n3.integrated.nomes))
unique(x = sapply(X = strsplit(x = colnames(x = alk3n3.integrated.nomes), split = "_"), FUN = "[", 1))
table(alk3n3.integrated.nomes$orig.ident)

# LINEAR DIMENSIONALITY REDUCTION ####
# this Seurat object was created and normalized using code outlined in: 2_alk3n3_seurat_object_setup_05202019.R
# This Seurat object is called alk3n3.integrated.nomes
# Without running the R files in sequence YOU WILL NOT get the data output. Running in sequence is essential.
# In continuation of the previous code, note that your default assay is set at "RNA", you must change to integrated for PCA calculations
# Look at your default assay
DefaultAssay(object = alk3n3.integrated.nomes)

# Change default assay to integrated, save information in the "integrated" assay
DefaultAssay(object = alk3n3.integrated.nomes) <- "integrated"

# PCA analysis data will be stored in the "reductions' slot
alk3n3.integrated.nomes <- RunPCA(object = alk3n3.integrated.nomes, features = VariableFeatures(object = alk3n3.integrated.nomes))

# Examine and visualize PCA results a few different ways
print(x = alk3n3.integrated.nomes[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object = alk3n3.integrated.nomes, dims = 1:2, reduction = "pca")
DimPlot(object = alk3n3.integrated.nomes, reduction = "pca")
DimHeatmap(object = alk3n3.integrated.nomes, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(object = alk3n3.integrated.nomes, dims = 1:15, cells = 500, balanced = TRUE)

# Using JacStraw plots to visualize LDR measurements
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
# alk3n3.integrated.nomes <- JackStraw(object = alk3n3.integrated.nomes, num.replicate = 100)
# alk3n3.integrated.nomes <- ScoreJackStraw(object = alk3n3.integrated.nomes, dims = 1:10)
# JackStrawPlot(object = alk3n3.integrated.nomes, dims = 1:10)

# Elbow plot
ElbowPlot(object = alk3n3.integrated.nomes)

# CLUSTER ANALYSIS ####
# Clustering, we run multiple permutations to allow clustree to analyze optimal clustering resolution.
alk3n3.integrated.nomes <- FindNeighbors(object = alk3n3.integrated.nomes, dims = 1:20)
alk3n3.integrated.nomes <- FindClusters(object = alk3n3.integrated.nomes, resolution = 0)
alk3n3.integrated.nomes <- FindClusters(object = alk3n3.integrated.nomes, resolution = 0.1)
alk3n3.integrated.nomes <- FindClusters(object = alk3n3.integrated.nomes, resolution = 0.2)
alk3n3.integrated.nomes <- FindClusters(object = alk3n3.integrated.nomes, resolution = 0.3)
alk3n3.integrated.nomes <- FindClusters(object = alk3n3.integrated.nomes, resolution = 0.4)
alk3n3.integrated.nomes <- FindClusters(object = alk3n3.integrated.nomes, resolution = 0.5)
alk3n3.integrated.nomes <- FindClusters(object = alk3n3.integrated.nomes, resolution = 0.6)

# Cluster-tree analysis, looking appropriate non-anomalous clustering resolution 
clustree(alk3n3.integrated.nomes, prefix = "integrated_snn_res.")

# We run our final, chosen clustering resolution, so that seurat_clusters an be populated for further analysis
alk3n3.integrated.nomes <- FindClusters(object = alk3n3.integrated.nomes, resolution = 0.4)

# Check data
table(alk3n3.integrated.nomes$seurat_clusters)
table(Idents(alk3n3.integrated.nomes), alk3n3.integrated.nomes$orig.ident)

# NON-LINEAR DIMENSIONALITY REDUCTION ####
# RunUMAP
alk3n3.integrated.nomes <- RunUMAP(alk3n3.integrated.nomes, dims = 1:20)

# DATA-VISUALIZATION ####
# Note: you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
# UMAP Visualization
DimPlot(object = alk3n3.integrated.nomes, reduction = "umap", pt.size = 2, cols = c("darkgreen",
                                                                                    "orange2",
                                                                                    "lightseagreen",
                                                                                    "royalblue1",
                                                                                    "red4",
                                                                                    "red",
                                                                                    "darkmagenta"
                                                                                    ))

# UMAP plot where cells are colored by replicate
DimPlot(alk3n3.integrated.nomes, group.by = c("orig.ident"), combine = FALSE, pt.size = 1, cols = c("seagreen", 
                                                                                                    "deepskyblue4", 
                                                                                                    "firebrick2"))

# Define an order of cluster identities remember after this step-
# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c(2, 0, 3, 4, 5, 1, 6)

# Re-level object@ident
alk3n3.integrated.nomes@active.ident <- factor(x = alk3n3.integrated.nomes@active.ident, levels = my_levels)

# Look at your default assay
DefaultAssay(object = alk3n3.integrated.nomes)

# Change default assay to SCT all analysis is based on the SCT slot
DefaultAssay(object = alk3n3.integrated.nomes) <- "RNA"

# Violin plots for gene expression. change gene value to look at diffrent genes
VlnPlot(object = alk3n3.integrated.nomes, 
        features = c("P2RY1"), 
        ncol = 1,
        slot = 'data',
        assay = 'RNA',
        cols = c("lightseagreen",
                 "darkgreen",
                 "royalblue",
                 "red4",
                 "red",
                 "orange2",
                 "darkmagenta"
                 ), y.max = 3, pt.size = 1)

# 2D UMAP plots for gene expression. change gene value to look at diffrent genes
FeaturePlot(object = alk3n3.integrated.nomes, 
            features = c("P2RY1"),
            pt.size = 0.5,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            max.cutoff = 0.5,
            order = TRUE)


# DIFFERENTIALLY EXPRESSED-GENE ANALYSIS ####
# Change default assay to integrated so that de analysis can be performed on all genes
# Look at your default assay
DefaultAssay(object = alk3n3.integrated.nomes)

# Change default assay to integrated, save information in the "integrated" assay
DefaultAssay(object = alk3n3.integrated.nomes) <- "SCT"

# Find markers for every cluster compared to all remaining cells, report only the positive ones
# Here we define a DE gene as a gene which has:
# Fold Change of >1.5x
# Atleast 10% of cells express that gene
# p value < 0.001
# Our analysis is based of only variable genes found in the SCT assay slot.
alk3n3.integrated.nomes.markers <- FindAllMarkers(object = alk3n3.integrated.nomes, 
                                                  features = VariableFeatures(alk3n3.integrated.nomes, assay = 'SCT'), 
                                                  only.pos = TRUE, 
                                                  min.pct = 0.1, 
                                                  logfc.threshold = 0.41, 
                                                  assay = 'SCT',
                                                  slot = c('data'))

alk3n3.integrated.nomes.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(alk3n3.integrated.nomes.markers, 'C:/Users/mxq52/Box/NPOD PAPER/Round 4/R_files/alk3n3.integrated.nomes.markers.csv')

# Look at your default assay
DefaultAssay(object = alk3n3.integrated.nomes)

# Change default assay to SCT, save information in the "SCT" assay
# You can toggle between integrated, SCT and RNA to see different expression profiles/different normalizations
DefaultAssay(object = alk3n3.integrated.nomes) <- "SCT"

# Create heatmap using doheatmap
top30.nomes <- alk3n3.integrated.nomes.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
DoHeatmap(object = alk3n3.integrated.nomes, 
          features = top30.nomes$gene, 
          disp.min = -2.5, 
          disp.max = 2.5) + scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", 
                                                                             "#fbfcbd", 
                                                                             "#ff0000"))(256))

# Dot plots - the size of the dot corresponds to the percentage of cells
# expressing the gene in each cluster. The color represents the average
# expression level
# Look at your default assay
DefaultAssay(object = alk3n3.integrated.nomes)

# Change default assay to RNA, save information in the "RNA" assay
# You can toggle between integrated, SCT and RNA to see different expression profiles/different normalizations
DefaultAssay(object = alk3n3.integrated.nomes) <- "RNA"

#Changing gene combinaion changes genes on the dotplot
feature.plot <- c("KRTAP2-3", "AKAP12",
                  "TFF1", "IGFBP3",
                  "CRP", "SPP1",
                  "WSB1", "XIST",
                  "CEL", "OLFM4",
                  "SYCN", "CPA2",
                  "SRGN", "IDO1")

DotPlot(object = alk3n3.integrated.nomes, 
        features = feature.plot, 
        dot.scale = 10, 
        col.min = 0,
        col.max = 2,
        ) + scale_colour_gradient2(low = "#0200ad", mid = "#ffe272", high = "#ff0000")

# Looking at acinar and ductal genes
feature.plot.ad <- c("SPP1", "CFTR",
                  "AQP1", "ALDH1A3",
                  "KRT19", "CRP",
                  "DEFB1", "CEACAM6",
                  "MMP7", "TSPAN8",
                  "ONECUT2", "LITAF",
                  "SOX4", "DAB2",
                  "CREB5", "HLA-DQB1",
                  "WWTR1", "PPARGC1A",
                  "PKHD1", "NFIB",
                  "PNLIP", "REG1B",
                  "PRSS1", "ALB",
                  "CPA2", "CTRB2",
                  "CEL", "PLA2G1B",
                  "CELA3A", "GATA4",
                  "MECOM", "NR5A2",
                  "ZFP36L1", "CEBPD",
                  "CREB3L1", "XBP1",
                  "LGR4", "NUPR1")

DotPlot(object = alk3n3.integrated.nomes, 
        features = feature.plot.ad, 
        dot.scale = 10, 
        col.min = 0,
        col.max = 2,
) + scale_colour_gradient2(low = "#0200ad", mid = "#ffe272", high = "#ff0000")

# Visualize co-expression of two features simultaneously
FeaturePlot(alk3n3.integrated.nomes, 
            features = c("PDX1", "P2RY1"), 
            blend = TRUE, 
            pt.size = 1, 
            order = TRUE, 
            blend.threshold = 0.1, 
            label.size = 6, 
            sort.cell = TRUE, 
            )

# Save RDS file
# saveRDS(alk3n3.integrated, file = "C:/Users/mxq52/Box/NPOD PAPER/Round 4/R_files/alk3n3.integrated.code3.rds")

################## #
################## #
# CODE END ####
