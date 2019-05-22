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
library(plot_ly)

# CONFIRM CORRECT INSTALL ####
# Confirm package version of Seurat and Monocle
packageVersion("Seurat")
packageVersion("monocle")

# SUBSET DATA ####
# In order to remove mesenchymal populations from the analysis
alk3n3.combined.nomes <- SubsetData(object = alk3n3.combined, ident.remove = c("6", "7"))

# looking at data after subsetting
FeatureScatter(object = alk3n3.combined.nomes, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(object = alk3n3.combined.nomes, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(object = alk3n3.combined.nomes, feature1 = "percent.mt", feature2 = "nFeature_RNA")

# CHECKING DATA ####
# Thresholding
head(x = colnames(x = alk3n3.combined.nomes))
tail(x = colnames(x = alk3n3.combined.nomes))
unique(x = sapply(X = strsplit(x = colnames(x = alk3n3.combined.nomes), split = "_"), FUN = "[", 1))
table(alk3n3.combined.nomes$orig.ident)


# LINEAR DIMENSIONALITY REDUCTION ####
# Calculating Highly variable genes
alk3n3.combined.nomes <- FindVariableFeatures(object = alk3n3.combined.nomes, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = alk3n3.combined.nomes), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = alk3n3.combined.nomes)
plot1
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
alk3n3.combined.nomes <- JackStraw(object = alk3n3.combined.nomes, num.replicate = 100)
alk3n3.combined.nomes <- ScoreJackStraw(object = alk3n3.combined.nomes, dims = 1:20)
JackStrawPlot(object = alk3n3.combined.nomes, dims = 1:15)

# Re-runing clustering analysis
alk3n3.combined.nomes <- FindNeighbors(object = alk3n3.combined.nomes, dims = 1:10, force.recalc = TRUE)
alk3n3.combined.nomes <- FindClusters(object = alk3n3.combined.nomes, resolution = 0.5, force.recalc = TRUE)

# Non-linear dimensionality reduction using tSNE
alk3n3.combined.nomes <- RunTSNE(object = alk3n3.combined.nomes, dims = 1:10)

# DATA-VISUALIZATION ####
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(object = alk3n3.combined.nomes, reduction = "tsne", pt.size = 2, cols = c("darkgreen",
                                                                                  "turquoise4",
                                                                                  "lightseagreen",
                                                                                  "red4",
                                                                                  "red",
                                                                                  "royalblue1",
                                                                                  "orange2",
                                                                                  "yellow3",
                                                                                  "royalblue4",
                                                                                  "darkorchid1",
                                                                                  "lawngreen",
                                                                                  "plum2",
                                                                                  "darkmagenta"))

# tSNE plot where cells are colored by replicate:  First, store the current
# identities in a new column of meta.data called CellType
alk3n3.combined.nomes$CellType <- Idents(object = alk3n3.combined.nomes)

# Next, switch the identity class of all cells to reflect replicate ID
Idents(object = alk3n3.combined.nomes) <- "orig.ident"
DimPlot(object = alk3n3.combined.nomes, reduction = "tsne", cols = c("seagreen",
                                                                     "deepskyblue4",
                                                                     "firebrick2"), pt.size = 2 )

table(Idents(object = alk3n3.combined.nomes), alk3n3.combined.nomes$seurat_clusters)

# Switch back to cell type labels
Idents(object = alk3n3.combined.nomes) <- "CellType"

# Define an order of cluster identities remember after this step-
# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c(7, 6, 4, 3, 1, 2, 5, 10, 8, 0, 9, 11, 12)

# Re-level object@ident
alk3n3.combined.nomes@active.ident <- factor(x = alk3n3.combined.nomes@active.ident, levels = my_levels)

# Violin plots for gene expression. change gene value to look at diffrent genes
VlnPlot(object = alk3n3.combined.nomes, features = c("CTNND1"), cols = c("yellow3",
                                                                         "orange2",
                                                                         "red",
                                                                         "red4",
                                                                         "turquoise4",
                                                                         "lightseagreen",
                                                                         "royalblue1",
                                                                         "lawngreen",
                                                                         "royalblue4",
                                                                         "darkgreen",
                                                                         "darkorchid",
                                                                         "plum2",
                                                                         "darkmagenta"))

# 2D tSNE plots for gene expression. change gene value to look at diffrent genes
FeaturePlot(object = alk3n3.combined.nomes, 
            features = c("CTNND1"),
            pt.size = 2,
            cols = c("darkgrey", "red"),
            min.cutoff = 0.2,
            max.cutoff = 1)

# DIFFERENTIALLY EXPRESSED-GENE ANALYSIS ####
# Find markers for every cluster compared to all remaining cells, report only the positive ones
# Here we defin a DE gene as a gene which has:
# Fold Change of >2
# Atleast 50% of cells express that gene
# p value < 0.001
# find markers for every cluster compared to all remaining cells, report only the positive ones
alk3n3.combined.markers.nomes <- FindAllMarkers(object = alk3n3.combined.nomes, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.7)
alk3n3.combined.markers.nomes %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(alk3n3.combined.markers.nomes, 'alk3n3.combined.markers.nomes.csv')

# Heatmap of all clusters, naming is from 0-12 
# Color Pallette https://htmlcolorcodes.com/
# Create a heatmap of top 10 most differentially expressed genes
top20 <- alk3n3.combined.markers.nomes %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC) 
DoHeatmap(object = alk3n3.combined.nomes, 
          features = top20$gene, 
          disp.min = -1.5, 
          disp.max = 1.5) + scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", 
                                                                             "#fbfcbd", 
                                                                             "#ff0000"))(256))

# Dot plots - the size of the dot corresponds to the percentage of cells
# expressing the gene in each cluster. The color represents the average
# expression level
#Changing gene combinaion changes genes on the dotplot
feature.plot <- c("CELA2B", "PGA5",
                  "CPA1", "CELA3B",
                  "OLFM4", "REG3A",
                  "SELENOM", "GSTM4",
                  "SPP1", "CRP",
                  "CXCL6", "CRYAB",
                  "TFPI2", "FGFBP1",
                  "INHBA", "KRT75",
                  "TRIM54", "TFF2",
                  "PGC", "RFLNA",
                  "IDO1", "CCL17",
                  "TYROBP", "FCER1G",
                  "VWF", "PECAM1"
)

DotPlot(object = alk3n3.combined.nomes, 
        features = feature.plot, 
        plot.legend = TRUE, 
        dot.scale = 8, 
        col.min = -1,
        col.max = 3,
        cex.use = 10
) + scale_colour_gradient2(low = "#0200ad", mid = "#ffe272", high = "#ff0000")


#Saving the entire file
saveRDS(alk3n3.combined, file = "~/alk3n3.nomes.combined.rds")

################## #
################## #
# CODE END ####
