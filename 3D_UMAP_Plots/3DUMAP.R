# DESCRIPTION ####
# The following code ASSUMES all dependencies in R  have been installed and previous R code have been run successfully see file(s): 
# 1_environment_setup.R 
# 5_alk3n3_integration_with_cellatlas.R
# 6_alk3n3_integration_visualization
# The purpose of this code is generate 3D tSNE visualizations of integrated the alk3n3/Cell atlas dataset
# Please visit Dragonmasetrx87 and the repository: https://github.com/Dragonmasterx87/3D-Plotting-in-Seurat-3.0.0
# For more details on how to make nice 3D plots using plot_ly

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

# CONFIRM CORRECT INSTALL ####
# Confirm package version of Seurat and Monocle
packageVersion("Seurat")
packageVersion("monocle")
packageVersion("plotly")

# 3D PLOTTING ####
# Construct a dataframe using data from your pre-clustered Seurat v3.0.0 object
# Here 'seurat_clusters' is list of numeric cluster identities, you can find it here: yourseuratobject[["seurat_cluster"]], 
# or yourseuratobject$seurat_clusters, where 'yourseuratobject' is a Seurat object created with Seurat v3.0.0
yourseuratobject <- pancreas.integrated

# Re-run UMAP s that you have accurate calculations for all tSNE(s)
yourseuratobject <- RunUMAP(yourseuratobject,
                            dims = 1:10,
                            n.components = 3L)

# Extract UMAP information from Seurat Object
umap_1 <- yourseuratobject[["umap"]]@cell.embeddings[,1]
umap_2 <- yourseuratobject[["umap"]]@cell.embeddings[,2]
umap_3 <- yourseuratobject[["umap"]]@cell.embeddings[,3]

# Visualize what headings are called so that you can extract them to form a dataframe
Embeddings(object = yourseuratobject, reduction = "umap")

# Prepare a dataframe for cell plotting
plotting.data <- FetchData(object = yourseuratobject, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "seurat_clusters"))

# Make a column of row name identities (these will be your cell/barcode names)
plotting.data$label <- paste(rownames(plotting.data))

# Plot your data, in this example my Seurat object had 21 clusters (0-20)
plot_ly(data = plotting.data, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~seurat_clusters, 
        colors = c("darkgreen",
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
                   "darkgreen",
                   "darkmagenta",
                   "deeppink2"
                   ),
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 3, width=2), # controls size of points
        text=~label, #This is that extra column we made earlier for which we will use
        hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names

# Say you wanto make a gene-expression 3D plot, where you can plot gene expression against a color scale
# Here using the same seurat object as above, we extract gene expression information for beta-actin 'ACTB'
# In our paper, we concentrate on SCT normalized data, or log normalized RNA NOT raw counts.
# In addition if you want, you may look at normalised-RNA, SCT or integrated slots, to look at gene expression
DefaultAssay(object = yourseuratobject)
DefaultAssay(object = yourseuratobject) <- "SCT"
DefaultAssay(object = yourseuratobject) <- "RNA"
DefaultAssay(object = yourseuratobject) <- "integrated"

# create a dataframe
plotting.data <- FetchData(object = yourseuratobject, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "THY1"), slot = 'data')

# Say you want change the scale, so that every cell having an expression >1 will be one color
# Basically, you are re-adjusting the scale here, so that any cell having a certain expression 
# will light up on your 3D plot

# First make another column in your dataframe, where all values above 1 are re-assigned a value of 1
# This information is stored in the 'changed' column of your dataframe
plotting.data$Expr. <- ifelse(test = plotting.data$THY1 <2, yes = plotting.data$THY1, no = 2)

# Add the label column, so that now the column has 'cellname-its expression value'
plotting.data$label <- paste(rownames(plotting.data)," - ", plotting.data$THY1, sep="")

# Plot your data, in this example my Seurat object had 21 clusters (0-20), and cells express a gene called ACTB
plot_ly(data = plotting.data, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~Expr., # you can just run this against the column for the gene as well using ~ACTB
        opacity = .5,
        colors = c('darkgreen', 'red'), 
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 3, width=2), 
        text=~label,
        hoverinfo="text"
)

################## #
################## #
# CODE END ####
