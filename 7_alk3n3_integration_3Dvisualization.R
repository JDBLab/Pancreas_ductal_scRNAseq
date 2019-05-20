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

# INTERACTIVE 3D tSNE PLOTING ####
# Plotting all cells in 3D, having the same color combinations for each celll as in the 2D tSNE analysis
plotting.data <- FetchData(object = yourseuratobject, vars = c("tSNE_1", "tSNE_2", "tSNE_3", "seurat_clusters"))
plotting.data$label <- paste(rownames(plotting.data))
plot_ly(data = plotting.data, 
        x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
        color = ~seurat_clusters, 
        opacity = .5,
        colors = c("lightseagreen",
                   "gray50",
                   "darkgreen",
                   "red4",
                   "red",
                   "turquoise4",
                   "black",
                   "yellow4",
                   "royalblue1",
                   "lightcyan3",
                   "peachpuff3",
                   "khaki3",
                   "gray20",
                   "orange2",
                   "royalblue4",
                   "yellow3",
                   "gray80",
                   "darkorchid1",
                   "lawngreen",
                   "plum2",
                   "darkmagenta"),
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 4, width=2), 
        text=~label,
        hoverinfo="text")

# Looking at gene expression for a particular gene and plotting that in 3D
plotting.data <- FetchData(object = yourseuratobject, vars = c("tSNE_1", "tSNE_2", "tSNE_3", "CTNND1"))
plotting.data$changed <- ifelse(test = plotting.data$CTNND1 <1, yes = plotting.data$CTNND1, no = 1)
plotting.data$label <- paste(rownames(plotting.data)," - ", plotting.data$CTNND1, sep="")

plot_ly(data = plotting.data, 
        x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
        color = ~changed,
        opacity = .5,
        colors = c('darkgreen', 'red'), 
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 4, width=2), 
        text=~label,
        hoverinfo="text"
)

################## #
################## #
# CODE END ####
