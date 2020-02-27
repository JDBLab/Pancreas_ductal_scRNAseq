# DESCRIPTION ####
# The following code ASSUMES all dependencies in R  have been installed and previous R code have been run successfully see file(s): 
# 1_environment_setup.R 
# 5_alk3n3_seurat_setup.R
# 3_alk3n3_analysis.R
# 4_alk3n3NOMES_analysis.R
# 5_alk3n3MuraroGrun.R
# 6_3DUMAP.R
# The purpose of this code is to perform monocle analysis of the alk3n3 dataset

# R v3.5.3 (x64 bit) and RStudio v1.2.1335 (x64 bit) were used for running this code (see readme on how to install)
# Seurat is a multimodal single Cell RNA seq analysis algorithm created by
# The Satija Lab. For more information please see: https://satijalab.org/seurat/
# Monocle is a multimodal single Cell RNA seq analysis algorithm created by
# The Trapnell lab. For more information please see: http://cole-trapnell-lab.github.io/monocle-release/

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

# Monocle
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("monocle")
install.packages('view')

# Load Monocle and Seurat
library(Seurat)
library(monocle)
packageVersion("Seurat")
packageVersion("monocle")

# CONFIRM CORRECT INSTALL ####
# Confirm package version of Seurat and Monocle
packageVersion("Seurat")
packageVersion("monocle")

# Output should be:
# >[1] '3.1.1'
# >[1] '2.10.1'

# SUBSET DATA ####
# In order to remove mesenchymal populations from the analysis
seurat.object <- subset(alk3n3.integrated.nomes, idents = c("0", "1", "2", "3", "4", "5"))

#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(seurat_object@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = seurat_object@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
alk3CDS <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       #lowerDetectionLimit = 0.5,
                       expressionFamily = uninormal())

#View data
pData(alk3CDS)
fData(alk3CDS)

#Run ordering algorithm
var_genes <- seurat_object[["integrated"]]@var.features
ordering_genes <- var_genes

alk3CDS <- setOrderingFilter(alk3CDS, ordering_genes)
print(dim(exprs(alk3CDS)))

## reduce dimension - do not normalize or include pseudo count. do use monocle scaling
alk3CDS <- reduceDimension(alk3CDS,
                        norm_method="none", 
                        reduction_method="DDRTree",
                        max_components=2,
                        pseudo_expr = 0,
                        relative_expr = TRUE,
                        scaling = TRUE,
                        verbose=TRUE
                        )

# First decide what you want to color your cells by
print(head(pData(alk3CDS)))

## order cells
alk3CDS <- orderCells(alk3CDS)
plot_cell_trajectory(alk3CDS, 
                     color_by = "seurat_clusters",
                     theta = -70,
                     show_branch_points = TRUE,
                     show_tree = TRUE,
                     cell_size = 3) + scale_color_manual(breaks = c("X", "Y", "Z"), 
                                                         values=c("darkgreen",
                                                                  "orange2",
                                                                  "lightseagreen",
                                                                  "royalblue1",
                                                                  "red4",
                                                                  "red"
                                                                  )) + theme(legend.position = "right")

plot_cell_trajectory(alk3CDS, color_by = "seurat_clusters", 
                     cell_size = 2, 
                     theta = -70,
                     show_tree = FALSE,
                     markers = "TFF1",
                     show_branch_points = FALSE,
                     use_color_gradient = TRUE,
                     markers_linear = FALSE) #+ scale_color_manual(breaks = c("X", "Y", "Z"), 
                                                                 values=c("darkgreen",
                                                                          "orange2",
                                                                          "lightseagreen",
                                                                          "royalblue1",
                                                                          "red4",
                                                                          "red",
                                                                          "darkmagenta") + theme(legend.position = "right") 

#Pseudotemporal lineage analysis
alk3CDS_filtered <- alk3CDS
alk3CDS_filtered_my_genes <- row.names(subset(fData(alk3CDS_filtered),
                             gene_short_name %in% c("KRT19")))
alk3CDS_filtered_cds_subset <- alk3CDS_filtered[alk3CDS_filtered_my_genes,]
plot_genes_in_pseudotime(alk3CDS_filtered_cds_subset, cell_size = 1, color_by = "seurat_clusters", y.max = 2 ) + scale_color_manual(breaks = c("X", "Y", "Z"), 
                                                                                                       values=c("darkgreen",
                                                                                                                "orange2",
                                                                                                                "lightseagreen",
                                                                                                                "royalblue1",
                                                                                                                "red4",
                                                                                                                "red",
                                                                                                                "darkmagenta"))

################## #
################## #
# CODE END ####
