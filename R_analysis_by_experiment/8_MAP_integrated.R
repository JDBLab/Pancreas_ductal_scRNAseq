# DESCRIPTION ####
# The following code ASSUMES all dependencies in R  have been installed and previous R code have been run successfully see file(s): 
# 1_environment_setup.R 
# 5_alk3n3_seurat_setup.R
# 3_alk3n3_analysis.R
# 4_alk3n3NOMES_analysis.R
# 5_alk3n3MuraroGrun.R
# 6_3DUMAP.R
# 7_MAP_integrated.R
# The purpose of this code is to perform monocle analysis of the alk3n3 dataset

# R v3.5.3 (x64 bit) and RStudio v1.2.1335 (x64 bit) were used for running this code (see readme on how to install)
# Seurat is a multimodal single Cell RNA seq analysis algorithm created by
# The Satija Lab. For more information please see: https://satijalab.org/seurat/
# Monocle is a multimodal single Cell RNA seq analysis algorithm created by
# The Trapnell lab. For more information please see: http://cole-trapnell-lab.github.io/monocle-release/

### monocle for timecourse data
### use seurat normalized and transformed data as input into monocle generated in prepare_for_monocle_v1.R

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

# LOAD DATA ####
#Load Seurat object
seurat_object_integrated <- pancreas.integrated

#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(seurat_object_integrated@assays$integrated@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = seurat_object_integrated@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
seurat_object_cds <- newCellDataSet(data,
                            phenoData = pd,
                            featureData = fd,
                            #lowerDetectionLimit = 0.5,
                            expressionFamily = uninormal())

#View data
pData(seurat_object_cds)
fData(seurat_object_cds)


#Run ordering algorithm
var_genes <- seurat_object_integrated[["integrated"]]@var.features
ordering_genes <- var_genes
seurat_object_cds <- setOrderingFilter(seurat_object_cds, ordering_genes)
print(dim(exprs(seurat_object_cds)))

## reduce dimension - do not normalize or include pseudo count. do use monocle scaling
seurat_object_cds <- reduceDimension(seurat_object_cds,
                            norm_method="none", 
                            reduction_method="DDRTree",
                            max_components=3,
                            pseudo_expr = 0,
                            #relative_expr = TRUE,
                            scaling = TRUE,
                            verbose=TRUE
                            )

# First decide what you want to color your cells by
print(head(pData(seurat_object_cds)))

## order cells
seurat_object_cds <- orderCells(seurat_object_cds)
plot_cell_trajectory(seurat_object_cds, 
                     color_by = "seurat_clusters",
                     theta = -10,
                     show_branch_points = FALSE,
                     show_tree = TRUE,
                     cell_size = 3) + scale_color_manual(breaks = c("X", "Y", "Z"), 
                                                         values=c("darkgreen",
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
                                                                  "deeppink2")) + theme(legend.position = "right")

plot_cell_trajectory(seurat_object_cds, color_by = "seurat_clusters", 
                     cell_size = 0, 
                     theta = -10,
                     show_tree = FALSE,
                     markers = "KRT19",
                     show_branch_points = FALSE,
                     markers_linear = TRUE) + scale_color_manual(breaks = c("X", "Y", "Z"), 
                                                                 values=c("darkgreen",
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
                                                                          "deeppink2")) + theme(legend.position = "right")


#Pseudotemporal lineage analysis
pancreas.integrated_filtered <- seurat_object_cds
my_genes <- row.names(subset(fData(pancreas.integrated_filtered),
                             gene_short_name %in% c("KRT19")))
cds_subset_integrated <- pancreas.integrated_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset_integrated, cell_size = 2, color_by = "seurat_clusters") + scale_color_manual(breaks = c("X", "Y", "Z"), 
                                                                                                       values=c("darkgreen",
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
                                                                                                                "deeppink2"))


# Pseudotemproal Heatmap select those genes showing nice seperation

genex <- c("MUCL3",
           "CRP",
           "CRYAB",
           "GSTM4",
           "REG3A",
           "TFPI2",
           "CPA1",
           "PGA5",
           "TRIM54",
           "INHBA"
)

sig_gene_names <- (genex)
head(sig_gene_names)
pseudotemporalplot <- plot_pseudotime_heatmap(pancreas.integrated[sig_gene_names],
                                              num_clusters = 9, 
                                              cores = 4,
                                              hmcols = NULL,
                                              show_rownames = T)


pseudotemporalplot
