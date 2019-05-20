# DESCRIPTION ####
# The following code ASSUMES all dependencies in R  have been installed and previous R code have been run successfully see file(s): 
# 1_environment_setup.R 
# 5_alk3n3_integration_with_cellatlas.R
# The purpose of this code is integrate the alk3n3 dataset with a human pancreas single cell atlas. This data has been utlized from 
# the following datasets: GSE81076 and GSE85241. A gene expression matrix of this dataset exist, as part of a Seurat v3.0.0 tutorial
# those files can be found here https://satijalab.org/seurat/v3.0/pancreas_integration_label_transfer.html
# These files are necessary for the analysis. We will attempt to upload them to GEO for asy access.

# R v3.5.3 (x64 bit) and RStudio v1.2.1335 (x64 bit) were used for running this code (see readme on how to install)
# Seurat is a multimodal single Cell RNA seq analysis algorithm created by
# The Satija Lab. For more information please see: https://satijalab.org/seurat/

# DATA VISUALIZATION ####
# Run the standard workflow for visualization and clustering
# It is important to note here, that nultiple permutations were run of optimal clustering, those clusters were chosen which 
# contextually made sense based on statistical tests, DE-genes/cluster amd biological relevance
integration.combined <- ScaleData(object = integration.combined, verbose = FALSE)
integration.combined <- RunPCA(object = integration.combined, npcs = 30, verbose = FALSE)
table(integration.combined$orig.ident)

# t-SNE and Clustering (do 30/0.8 for analysis)
integration.combined <- RunTSNE(object = integration.combined, reduction = "pca", dims = 1:30)
integration.combined <- FindNeighbors(object = integration.combined, reduction = "pca", dims = 1:30)
integration.combined <- FindClusters(integration.combined, resolution = 0.8)

# Visualization in tSNE plots
DimPlot(object = integration.combined, reduction = "tsne", group.by = "sample", pt.size = 2)
DimPlot(object = integration.combined, reduction = "tsne", group.by = "exp", pt.size = 2)
DimPlot(object = integration.combined, 
        reduction = "tsne", 
        label = FALSE,
        group.by = "seurat_clusters",
        cols =  c("lightseagreen",
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
                  "darkmagenta"
        ),
        pt.size = 2)


# Visualization in tSNE plots
FeaturePlot(object = integration.combined, features = c("MIR7-1"), min.cutoff =0, max.cutoff = 2, 
            cols = c("grey", "red"), pt.size = 2)


# Visualization in Violin plots, we change data to RNA, as scalled data ranges from -ve to +ve values.
DefaultAssay(integration.combined) <- "RNA"
VlnPlot(object = integration.combined, features = c("PPY"),  cols =  c("lightseagreen",
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
                                                                       "darkmagenta"
))

# Return data to now be based of integrated analysis
DefaultAssay(integration.combined) <- "integrated"

# DIFFERENTIAL GENE EXPRESSION ####
# find markers for every cluster compared to all remaining cells, report only the positive ones
# remember to run this analysis against default RNA data
DefaultAssay(integration.combined) <- "RNA"
integration.combined.markers <- FindAllMarkers(object = integration.combined, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.7)
integration.combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(integration.combined.markers, 'integration_combined.csv')
getwd()

# I have re-analyzed for var and scaled as some genes in the heatmap, may not be variable genes.
DefaultAssay(integration.combined) <- "RNA"
integration.combined <- FindVariableFeatures(object = integration.combined, selection.method = "vst", nfeatures = 30000)
integration.combined <- ScaleData(object = integration.combined, verbose = FALSE)

# Heatmap visualization
top20 <- integration.combined.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC) 
DoHeatmap(object = integration.combined, 
          features = top20$gene, 
          disp.min = -2.5, 
          disp.max = 2.5) + scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", 
                                                                             "#fbfcbd", 
                                                                             "#ff0000"))(256))

################## #
################## #
# CODE END ####
