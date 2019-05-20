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

# OBJECT SETUP AND NORMALIZATION####
# Load the Human Pancreas Single-cell RNAseq Atlas.
pancreas.data <- readRDS(file = "C:/Users/mxq52/Box/FAHD SHARED JDB-RLP/NEW R files/pancreas_v3_files/pancreas_expression_matrix.rds")
metadata <- readRDS(file = "C:/Users/mxq52/Box/FAHD SHARED JDB-RLP/NEW R files/pancreas_v3_files/pancreas_metadata.rds")

# STEP 1: Create a integrated Seurat Object of selected datasets
pancreas <- CreateSeuratObject(counts = pancreas.data, meta.data = metadata)
pancreas.list <- SplitObject(object = pancreas, split.by = "tech")

# Extract and Organise Celseq dataset
celseq <- pancreas.list[["celseq"]]

# Sample specific Metadata addition
celseq$sample <- "celseq" 

# Extract and Organise Celseq2 dataset
celseq2 <- pancreas.list[["celseq2"]]

# Sample specific Metadata addition
celseq2$sample <- "celseq2" 

# Merge refrence datasets
muraro.dat <- merge(x = celseq, 
                    y = celseq2
)

#Experiment specific Metadata addition
muraro.dat$exp <- "Muraro"
table(muraro.dat$orig.ident)

# Visualize QC metrics as a violin plot
VlnPlot(object = muraro.dat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# STEP 2: THRESHOLDING: RNA and Mitochondria based cell thresholding
muraro.dat <- subset(x = muraro.dat, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)

# STEP 3: NORMALIZATION: normalizing according to S3.0.0 this also normalize against 
# i) TOTAL gene expression
# ii) x10,000 all Data
# iii) Log transforms data
muraro.dat <- NormalizeData(object = muraro.dat)

# STEP 4: VARIABLE GENES: find highly variable genes
muraro.dat <- FindVariableFeatures(object = muraro.dat, selection.method = "vst", nfeatures = 2000)

# Visualize QC results
VlnPlot(object = muraro.dat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# Pre-processing alk3n3 Dataset
# Data Preprocessing for alk3n3 dataset
# STEP 1: Create a integrated Seurat Obj of the alk3n3 datasets

# Create Seurat obj: HPD1
HPD1.data <- Read10X(data.dir = "C:/Users/mxq52/Box/Fahd Shared JUAN/Dominguez-Bendala, Juan scRNAseq prj/10x Gen input files/HPD1/filtered_feature_bc_matrix")
HPD1 <- CreateSeuratObject(counts = HPD1.data, 
                             project = "HPD1"
)

# Sample specific Metadata addition
HPD1$sample <- "HPD1"

# Create Seurat obj: HPD3
HPD2.data <- Read10X(data.dir = "C:/Users/mxq52/Box/Fahd Shared JUAN/Dominguez-Bendala, Juan scRNAseq prj/10x Gen input files/HPD2/filtered_feature_bc_matrix")
HPD2 <- CreateSeuratObject(counts = HPD2.data,
                             project = "HPD2"
)

# Sample specific Metadata addition
HPD2$sample <- "HPD2"

# Create Seurat obj: HPD3
HPD3.data <- Read10X(data.dir = "C:/Users/mxq52/Box/Fahd Shared JUAN/Dominguez-Bendala, Juan scRNAseq prj/10x Gen input files/HPD3/filtered_feature_bc_matrix")
HPD3 <- CreateSeuratObject(counts = HPD3.data,
                             project = "HPD3"
)

# Sample specific Metadata addition
HPD3$sample <- "HPD3"

# Merge refrence-qadir datasets
alk3n3.combined <- merge(x = HPD1, y = c(HPD2, HPD3), add.cell.ids = c("HPD1", "HPD2", "HPD3"), project = "ALK3n3")
alk3n3.combined

#Experiment specific Metadata addition
alk3n3.combined$exp <- "Qadir" 

# The operator can add columns to object metadata. This is a great place to stash QC stats
alk3n3.combined[["percent.mt"]] <- PercentageFeatureSet(object = alk3n3.combined, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(object = alk3n3.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# RNA based cell thresholding
alk3n3.combined <- subset(x = alk3n3.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 5)

# Data normalization
alk3n3.combined <- NormalizeData(object = alk3n3.combined)

#highly variable genes
alk3n3.combined <- FindVariableFeatures(object = alk3n3.combined, selection.method = "vst", nfeatures = 2000)

# Visualize QC results
VlnPlot(object = alk3n3.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# DATA INTEGRATION ####
# Merging the experiment with the refrence dataset
integration <- merge(x = muraro.dat, 
                     y = alk3n3.combined
)

# Perfrom integration (leave this at 20)
integration.anchors <- FindIntegrationAnchors(object.list = list(muraro.dat, alk3n3.combined), dims = 1:20)
integration.combined <- IntegrateData(anchorset = integration.anchors, dims = 1:20)

# Perform an Integrated Analysis
DefaultAssay(object = integration.combined) <- "integrated"

################## #
################## #
# CODE END ####
