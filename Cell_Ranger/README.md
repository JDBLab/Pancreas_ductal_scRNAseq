# What is this?
In the attached file you will fnd the length of code used for analyzing fastq files using Cellranger v3.0.2.
For more on tutorials and using Cellranger for scRNAseq analysis visit the primary readme file for details.
Code to process FASTQ files to 10X filtered barcodes used in further processing

/cellranger-3.0.2/cellranger count --localcores=16 --id=SampleID --transcriptome=/refdata-cellranger-GRCh38-3.0.0/ --fastqs=/fastq_dir/ --sample=SampleID --expect-cells=1000
