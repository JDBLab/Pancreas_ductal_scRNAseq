Code to process FASTQ files to 10X filtered barcodes used in further processing

/cellranger-3.0.2/cellranger count --localcores=16 --id=SampleID --transcriptome=/refdata-cellranger-GRCh38-3.0.0/ --fastqs=/fastq_dir/ --sample=SampleID --expect-cells=1000
