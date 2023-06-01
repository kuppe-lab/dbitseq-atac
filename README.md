# Spatial ATAC-seq
A snakemake pipeline to process spatial ATAC-seq raw data

Modified from https://github.com/dyxmvp/Spatial_ATAC-seq/blob/main/Data_preprocessing/README.md

## Dependiencies

* [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html). 
* [Biopython](https://biopython.org/docs/1.75/api/index.html).
* [Cell Ranger ATAC](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/installation)
* [BBMap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/).

## Run the pipeline
1. Replace the cellranger-atac-2.1.0/lib/python/atac/barcodes/737K-cratac-v1.txt.gz with the new barcode file in this folder.
2. Configure Snakefile
3. Configure cluster.json
4. Configure Snakemake.sh
5. Place paired FASTQ's in input directory. Folder names should contain "Mouse" or "Human" for correct genome alignment.
6. To run the pipeline, use the command:
```
sbatch Snakemake.sh
```
