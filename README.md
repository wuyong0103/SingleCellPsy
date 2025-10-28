# Summary
In this study, we collected and meta-analyzed ~3M nuclei transcriptomics from 22 published studies. Using canonical cell markers, we annotated all cell types of BICCN. We also mapped the cells to 10 age stages from first trimester to elder (60+ years). Integrating these snRNA-seq data with statistics of genome-wide association study (GWAS) and whole-exome sequencing (WES), we tried to identify cell types associated with schizophrenia (SCZ) and bipolar disorder (BD) in different age stages. <br>

## Prepare the snRNA-seq data
We downloaded the snRNA-seq data from 22 published studies. For more details, please refer to our Manuscript and Supplementat Table 1.
1. Create Seurat object for each study using scripts in the `CreateSeuratObject` directory.
2. Integrate all the Seurat objects using script `01_Integrate.R`
3. Reclassify excitatory neurons using script `02_ExNeu.R`
4. Reclassify inhibitory neurons using script `03_InNeu.R`

## Downsample the snRNA-seq data
Since the statistical power of these analyses depends on the accuracy of gene expression estimates, and cell types with fewer cells typically have less reliable average expression values, we applied a downsampling strategy to ensure comparable cell numbers across types.
1. Downsample the snRNA-seq data using script `04_Downsample.R`

## Calculate cellular specificity score
Expression specificity for each gene was calculated as the ratio of its expression in a given cell type to the sum of its expression across all cell types, yielding a score between 0 and 1, where 1 indicates complete specificity and 0 indicates no expression in that cell type.
1. Calculate expression specificity using script `04_Downsample.R` (Integrated in downsample script).

## Create EWCE object
We employed the Expression Weighted Cell Type Enrichment ([EWCE](https://github.com/NathanSkene/EWCE)) R packages to conduct the enrichment of rare variants in cell types. EWCE detects the association between trait and cell type by evaluating whether the expression of a set of genes associated with trait in a particular cell type were higher than that of randomly selected genes.
1. Create EWCE object using script `05_EwceCtd.R`
   
## Prepare the GWAS data and WES data
1. GWAS data of SCZ and BD was downloaded from [Psychiatric Genomic Consortium](https://pgc.unc.edu).
2. WES data of SCZ was downloaded from [Schizophrenia Exome Sequencing Meta-analysis (SCHEMA) consortium](https://schema.broadinstitute.org/).
3. WES data of BD was downloaded from [Bipolar Exome (BipEx) sequencing project](https://bipex.broadinstitute.org/).
4. Preprocess GWAS data using [LDSC](https://github.com/bulik/ldsc) and [MAGMA](https://cncr.nl/research/magma/).

### LDSC and MAGMA analysis
We conducted the LDSC and MAGMA analysis according to the publication of [Bryois et al., 2020, Nat Genet.](https://github.com/jbryois/scRNA_disease/tree/master).
1. Please refer to `00_Process.sh` for more running details.

## EWCE analysis
1. Conduct EWCE analysis using script `EWCE_downsample500.R`

## SCENIC analysis
We used [SCENIC](https://github.com/aertslab/SCENIC) to identify age stage-specific transcription factors and regulons.
1. Conduct EWCE analysis using script `06_SCENIC.sh`

## ClusterGVis
We identified differentially expressed genes using `FindAllMarkers` funciton of [Seurat](https://satijalab.org/seurat/) and visualize the result using [ClusterGVis](https://github.com/junjunlab/ClusterGVis)
1. Conduct differential expression analysis and visualize the result using script `09_DEAndClusterGVis.R`.
2. ClusterGVis uses [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) to conduct functional enrichment analysis.

## Visualization
1. Visualize the results of LDSC and MAGMA using script `07_Figure-LDSC-MAGMA.R`.
2. isualize the results of LEWCE using script `08_Figure-EWCE.R`.
