# scDEA
R package applies ensemble learning for single-cell differential expression analysis on single-cell RNA-seq dataset.

The scDEA package has the following R-pakage dependencies: BPSC, DEsingle, DESeq2, edgeR, MAST, monocle, scDD, limma, Seurat, zingeR, SingleCellExperiment, scater, aggregation. In the dependencies, you can load on most of dependency packages on your R. If the dependencies are not installed correctly, please install them by yourself.

## BPSC (https://github.com/nghiavtr/BPSC)
library(devtools)

install_github("nghiavtr/BPSC")

## DEsingle (http://www.bioconductor.org/packages/release/bioc/html/DEsingle.html)
if (!requireNamespace("BiocManager", quietly = TRUE))

    install.packages("BiocManager")

BiocManager::install("DEsingle")

## DESeq2 (http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html)
if (!requireNamespace("BiocManager", quietly = TRUE))

    install.packages("BiocManager")

BiocManager::install("DESeq2")

## edgeR (http://www.bioconductor.org/packages/release/bioc/html/edgeR.html)
if (!requireNamespace("BiocManager", quietly = TRUE))

    install.packages("BiocManager")

BiocManager::install("edgeR")

## MAST (http://www.bioconductor.org/packages/release/bioc/html/MAST.html)
if (!requireNamespace("BiocManager", quietly = TRUE))

    install.packages("BiocManager")

BiocManager::install("MAST")

## monocle (http://www.bioconductor.org/packages/release/bioc/html/monocle.html)
if (!requireNamespace("BiocManager", quietly = TRUE))

    install.packages("BiocManager")

BiocManager::install("monocle")

## scDD (http://www.bioconductor.org/packages/release/bioc/html/scDD.html)
if (!requireNamespace("BiocManager", quietly = TRUE))

    install.packages("BiocManager")

BiocManager::install("scDD")

## limma (http://www.bioconductor.org/packages/release/bioc/html/limma.html)
if (!requireNamespace("BiocManager", quietly = TRUE))

    install.packages("BiocManager")

BiocManager::install("limma")

## Seurat (http://www.bioconductor.org/packages/release/bioc/html/Seurat.html)
if (!requireNamespace("BiocManager", quietly = TRUE))

    install.packages("BiocManager")

BiocManager::install("Seurat")

## zingeR (https://github.com/statOmics/zingeR/)
library(devtools)

install_github("statOmics/zingeR")

## SingleCellExperiment (http://www.bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)
if (!requireNamespace("BiocManager", quietly = TRUE))

    install.packages("BiocManager")

BiocManager::install("SingleCellExperiment")

## scater (http://www.bioconductor.org/packages/release/bioc/html/scater.html)
if (!requireNamespace("BiocManager", quietly = TRUE))

    install.packages("BiocManager")

BiocManager::install("scater")

## aggregation (https://cran.r-project.org/web/packages/aggregation/)
install.packages("aggregation")

 After upper step, you can use the following commands to install scDEA from Github.
 
 # Step 1. Install the devtools package. Invoke R and then type
 
 install.packages("devtools")
 
 # Step 2. Load the devtools package.
 
 library("devtools")
 
 # Step 3. Install the scDEA package from Github.
 
 devtools::install_github("Zhangxf-ccnu/scDEA")
 
 
 
# Step 4.  Useage
 Load the library scDEA in the R console, by running

 library("scDEA")
 
 
 # Step 5. Simply run the scDEA on the Grun datasets.

 data("Grun.counts.matrix")

 data("Grun.group.information")

 Pvals <- scDEA_individual_methods(raw.count = Grun.counts.matrix,

 cell.label = Grun.group.information, DESeq2.parallel = FALSE, scDD = FALSE, monocle.parallel = FALSE)

 combination.Pvals <- lancaster.combination(Pvals, weight = TRUE, trimmed = 0.2)

 adjusted.Pvals <- p.adjust(combination.Pvals, adjusted.method = "bonferroni")
 
Please do not hesitate to contact Dr. Xiao-Fei Zhang at zhangxf@mail.ccnu.edu.cn to seek any clarifications regarding any contents or operation of the archive.
 
 
