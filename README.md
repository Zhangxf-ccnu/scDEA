# scDEA
R package applies ensemble learning for single-cell differential expression analysis on single-cell RNA-seq dataset. 
Besides, we also develop a Shiny Application for users. The source code of Shiny Application are published on the (https://github.com/Zhangxf-ccnu/scDEA-shiny) .
In the further, we will deploy a web server for users who do not have a strong programming background. 


The scDEA package has the following R-package dependencies: BPSC, DEsingle, DESeq2, edgeR, MAST, monocle, scDD, limma, Seurat, zingeR, SingleCellExperiment, scater, aggregation. 
For the R-package dependencies, you can load on most of dependencies packages on your R when install the scDEA R package. However,
If the dependencies are not installed correctly, please install them by yourself. In my experience, the scDEA can't install BPSC, DEsingle and zingeR automatically on Windows system but can load all on Ubuntu 16.04.
 Hence, we advise, users should install BPSC, DEsingle and zingeR firstly by the following instruction when you use scDEA on windows system.

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
 
 devtools::install_github("Zhangxf-ccnu/scDEA", subdir = "pkg")
 
 
# Step 4.  Useage
 Load the library scDEA in the R console, by running

 library("scDEA")
 
 
 # Step 5. Simply run the scDEA on the Grun datasets.

 data("Grun.counts.hvg")
 
 data("Grun.group.information")
 
 Pvals <- scDEA_individual_methods(raw.count = Grun.counts.hvg, cell.label = Grun.group.information)
 
 combination.Pvals <- lancaster.combination(Pvals, weight = TRUE, trimmed = 0.2)
 
 adjusted.Pvals <- p.adjust(combination.Pvals, adjusted.method = "bonferroni")
 
 # Step 6. Run Shiny application locally.  
 
 Load on https://github.com/Zhangxf-ccnu/scDEA-shiny ,download the scDEA-shiny-main.zip, run the shiny application  
 
 install.packages("shiny")  
 
 library(shiny)  
 
 shiny::runApp('scDEA-shiny-main')  
 
 
Please do not hesitate to contact Dr. Xiao-Fei Zhang at zhangxf@mail.ccnu.edu.cn to seek any clarifications regarding any contents or operation of the archive.
 
 
