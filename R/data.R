#' Gene expression count matrix  of Grun
#'
#' We obtain the preprocessing Grun data from Soneson. Two group cell types, WG and YPS, are
#' selected to perform DE analysis. To demonstrate the effectiveness of the scDEA, We use the function FindVariableFeatures in Seurat R package
#' to select 2000 highly variable genes. In this study, the Grun data contains 2000 genes, 338 cells of WG and
#' 378 cells of YPS.
#'
#' @name Grun.counts.hvg
#'
#' @docType data
#'
#' @usage data(Grun.counts.hvg)
#'
#' @keywords datasets
#' @format a large matrix
#' @examples
#' data(Grun.counts.hvg)
"Grun.counts.hvg"

#' Cell type lebals of Grun data
#'
#' @name Grun.group.information
#' @docType data
#' @usage data(Grun.group.information)
#' @keywords datasets
#' @format a vector
#' @examples
#' data(Grun.group.information)
"Grun.group.information"

