#' Gene expression count matrix  of Grun
#'
#' We obtain the preprocessing Grun data from Soneson. Two group cell types, WG and YPS, are
#' selected to perform DE analysis.  We use the function FindVariableFeatures in Seurat R package
#' to select 20000 highly variable genes. In this study, the Grun data contains 20000 genes, 338 cells of WG and
#' 378 cells of YPS.
#'
#' @name Grun.counts.matrix
#'
#' @docType data
#'
#' @usage data(Grun.counts.matrix)
#'
#' @keywords datasets
#' @format a large matrix
#' @examples
#' data(Grun.counts.matrix)
"Grun.counts.matrix"

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

