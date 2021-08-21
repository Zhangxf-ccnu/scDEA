#' Data process
#'
#' This function focus on dealing various single-cell RNA-seq input and unifying output format. \cr
#'
#'
#' @importFrom monocle relative2abs
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#'
#'
#' @param Data single-cell RNA-seq matrix. The format could be raw-counts, FPKM/RPKM, TPM or UMI-counts. The matrix need include gene names and cell names.
#' @param group group information. The cell need be divided into two category.
#' @param norm.form character item. We provide several normalized method for raw-counts data. The method include "TMM","RLE", "CPM". The default is "CPM".
#' @param is.normalized logical. A logical flag to determin whether or not the input dataset normalizes. If TRUE, we will take the Data as normcounts and input for downstream analysis. If not, we provide method for the process.
#'
#' @return
#' \itemize{
#'   \item \strong{sce} :  A \code{\linkS4class{SingleCellExperiment}} item. The object include expression matrix, group information. The expression matrix contains counts and normcounts.
#' }
#' @examples
#' data("Grun.counts.hvg")
#' data("Grun.group.information")
#' sce <- data_process(Data = Grun.counts.hvg, group = Grun.group.information)
#' @details
#' We take \code{\link[monocle]{relative2abs}} transfering relative expression values into absolute transcript counts.
#' However, the process maybe break the original dataset statistical properties. Hence, we advise user don't normalize firstly.
#'
#' @author Huisheng, Li,  <lihs@mails.ccnu.edu.cn>
#' @export
data_process <- function(Data, group, norm.form = "CPM",  is.normalized = FALSE){
  options(warn = -1)

  if(is.normalized){

    normcounts <- Data
    gene_df <- data.frame(Gene = rownames(Data))
    cell_df <- data.frame(cell = colnames(Data))
    # pd <- new("AnnotatedDataFrame", data = cell_df)
    # fd <- new("AnnotatedDataFrame", data = gene_df)
    # transfer <- new("CellDataSet", exprs = as.matrix(Data))
    cds <- monocle::newCellDataSet(cellData = Data, expressionFamily = tobit())
    counts_relative <- monocle::relative2abs(cds)
    counts_relative <- floor(counts_relative)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts_relative, normcounts = normcounts))
    gene_df <- DataFrame(Gene = rownames(sce))
    cell_df <- DataFrame(label = group, cell = colnames(sce))
    rownames(gene_df) <- gene_df$Gene
    rownames(cell_df) <- cell_df$cell
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts_relative, normcounts = normcounts),
                                                      colData = cell_df,
                                                      rowData = gene_df)
    # sce <- scater::calculateQCMetrics(sce)
  } else
  {
    normcounts <- normalized(Data, method = norm.form)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = Data, normcounts = normcounts))
    gene_df <- DataFrame(Gene = rownames(sce))
    cell_df <- DataFrame(label = group, cell = colnames(sce))
    rownames(gene_df) <- gene_df$Gene
    rownames(cell_df) <- cell_df$cell
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = Data, normcounts = normcounts),
                                                      colData = cell_df,
                                                      rowData = gene_df)
    # sce <- scater::calculateQCMetrics(sce)
  }

  return(sce)
}


#' Normalized process
#'
#' The function provide several normalized methods
#' @importFrom edgeR calcNormFactors
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom stats median
#'
#'
#' @param counts.matrix count expression matrix
#' @param method character. "TMM", "RLE", "CPM". The default value is "CPM".
#'
#'
#'
normalized <- function(counts.matrix, method = "CPM"){
  if(method == "TMM"){
    norm_factor <- edgeR::calcNormFactors(counts.matrix, method = method)
    norm.item <- t(t(counts.matrix)/norm_factor)
    return(norm.item)
  }
  if(method == "RLE"){
    geomeans <- exp(rowMeans(log(counts.matrix)))
    SF <-function(cnts){
      stats::median((cnts/geomeans)[is.finite(geomeans) & geomeans >0])
    }
    norm_factor <- apply(counts.matrix, 2, SF)
    norm.item <- t(t(counts.matrix)/norm_factor)
    return(norm.item)
  }
  if(method == "CPM"){
    # sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts.matrix))
    # gene_df <- DataFrame(Gene = rownames(sce))
    # cell_df <- DataFrame(cell = colnames(sce))
    # rownames(gene_df) <- gene_df$Gene
    # rownames(cell_df) <- cell_df$cell
    # sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts.matrix),
    #                                                   colData = cell_df,
    #                                                   rowData = gene_df)
    # norm.item <- scater::calculateCPM(sce)
    norm_factor <- colSums(counts.matrix)
    norm.item <- t(t(counts.matrix)/norm_factor) * 1e6
    return(norm.item)
  }
  # if(method == "TPM"){
  #   sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts.matrix))
  #   gene_df <- DataFrame(Gene = rownames(sce))
  #   cell_df <- DataFrame(cell = colnames(sce))
  #   rownames(gene_df) <- gene_df$Gene
  #   rownames(cell_df) <- cell_df$cell
  #   sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts.matrix),
  #                                                     colData = cell_df,
  #                                                     rowData = gene_df)
  #   norm.item <- scater::calculateTPM(sce, exprs_values = "counts")
  #   return(norm.item)
  # }
}

# perform BPSC
Execute_BPSC <- function(object, coef = 2, BPSC.parallel = TRUE,...){

  object_BPSC <- as.matrix(SummarizedExperiment::assay(object, "normcounts"))

  controllds <- which(object$label == levels(factor(object$label)[1]))
  design <- model.matrix(~ object$label)
  resbp <- BPSC::BPglm(data = object_BPSC, controlIds = controllds,
                       design = design, coef = coef, estIntPar = FALSE, useParallel = BPSC.parallel,...)
  FDR <- p.adjust(resbp$PVAL, method = "BH")
  result_BPSC <- list(gene_names = names(resbp$PVAL), pvalue = resbp$PVAL, FDR = FDR)
  return(result_BPSC)
}


# perform DEsingle
Execute_DEsingle <- function(object, DEsingle.parallel = TRUE){

  object_DEsingle <- as.matrix(SummarizedExperiment::assay(object, "counts"))
  results_DEsingle <- DEsingle::DEsingle(counts = object_DEsingle,
                                         group = factor(object$label),
                                         parallel = DEsingle.parallel)
  result_DEsingle_DE <- list(gene_names = row.names(results_DEsingle), pvalue = results_DEsingle$pvalue,
                             FDR = results_DEsingle$pvalue.adj.FDR)
  return(result_DEsingle_DE)
}


# perform DESeq2
Execute_DESeq2 <- function(object, DESeq2.test = "Wald", beta.Prior = TRUE, DESeq2.parallel = TRUE, DESeq2.fitType = "parametric",...){
  options(warn = -1)
  object_DESeq2 <- as.matrix(SummarizedExperiment::assay(object, "counts"))
  object_DESeq2 <- DESeq2::DESeqDataSetFromMatrix(countData = round(object_DESeq2 + 1),
                                                  colData = data.frame(condition = factor(object$label)),
                                                  design = ~ condition)
  if(DESeq2.test == "LRT"){
    object_DESeq2 <- DESeq2::DESeq(object_DESeq2, test = DESeq2.test, parallel = DESeq2.parallel, betaPrior = FALSE, fitType = DESeq2.fitType, reduced = ~ 1)
  } else{
    object_DESeq2 <- DESeq2::DESeq(object_DESeq2, test = DESeq2.test, parallel = DESeq2.parallel, betaPrior = beta.Prior, fitType = DESeq2.fitType)
  }
  # object_DESeq2 <- DESeq2::DESeq(object_DESeq2, test = DESeq2.test, parallel = DESeq2.parallel, ...)
  res_cpm <- DESeq2::results(object_DESeq2)
  result_DESeq2 <- list(gene_names = rownames(res_cpm),
                        pvalue = res_cpm$pvalue,
                        FDR = res_cpm$padj)
  return(result_DESeq2)
}


#perform MAST
Execute_MAST <- function(object, method_MAST = "bayesglm", MAST.parallel = TRUE, ...){

  options(warn = -1)
  object_MAST <- as.matrix(SummarizedExperiment::assay(object, "normcounts"))

  grp <- object$label
  names(grp) <- colnames(object_MAST)

  sca <- MAST::FromMatrix(exprsArray = log2(object_MAST + 1),
                          cData = data.frame(wellKey = names(grp),
                                             grp = grp))
  zlmdata <- MAST::zlm(~ grp, sca, method = method_MAST, parallel = MAST.parallel, ...)
  mast <- MAST::lrTest(zlmdata, "grp")
  FDR  <- p.adjust(mast[, "hurdle", "Pr(>Chisq)"], method = "BH")

  result_MAST <- list(gene_names = names(mast[, "hurdle", "Pr(>Chisq)"]),
                      pvalue = mast[, "hurdle", "Pr(>Chisq)"],
                      FDR = FDR)
  return(result_MAST)
}


# perform monocle
Execute_monocle <- function(object, monocle.cores = 1,...){
  options(warn = -1)
  objects_monocle <- as.matrix(SummarizedExperiment::assay(object, "counts"))


  counts <-  floor(objects_monocle)
  object_monocle <- monocle::newCellDataSet(as.matrix(counts),
                                            phenoData = new("AnnotatedDataFrame",
                                                            data = data.frame(condition = object$label,
                                                                              row.names = colnames(counts))),
                                            lowerDetectionLimit = 0.5,  expressionFamily = VGAM::negbinomial.size())


  object_monocle <- BiocGenerics::estimateSizeFactors(object_monocle)
  object_monocle <- BiocGenerics::estimateDispersions(object_monocle)

  result_monocle <- monocle::differentialGeneTest(object_monocle, fullModelFormulaStr = "~ condition", cores = monocle.cores,...)


  result_monocle <- list(gene_names = rownames(result_monocle),
                         pvalue = result_monocle$pval,
                         FDR = result_monocle$qval)
  return(result_monocle)
}



# perform scDD
Execute_scDD <- function(object, alpha1 = 0.01, mu0 = 0, s0 = 0.01, a0 = 0.01, b0 = 0.01, scDD.permutation = 0, ...){

  object_scDD <- as.matrix(SummarizedExperiment::assay(object, "normcounts"))

  ### creat scDD input
  condition <- object$label
  names(condition) <- colnames(object)

  SDSumExp <- SingleCellExperiment(assays = list(normcounts = object_scDD),
                                   colData = data.frame(condition))

  prior_param= list(alpha = alpha1, mu0 = mu0, s0 = s0, a0 = a0, b0 = b0)
  scDatExSim <- scDD::scDD(SDSumExp, prior_param=prior_param, condition = "condition", min.nonzero = NULL, min.size = 3, permutations = scDD.permutation, ...)
  res <- scDD::results(scDatExSim)
  result_scDD <- list(gene_names = as.vector(res$gene),
                      pvalue = res$combined.pvalue,
                      FDR = res$combined.pvalue.adj)
  return(result_scDD)
}






# perform Wilcoxon test
Execute_Wilcoxon <- function(object){

  cpm <- as.matrix(SummarizedExperiment::assay(object, "normcounts"))
  groups <- factor(object$label)
  idx <- 1:nrow(cpm)
  names(idx) <- rownames(cpm)
  wilcox_p <- sapply(idx, function(i){
    wilcox.test(cpm[i,]~ groups)$p.value
  })
  FDR <- p.adjust(wilcox_p, method = "BH")
  result_Wilcoxon <- list(gene_names = names(wilcox_p),
                          pvalue = wilcox_p,
                          FDR = FDR)
  return(result_Wilcoxon)
}


# perform T-test
Execute_Ttest <- function(object){

  cpm <- as.matrix(SummarizedExperiment::assay(object, "normcounts"))
  logcpm <- log2(cpm + 1)
  groups <- factor(object$label)
  idx <- seq_len(nrow(logcpm))
  names(idx) <- rownames(logcpm)
  ttest_p <- sapply(idx, function(i){
    t.test(logcpm[i,]~ groups)$p.value
  })
  FDR <- p.adjust(ttest_p, method = "BH")
  result_Ttest <- list(gene_names = names(ttest_p),
                       pvalue = ttest_p,
                       FDR = FDR)
  return(result_Ttest)
}


# perform edgeR
Execute_edgeR <- function(object, Test = "QLFT"){

  object_edgeR <- as.matrix(SummarizedExperiment::assay(object, "counts"))

  dge <- edgeR::DGEList(object_edgeR, group = factor(object$label))
  dge <- edgeR::calcNormFactors(dge)
  groups <- object$label
  design <- stats::model.matrix(~groups)
  dge <- edgeR::estimateDisp(dge, design = design)
  if (Test == "LRT"){
    fit <- edgeR::glmFit(dge, design = design)
    lrt <- edgeR::glmLRT(fit)
  }
  if(Test == "QLFT"){
    fit <- edgeR::glmQLFit(dge, design = design)
    lrt <- edgeR::glmQLFTest(fit)
  }
  tt <- edgeR::topTags(lrt, n = nrow(object_edgeR))
  result_edgeR <- list(gene_names = rownames(tt$table),
                       pvalue = tt$table$PValue,
                       FDR = tt$table$FDR)
  return(result_edgeR)
}





# Execute Limma differential expression analysis
Execute_limma <- function(object, method.fit = "ls", limma.trend =TRUE, limma.robust = TRUE){

  object_Limma <- as.matrix(SummarizedExperiment::assay(object, "counts"))

  dge <- edgeR::DGEList(object_Limma, group = object$label)
  dge <- edgeR::calcNormFactors(dge)
  design <- stats::model.matrix(~object$label)
  y <- methods::new("EList")
  y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  fit <- limma::lmFit(y, design = design, method = method.fit)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  tt <- limma::topTable(fit, n = Inf, adjust.method = "BH")
  result_limma <- list(gene_names = rownames(tt),
                       pvalue = tt$P.Value,
                       FDR = tt$adj.P.Val)
  return(result_limma)
}




# perform Seurat
Execute_Seurat <- function(object,  method.Seurat = "bimod"){

  object_Seurat <- as.matrix(SummarizedExperiment::assay(object, "counts"))

  tmp <- object$label
  names(tmp) <- colnames(object_Seurat)
  meta.data <- data.frame(groups = tmp)

  Seurat.input <- Seurat::CreateSeuratObject(counts = object_Seurat, meta.data = meta.data)
  Seurat.input <- Seurat::NormalizeData(object = Seurat.input)
  res <- Seurat::FindMarkers(object = Seurat.input, ident.1 = levels(as.factor(tmp))[1],
                             ident.2 = levels(as.factor(tmp))[2], group.by = 'groups',
                             logfc.threshold = -Inf, test.use = method.Seurat,
                             only.pos = FALSE, verbose = FALSE)
  results_Seurat <- list(gene_names = row.names(res),
                         pvalue = res$p_val,
                         FDR = res$p_val_adj)

  return(results_Seurat)
}

# perform zingeR.edgeR
Execute_zingeR.edgeR <- function(object, maxit.EM = 200){

  object_zingerR.edgeR <- as.matrix(SummarizedExperiment::assay(object, "counts"))

  dge <- edgeR::DGEList(round(object_zingerR.edgeR), group = factor(object$label))
  dge <- calcNormFactors(dge)
  groups <- object$label
  design <- stats::model.matrix(~groups)
  weights <- zeroWeightsLS(counts=dge$counts, design=design, maxit = maxit.EM, normalization="TMM")
  dge$weights <- weights

  dge = edgeR::estimateDisp(dge, design)
  fit = edgeR::glmFit(dge, design)
  lrt = zingeR::glmWeightedF(fit, coef=2, independentFiltering = TRUE)

  result_zingerR.edgeR <- list(gene_names = rownames(lrt$table),
                               pvalue = lrt$table$PValue,
                               FDR = lrt$table$padjFilter)
}


#### measure the jaccard index of differential expression analysis
jaccard_index <- function(results1, results2, cutoff = c(200, 400, 600, 1000)){
  if(length(results1) != length(results2))
    stop("The results of differential expression analysis need same length!")
  jaccard.results <- matrix(0, nrow = length(cutoff), ncol = length(results1))
  rownames(jaccard.results) <- paste0("cutoff", cutoff)
  colnames(jaccard.results) <- names(results1)
  for (i in 1:length(cutoff)) {
    for (j in 1:length(results1)) {
      index_smmoothing1 <- sort(results1[[j]],  decreasing = FALSE)
      gene_smoothing1 <- names(index_smmoothing1)[1:cutoff[i]]

      index_smmoothing2 <- sort(results2[[j]],  decreasing = FALSE)
      gene_smoothing2 <- names(index_smmoothing2)[1:cutoff[i]]


      # jaccard.results[i,j] <- length(intersect(gene_smoothing1, gene_smoothing2)) / length(union(gene_smoothing1, gene_smoothing2))
      jaccard.results[i,j] <- sum(gene_smoothing1 %in% gene_smoothing2) / length(union(gene_smoothing1, gene_smoothing2))
      # jaccard.results[i,j] <- jaccard::jaccard(gene_smoothing1, gene_smoothing2)
    }

  }
  return(jaccard.results)
}



#' Run individual DE analysis methods to perform DE analysis on scRNA-seq datasets.
#'
#' This function is implemented to perform individual DE analysis methods.  The current implementation of
#' scDEA integrates tweleve state-of-the-art methods:  Beta-poisson mixture model (BPSC),
#' DEsingle, DESeq2, edgeR, Model-based analysis of single-cell transcriptomics (MAST), monocle, scDD,
#' T-test, Wilcoxon rank sum test (Wilcoxon test), limma, Seurat and zingeR.edgeR.
#' This function depends on the follwing R package: BPSC, DEsingle, DESeq2, edgeR, MAST, monocle, scDD,
#' limma, Seurat, zingeR, SingleCellExperiment, dplyr.
#' These packages will be automatically installed along
#' with scDEA. \cr
#'
#' @importFrom monocle relative2abs
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom stats median p.adjust model.matrix wilcox.test
#' @importFrom BPSC BPglm
#' @importFrom DEsingle DEsingle DEtype
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom MAST FromMatrix zlm lrTest
#' @importFrom BiocGenerics estimateSizeFactors estimateDispersions
#' @importFrom methods new
#' @importFrom VGAM negbinomial.size tobit
#' @importFrom scDD preprocess scDD
#' @importFrom S4Vectors DataFrame
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmFit glmLRT glmQLFit glmQLFTest topTags dimnames.TopTags
#' @importFrom limma lmFit eBayes topTable
#' @importFrom Seurat CreateSeuratObject FindMarkers
#' @importFrom zingeR glmWeightedF zeroWeightsLS
#' @importFrom stats cor t.test
#' @importFrom SummarizedExperiment assay
#'
#' @param raw.count single-cell RNA-seq matrix. The format could be raw read count or normalized matrix. The rows correspond to genes and the columns.
#' @param cell.label  cell labels information. The cells need be divided into two categories.
#' @param is.normalized a boolean variable that defines whether the input raw.count has been normalized? Default is FALSE.
#' @param verbose a boolean variable that defines whether to save the DE analysis results and name "Results_DE_individual.RData" in the current working directory.
#' @param BPSC a boolean variable that defines whether to perform DE analysis using the BPSC method. Default is TRUE.
#' @param DEsingle a boolean variable that defines whether to perform DE analysis using the DEsingle method. Default is TRUE.
#' @param DESeq2 a boolean variable that defines whether to perform DE analysis using the DESeq2 method. Default is TRUE.
#' @param edgeR a boolean variable that defines whether to perform DE analysis using the edgeR method. Default is TRUE.
#' @param MAST a boolean variable that defines whether to perform DE analysis using the MAST method. Default is TRUE.
#' @param monocle a boolean variable that defines whether to perform DE analysis using the MONOCLE method. Default is TRUE.
#' @param scDD a boolean variable that defines whether to perform DE analysis using the scDD method. Default is TRUE.
#' @param Ttest a boolean variable that defines whether to perform DE analysis using the T-test method. Default is TRUE.
#' @param Wilcoxon a boolean variable that defines whether to perform DE analysis using the Wilcoxon method. Default is TRUE.
#' @param limma a boolean variable that defines whether to perform DE analysis using the limma method. Default is TRUE.
#' @param Seurat a boolean variable that defines whether to perform DE analysis using the Seurat method. Default is TRUE.
#' @param zingeR.edgeR a boolean variable that defines whether to perform DE analysis using the zingeR.edgeR method. Default is TRUE.
#' @param BPSC.coef an integer to point out the column index corresponding to the coefficient for the generalized linear mode (GLM) testing in BPSC. Default value is 2.
#' @param BPSC.normalize a string variable specifying the type of size factor estimation in BPSC method. Possible values: "TMM", "RLE", "CPM". Default is "CPM".
#' @param BPSC.parallel a boolean variable that defines whether to execute parallel computing for BPSC method. Default is TRUE.
#' @param DEsingle.parallel a boolean variable that defines whether to execute parallel computing for DEsingle method. Default is TRUE.
#' @param DEsingle.normalize a string variable specifying the type of size factor estimation in DEsingle method. Possible values: "TMM", "RLE", "CPM". Default is "CPM".
#' @param DESeq2.test  a string variable specifying the type of test the difference in deviance between a full and reduced model formula in DESeq2 method. Possible values: "Wald" or "LRT". The values represent Wald tests or likelihood ratio test. Default is "Wald".
#' @param DESeq2.parallel a boolean variable that defines whether to execute parallel computing for DESeq2 method. Default is TRUE. The parallel computing may fail on Windows system. Default is TRUE.
#' @param DESeq2.beta.prior a boolean variable that defines whether or not to put a zero-mean normal prior on the non-intercept coefficient in DESeq2 method. Default is TRUE.
#' @param DESeq2.fitType a string variable specifying the type of fitting of dispersions to the mean intensity in DESeq2 method. Possible values: "parametric", "local", "mean". Default is "parametric".
#' @param DESeq2.normalize a string variable specifying the type of size factor estimation in DESeq2 method. Possible values: "TMM", "RLE", "CPM". Default is "CPM".
#' @param edgeR.Test a string variable specifying the type of fitting distribution to count data for each gene. Possible values: "LRT", "QLFT". The values represent negative binomial generalized log-linear model and quasi-likelihood negative binomial generalized log-linear model. Default is "QLFT".
#' @param edgeR.normalize a string variable specifying the type of size factor estimation in edgeR method. Possible values: "TMM", "RLE", "CPM". Default is "CPM".
#' @param limma.method.fit a string variable specifying the type of fitting method in limma method. Possible values: "ls", "robust". The values represent least squares and robust regression. Default is "ls".
#' @param limma.trend a boolean variable that defines whether or not to allow an intensity-trend for the prior variance in limma method. Default is TRUE.
#' @param limma.robust a boolean variable that defines whether or not to estimate defined prior information and variance prior against outlier sample variances in limma method. Default is TRUE.
#' @param limma.normalize  a string variable specifying the type of size factor estimation in limma method. Possible values: "TMM", "RLE", "CPM". Default is "CPM".
#' @param Seurat.normalize  a string variable specifying the type of size factor estimation in Seurat method. Possible values: "TMM", "RLE", "CPM". Default is "CPM".
#' @param Seurat.method a string variable specifying the type of test method in Seurat method. Possible values: "LR", "bimod", "roc". The values represent likelihood-ratio test, negative binomial generalized linear model, ROC analysis. Default is "bimod".
#' @param MAST.method a string variable specifying the type of test method in MAST method. Possible values: "glm", "glmer", "bayesglm". Default is "bayesglm".
#' @param MAST.normalize a string variable specifying the type of size factor estimation in Seurat method. Possible values: "TMM", "RLE", "CPM". Default is "CPM".
#' @param MAST.parallel a boolean variable that defines whether to execute parallel computing for MAST method. Default is TRUE.
#' @param monocle.cores the number of cores to be used while testing each gene for differential expression.. Default is 1.
#' @param monocle.normalize a string variable specifying the type of size factor estimation in monocle method. Possible values: "TMM", "RLE", "CPM". Default is "CPM".
#' @param scDD.alpha1 prior parameter value to be used to model each gene as a mixture of DP normals in scDD method. Default is 0.01.
#' @param scDD.mu0 prior parameter values to be used to model each gene as a mixture of DP normals in scDD method. Default is 0.
#' @param scDD.s0 prior parameter values to be used to model each gene as a mixture of DP normals in scDD method. Default is 0.01.
#' @param scDD.a0 prior parameter values to be used to model each gene as a mixture of DP normals in scDD method. Default is 0.01.
#' @param scDD.b0 prior parameter values to be used to model each gene as a mixture of DP normals in scDD method. Default is 0.01.
#' @param scDD.normalize  a string variable specifying the type of size factor estimation in scDD method. Possible values:  "TMM", "RLE", "CPM". Default is "CPM".
#' @param scDD.permutation the number of permutations to be used in calculating empirical p-values in scDD method. If the parameter value is set to 0, the full Bayes Factor will not be performed. Else, scDD method takes the nonparametric Kolmogorove-Smirnov test to identify DGEs. Default is 0.
#' @param Ttest.normalize a string variable specifying the type of size factor estimation in t-test method. Possible values: "TMM", "RLE", "CPM". Default is "CPM".
#' @param Wilcoxon.normalize a string variable specifying the type of size factor estimation in Wilcoxon method. Possible values: "TMM", "RLE", "CPM". Default is "CPM".
#' @param zingeR.edgeR.normalize a string variable specifying the type of size factor estimation in zingeR.edgeR method. Possible values: "TMM", "RLE", "CPM". Default is "CPM".
#' @param zingeR.edgeR.maxit.EM The number of iterations for EM-algorithm in zingeR.edgeR method. If the EM-algorithm does not stop automatically, then, the algorithm may not be convergence. The user need set a larger value. Default is 100.
#'
#' @return a p-values matrix contains the p-values of each differential expression anlysis methods.
#'
#' @examples
#' data("Grun.counts.hvg")
#' data("Grun.group.information")
#' # scDD is very slow
#' Pvals <- scDEA_individual_methods(raw.count = Grun.counts.hvg,
#' cell.label = Grun.group.information,  verbose = FALSE)
#' combination.Pvals <- lancaster.combination(Pvals, weight = TRUE, trimmed = 0.2)
#' adjusted.Pvals <- scDEA.p.adjust(combination.Pvals, adjusted.method = "bonferroni")
#' @export
#'
#' @author Huisheng, Li,  <lihs@mails.ccnu.edu.cn>
#'


scDEA_individual_methods <- function(raw.count, cell.label,
                                     is.normalized = FALSE, verbose = TRUE,
                                     BPSC = TRUE, DEsingle = TRUE, DESeq2 = TRUE, edgeR = TRUE,
                                     MAST = TRUE, monocle = TRUE, scDD = TRUE, Ttest = TRUE, Wilcoxon = TRUE,
                                     limma = TRUE, Seurat = TRUE, zingeR.edgeR = TRUE,
                                     BPSC.coef = 2, BPSC.normalize = "CPM", BPSC.parallel = TRUE,
                                     DEsingle.parallel = TRUE, DEsingle.normalize = "CPM",
                                     DESeq2.test = "LRT", DESeq2.parallel = TRUE, DESeq2.beta.prior = TRUE, DESeq2.fitType = "parametric", DESeq2.normalize = "CPM",
                                     edgeR.Test = "QLFT", edgeR.normalize = "TMM",
                                     limma.method.fit = "ls", limma.trend = TRUE, limma.robust = TRUE, limma.normalize = "CPM",
                                     Seurat.normalize = "CPM", Seurat.method = "bimod",
                                     MAST.method = "bayesglm", MAST.normalize = "CPM", MAST.parallel = TRUE,
                                     monocle.cores = 1, monocle.normalize = "CPM",
                                     scDD.alpha1 = 0.01, scDD.mu0 = 0, scDD.s0 = 0.01, scDD.a0 = 0.01, scDD.b0 = 0.01, scDD.normalize = "CPM", scDD.permutation = 0,
                                     Ttest.normalize = "CPM",
                                     Wilcoxon.normalize = "CPM",
                                     zingeR.edgeR.normalize = "CPM", zingeR.edgeR.maxit.EM = 100
){
  if(is.normalized)
    cat("Some differential expression analysis methods are based on counts.
          Hence, we don't advise user normalizes firstly. However, if you only have the dataset,
          in order to keep the unity of input data form, we takes relative2abs() function transfer to
          mRNA counts. Cautious, the results may be not precise. \n")

  Methods <- c("BPSC", "DEsingle", "DESeq2", "edgeR", "MAST", "monocle", "scDD","Ttest", "Wilcoxon", "limma", "Seurat", "zingeR.edgeR")
  Methods.idx <- c(BPSC,  DEsingle, DESeq2, edgeR, MAST, monocle, scDD, Ttest, Wilcoxon, limma, Seurat, zingeR.edgeR)


  Methods.used = Methods[Methods.idx]

  K <- length(Methods.used)
  p <- dim(raw.count)[1]

  Results.DE.individual <- list()

  k <- 1
  #BPSC
  if(BPSC == TRUE){
    cat("Execute BPSC analysis....\n")
    sce.BPSC <- data_process(Data = raw.count, group = cell.label, norm.form = BPSC.normalize,  is.normalized = is.normalized)

    pred.BPSC <- Sys.time()
    results.BPSC <- Execute_BPSC(object = sce.BPSC, coef = BPSC.coef, BPSC.parallel = BPSC.parallel)
    end.BPSC <- Sys.time()
    time.BPSC <- difftime(end.BPSC, pred.BPSC, units = "mins")
    cat("Run time for BPSC: ", time.BPSC, "min","\n")

    temp <- results.BPSC$pvalue
    names(temp) <- results.BPSC$gene_names
    Results.DE.individual$BPSC <- temp
    k <- k + 1
  }
  #DEsingle
  if(DEsingle == TRUE){
    cat("Execute DEsingle analysis....\n")
    sce.DEsingle <- data_process(Data = raw.count, group = cell.label, norm.form = DEsingle.normalize,  is.normalized = is.normalized)

    pred.DEsingle <- Sys.time()
    results.DEsingle <- Execute_DEsingle(object = sce.DEsingle, DEsingle.parallel = DEsingle.parallel)
    end.DEsingle <- Sys.time()
    time.DEsingle <- difftime(end.DEsingle, pred.DEsingle, units = "mins")
    cat("Run time for DEsingle: ", time.DEsingle, "min","\n")

    temp <- results.DEsingle$pvalue
    names(temp) <- results.DEsingle$gene_names
    Results.DE.individual$DEsingle <- temp
    k <- k + 1
  }
  #DESeq2
  if(DESeq2 == TRUE){
    cat("Execute DESeq2 analysis....\n")
    sce.DESeq2 <- data_process(Data = raw.count, group = cell.label, norm.form = DESeq2.normalize,  is.normalized = is.normalized)

    pred.DESeq2 <- Sys.time()
    results.DESeq2 <- Execute_DESeq2(sce.DESeq2, DESeq2.test = DESeq2.test, DESeq2.parallel = DESeq2.parallel, beta.Prior = DESeq2.beta.prior, DESeq2.fitType = DESeq2.fitType)
    end.DESeq2 <- Sys.time()
    time.DESeq2 <- difftime(end.DESeq2, pred.DESeq2, units = "mins")
    cat("Run time for DESeq2: ", time.DESeq2, "min","\n")

    temp <- results.DESeq2$pvalue
    names(temp) <- results.DESeq2$gene_names
    Results.DE.individual$DESeq2 <- temp
    k <- k + 1
  }
  #edgeR
  if(edgeR == TRUE){
    cat("Execute edgeR analysis....\n")
    sce.edgeR <- data_process(Data = raw.count, group = cell.label, norm.form = edgeR.normalize,  is.normalized = is.normalized)

    pred.edgeR <- Sys.time()
    results.edgeR <- Execute_edgeR(sce.edgeR, Test = edgeR.Test)
    end.edgeR <- Sys.time()
    time.edgeR <- difftime(end.edgeR, pred.edgeR, units = "mins")
    cat("Run time for edgeR: ", time.edgeR, "min","\n")

    temp <- results.edgeR$pvalue
    names(temp) <- results.edgeR$gene_names
    Results.DE.individual$edgeR <- temp
    k <- k + 1
  }
  #MAST
  if(MAST == TRUE){
    cat("Execute MAST analysis....\n")
    sce.MAST <- data_process(Data = raw.count, group = cell.label, norm.form = MAST.normalize,  is.normalized = is.normalized)

    pred.MAST <- Sys.time()
    results.MAST <- Execute_MAST(sce.MAST,  method_MAST = MAST.method, MAST.parallel = MAST.parallel)
    end.MAST <- Sys.time()
    time.MAST <- difftime(end.MAST, pred.MAST, units = "mins")
    cat("Run time for MAST: ", time.MAST, "min","\n")


    temp <- results.MAST$pvalue
    names(temp) <- results.MAST$gene_names
    Results.DE.individual$MAST <- temp
    k <- k + 1
  }
  #MONOCLE

  if(monocle == TRUE){
    cat("Execute monocle analysis....\n")
    sce.monocle <- data_process(Data = raw.count, group = cell.label, norm.form = monocle.normalize,  is.normalized = is.normalized)

    pred.monocle <- Sys.time()
    results.monocle <- Execute_monocle(sce.monocle, monocle.cores = monocle.cores)
    end.monocle <- Sys.time()
    time.monocle <- difftime(end.monocle, pred.monocle, units = "mins")
    cat("Run time for monocle: ", time.monocle, "min","\n")

    temp <- results.monocle$pvalue
    names(temp) <- results.monocle$gene_names
    Results.DE.individual$monocle <- temp
    k <- k + 1
  }
  #scDD
  if(scDD == TRUE){
    cat("Execute scDD analysis....\n")
    sce.scDD <- data_process(Data = raw.count, group = cell.label, norm.form = scDD.normalize,  is.normalized = is.normalized)

    pred.scDD <- Sys.time()
    results.scDD <- Execute_scDD(sce.scDD, alpha1 = scDD.alpha1, mu0 = scDD.mu0, s0 = scDD.s0, a0 = scDD.a0, b0 = scDD.b0, scDD.permutation = scDD.permutation)
    end.scDD <- Sys.time()
    time.scDD <- difftime(end.scDD, pred.scDD, units = "mins")
    cat("Run time for scDD: ", time.scDD, "min","\n")

    temp <- results.scDD$pvalue
    names(temp) <- results.scDD$gene_names
    Results.DE.individual$scDD <- temp
    k <- k + 1
  }
  #Ttest
  if(Ttest == TRUE){
    cat("Execute Ttest analysis....\n")
    sce.Ttest <- data_process(Data = raw.count, group = cell.label, norm.form = Ttest.normalize,  is.normalized = is.normalized)

    pred.Ttest <- Sys.time()
    results.Ttest <- Execute_Ttest(sce.Ttest)
    end.Ttest <- Sys.time()
    time.Ttest <- difftime(end.Ttest, pred.Ttest, units = "mins")
    cat("Run time for Ttest: ", time.Ttest, "min","\n")


    temp <- results.Ttest$pvalue
    names(temp) <- results.Ttest$gene_names
    Results.DE.individual$Ttest <- temp
    k <- k + 1
  }
  #Wilcoxon
  if(Wilcoxon == TRUE){
    cat("Execute Wilcoxon analysis....\n")
    sce.Wilcoxon <- data_process(Data = raw.count, group = cell.label, norm.form = Wilcoxon.normalize,  is.normalized = is.normalized)

    pred.Wilcoxon <- Sys.time()
    results.Wilcoxon <- Execute_Wilcoxon(sce.Wilcoxon)
    end.Wilcoxon <- Sys.time()
    time.Wilcoxon <- difftime(end.Wilcoxon, pred.Wilcoxon, units = "mins")
    cat("Run time for Wilcoxon: ", time.Wilcoxon, "min","\n")

    temp <- results.Wilcoxon$pvalue
    names(temp) <- results.Wilcoxon$gene_names
    Results.DE.individual$Wilcoxon <- temp
    k <- k + 1
  }
  #limma
  if(limma == TRUE){
    cat("Execute limma analysis....\n")
    sce.limma <- data_process(Data = raw.count, group = cell.label, norm.form = limma.normalize,  is.normalized = is.normalized)

    pred.limma <- Sys.time()
    results.limma <- Execute_limma(sce.limma, method.fit = limma.method.fit, limma.trend = limma.trend, limma.robust = limma.robust)
    end.limma <- Sys.time()
    time.limma <- difftime(end.limma, pred.limma, units = "mins")
    cat("Run time for limma: ", time.limma, "min","\n")

    temp <- results.limma$pvalue
    names(temp) <- results.limma$gene_names
    Results.DE.individual$limma <- temp
    k <- k + 1
  }
  #Seurat
  if(Seurat == TRUE){
    cat("Execute Seurat analysis....\n")
    sce.Seurat <- data_process(Data = raw.count, group = cell.label, norm.form = Seurat.normalize,  is.normalized = is.normalized)

    pred.Seurat <- Sys.time()
    results.Seurat <- Execute_Seurat(sce.Seurat, method.Seurat = Seurat.method)
    end.Seurat <- Sys.time()
    time.Seurat <- difftime(end.Seurat, pred.Seurat, units = "mins")
    cat("Run time for Seurat: ", time.Seurat, "min","\n")

    temp <- results.Seurat$pvalue
    names(temp) <- results.Seurat$gene_names
    Results.DE.individual$Seurat <- temp
    k <- k + 1
  }
  #zingeR.edgeR
  if(zingeR.edgeR == TRUE){
    cat("Execute zingeR.edgeR analysis....\n")
    sce.zingeR.edgeR <- data_process(Data = raw.count, group = cell.label, norm.form = zingeR.edgeR.normalize,  is.normalized = is.normalized)

    pred.zingeR.edgeR <- Sys.time()
    results.zingeR.edgeR <- Execute_zingeR.edgeR(sce.zingeR.edgeR,  maxit.EM = zingeR.edgeR.maxit.EM)
    end.zingeR.edgeR <- Sys.time()
    time.zingeR.edgeR <- difftime(end.zingeR.edgeR, pred.zingeR.edgeR, units = "mins")
    cat("Run time for zingeR.edgeR: ", time.zingeR.edgeR, "min","\n")

    temp <- results.zingeR.edgeR$pvalue
    names(temp) <- results.zingeR.edgeR$gene_names
    Results.DE.individual$zingeR.edgeR <- temp
    k <- k + 1
  }
  gene.names <- rownames(raw.count)
  for (l in 1:length(Results.DE.individual)) {
    Results.DE.individual[[l]][is.na(Results.DE.individual[[l]])] <- 1

    index <- gene.names %in% names(Results.DE.individual[[l]])
    if(sum(index) == length(gene.names)){
      Results.DE.individual[[l]] <- Results.DE.individual[[l]][gene.names]
    } else{
      tmp <- matrix(1, nrow = sum(!index))
      lossing.item <- as.vector(tmp)
      names(lossing.item) <- gene.names[!index]

      integrated <- c(Results.DE.individual[[l]], lossing.item)

      Results.DE.individual[[l]] <- integrated[gene.names]
    }
    Results.DE.individual[[l]][is.na(Results.DE.individual[[l]])] <- 1
  }


  p.values <- matrix(0, nrow = length(Results.DE.individual[[1]]), ncol = length(Results.DE.individual))
  rownames(p.values) <- names(Results.DE.individual[[1]])
  colnames(p.values) <- names(Results.DE.individual)
  for (i in 1:ncol(p.values)) {
    temp <-  Results.DE.individual[[i]][rownames(p.values)]
    p.values[,i] <- temp
  }
  Pvals <- p.values


  if(verbose)
  {
    cat("Saving individual DE analysis results...\n")

    save(Results.DE.individual, file = "Results_DE_individual.RData")
  }
  return(Pvals)
}


#' Combining the p-values identified by DE ananlysis methods
#'
#' This function is used to combine the results of differential expression anlaysis and obtain an uniform result.
#'
#' @importFrom Matrix rowSums
#' @importFrom aggregation lancaster
#'
#' @param Pvals a p-value matrix. The rows represent genes and the columns correspond to the individual DE analysis methods.
#' @param weight a boolean variable that defines whether to use spearman correlation measure the similarity between different DE analysis methods. Default is "TRUE".
#' @param trimmed  a real number between 0 and 0.5 to trim p-values. Default value is 0.2.
#'
#' @return a vector represents the combined p-value for each gene.
#'
#' @export
#'
#' @author Huisheng, Li,  <lihs@mails.ccnu.edu.cn>

lancaster.combination <- function(Pvals, weight = TRUE, trimmed = 0.2){

  #### check if there is NA
  if (sum(is.na(Pvals)) > 0){
    stop("Cannot have NAs in the p-values!")
  }
  #### check if Pvals are between 0 and 1
  if ((sum(Pvals < 0)+sum(Pvals > 1)) > 0){
    stop("P-values must be between 0 and 1!")
  }


  ####### check the validity of the user supplied weights and standadize them.
  if(weight){
    weight.matrix <- spearman.matrix(p.values = Pvals)
    temp <- Matrix::rowSums(weight.matrix)/ncol(weight.matrix)
    weight.pval <- temp
  }else{
    num_pval <- ncol(Pvals)
    weight.pval <- rep(2, num_pval)
    names(weight.pval) <- colnames(Pvals)
  }
  ###### calculate combanation p.values for each gene
  combination.pvalues.lancaster <- function(x, weight){
    location <- x == 1
    temp <- sort(x[! location], decreasing = FALSE)
    num <- floor(length(temp) * trimmed)
    Pvals.trimmed <- temp[(num + 1) : (length(temp) - num)]
    weight.use <- weight[names(Pvals.trimmed)]
    weight.use.norm <- weight.use/sum(weight.use)
    ####  Replace extreme p-values with 10e-320 to obtain an upper bound for the aggregated p-value. ####
    if (any(Pvals.trimmed < 9.99988867182683e-320)) {
      Pvals.trimmed[Pvals.trimmed < 9.99988867182683e-320] <- 9.99988867182683e-320
    }
    pval <- aggregation::lancaster(pvalues = Pvals.trimmed, weights = weight.use.norm)
  }
  pvals <- apply(Pvals, 1, combination.pvalues.lancaster, weight = weight.pval)
  combination.Pvals <- pvals
  return(combination.Pvals)
}

#' adjusting the p-values identified by individual DE ananlysis methods and scDEA
#'
#' This function is used to adjust the p-values generated by individual DE analysis methods and scDEA.
#'
#' @param combination.Pvals a vector of combined p-value for each gene.
#' @param adjusted.method a string variable specifying the type of p.adjust method used in t-test method.
#' Possible values: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". Default is "bonferroni".
#' @return adjusted p-value ,
#'
#'
#' @export

scDEA.p.adjust <- function(combination.Pvals, adjusted.method = "bonferroni"){
  # adjusted.Pvals <- matrix(0, nrow = nrow(Pvals), ncol = (ncol(Pvals) + 1))
  # rownames(adjusted.Pvals) <- rownames(Pvals)
  # colnames(adjusted.Pvals) <- c(colnames(Pvals), "scDEA")
  # for (i in 1:ncol(Pvals)) {
  #   adjusted.Pvals[,i] <- p.adjust(p = Pvals[,i], method = adjusted.method)
  # }
  adjusted.Pvals <- stats::p.adjust(p = combination.Pvals, method = adjusted.method)
  # adjusted.Pvals[,ncol(adjusted.Pvals)] <- temp[rownames(adjusted.Pvals)]
  return(adjusted.Pvals)
}


spearman.matrix <- function(p.values){
  methods.names <- colnames(p.values)
  spearman.correlation <- matrix(0, nrow = ncol(p.values), ncol = ncol(p.values))
  rownames(spearman.correlation) <- methods.names
  colnames(spearman.correlation) <- methods.names
  for (i in 1:ncol(spearman.correlation)) {
    temp.i <- p.values[,i]
    for (j in 1:ncol(spearman.correlation)) {
      temp.j <- p.values[,j]
      spearman.correlation[i,j] <- cor(x = temp.i, y = temp.j, method = c("spearman"))
    }
  }
  diag(spearman.correlation) <- 0
  spearman.correlation[spearman.correlation < 0] <- 0
  return(spearman.correlation)
}
