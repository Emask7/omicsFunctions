#' @title DESeq2 wrapper
#'
#' @description DESeq2 wrapper
#' @param geneCounts data.frame containing 6 columns of raw gene or transcript counts.
#' @param sample_data data.frame containing experimental factrors, with the number of rows equal to the number of samples
#' @param expt_design a formula indicating the factors of the experimental design
#' @param coefficient a number indicating the coefficient to be tested
#' @param count_type string. Either "gene" or "transcript"
#' @export
#' @examples
#' run_DESeq2(geneCounts, sample_data, expt_design, coefficient, count_type)

run_DESeq2 <- function(geneCounts, sample_data, expt_design, coefficient, count_type) {
  # Make a DESeqDataSet object ------------------------------------------------
    geneCounts <- as.matrix(geneCounts)
    geneCounts[, 1:6] <- as.integer(geneCounts)

    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = geneCounts, colData = sample_data, design = expt_design
    )

  # Wald test with LFC shrinkage ----------------------------------------------
    dds <- DESeq2::DESeq(dds, quiet = TRUE)

    shrink <- DESeq2::lfcShrink(
      dds, type = "apeglm", lfcThreshold = 1, coef = coefficient, quiet = TRUE
    )

    print("DESeq2 Wald test with LFC shrinkage DEG summary:", quote = FALSE)
    DESeq2::summary(shrink, 0.05)

    res <- data.frame(rownames(shrink), shrink[, c(2, 4)])
    if (count_type == "gene") colnames(res) <- c("Gene_ID", "LFC", "padj")
    else if (count_type == "transcript") {
      colnames(res) <- c("Transcript_ID", "LFC", "padj")
    } else print("Error: count_type must be either 'gene' or 'transcript'")

    res_filt <- subset(res, res$padj <= 0.05 & abs(res$LFC) >= 1)

  # Return a list object with all of the data ---------------------------------
    list(DESeq_object = dds, results = res, filtered_results = res_filt)
}

#' @title edgeR wrapper
#'
#' @description edgeR wrapper
#' @param gene_counts A data frame containing 6 columns of raw gene counts.
#' @param ortho_list A data frame containing two columns (Gene_ID and Human_Gene_ID)
#' @param LFC_filter logical. If TRUE, DEG lists are filtered by LFC values in addition to FDR values
#' @import edgeR
#' @export
#' @examples
#' run_edgeR(gene_counts, ortho_list)

run_edgeR <- function(gene_counts) {

  # Specify experimental design factors ---------------------------------------
  animal <- factor(c(rep(c("15979", "30760", "31151"), 2)))
  timepoint <- factor(c(rep("T-7", 3), rep("T15", 3)))
  design <- stats::model.matrix(~ animal + timepoint)

  # Make a DGEList object from the raw count data and filter for low reads ----
  dgeList <- DGEList(counts = gene_counts, group = timepoint)
  keep <- filterByExpr(dgeList, min.count = 5)
  dgeList <- dgeList[keep, , keep.lib.sizes = FALSE]

  # Calculate normalization factors -------------------------------------------
  dgeList <- calcNormFactors(dgeList, method = "TMM")

  # # Data exploration: multidimensional scaling plot ---------------------------
  #   if (show_plots) limma::plotMDS(dgeList, main = "MDS Plot")

  # Estimate dispersion -------------------------------------------------------
  dgeList <- estimateDisp(dgeList, design, robust = TRUE)

  # # Data exploration: biological coefficient of variation plot ----------------
  #   if (show_plots) plotBCV(dgeList, main = "BCV Plot")

  # DE testing ----------------------------------------------------------------
  fit <- glmFit(dgeList, design)
  lrtTreat <- glmTreat(fit, lfc = 1, null = "interval")

  # if (show_plots) limma::plotMD(lrtTreat, main = "Mean-Difference Plot")

  res <- data.frame(
    rownames(lrtTreat$table), lrtTreat$table$logFC,
    stats::p.adjust(lrtTreat$table$PValue, method = "fdr")
  )
  colnames(res) <- c("Gene_ID", "LFC", "padj")

  res_filt <- subset(res, res$padj <= 0.05 & abs(res$LFC) >= 1)

  list(
    DGE_list = dgeList,
    res = res_filt
  )
}
