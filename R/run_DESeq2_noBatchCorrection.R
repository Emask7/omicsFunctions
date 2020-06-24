#' @title DESeq2 wrapper
#'
#' @description DESeq2 wrapper
#' @param geneCounts A data frame containing 6 columns of raw gene counts.
#' @param orthoList A data frame containing two columns (Gene_ID and Human_Gene_ID)
#' @param LFC_filter logical. If TRUE, DEG lists are filtered by LFC values in addition to FDR values
#' @import DESeq2
#' @export
#' @examples
#' run_DESeq2_noBatchCorrection(geneCounts, orthoList, LFC_filter)

run_DESeq2_noBatchCorrection <- function(geneCounts, orthoList, LFC_filter) {
  # Specify experimental design factors ---------------------------------------
  timepoint <- factor(c(rep("T-7", 3), rep("T15", 3)))
  sampleData <- data.frame(timepoint)

  # Make a DESeqDataSet object ------------------------------------------------
  geneCounts <- as.matrix(geneCounts)
  geneCounts[, 1:6] <- as.integer(geneCounts)

  dds <- DESeqDataSetFromMatrix(
    countData = geneCounts, colData = sampleData,
    design = ~ timepoint
  )

  # Wald test -----------------------------------------------------------------
  waldDDS <- DESeq(dds, quiet = TRUE)
  # wald <- DESeq2::results(waldDDS, lfcThreshold = 1)
  #
  # wald_res <- data.frame(rownames(wald), wald[, c(2, 6)])
  # colnames(wald_res) <- c("Gene_ID", "LFC", "padj")
  # wald_res <- filter_DEG_table(wald_res, orthoList, LFC_filter)
  #
  # print("DESeq2 Wald test DEG summary:", quote = FALSE)
  # DESeq2::summary(wald, 0.05)
  # print(" ", quote = FALSE)

  # Wald test with LFC shrinkage ----------------------------------------------
  shrink <- lfcShrink(
    waldDDS, type = "apeglm", lfcThreshold = 1, coef = (2), quiet = TRUE
  )

  shrink_res <- data.frame(rownames(shrink), shrink[, c(2, 4)])
  colnames(shrink_res) <- c("Gene_ID", "LFC", "padj")
  shrink_res <- filter_DEG_table(shrink_res, orthoList, LFC_filter)

  print("DESeq2 Wald test with LFC shrinkage DEG summary:", quote = FALSE)
  DESeq2::summary(shrink, 0.05)

  # Return a list object with all of the data ---------------------------------
  list(res = shrink_res, dds = waldDDS)
}
