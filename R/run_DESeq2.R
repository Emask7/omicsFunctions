#' @title DESeq2 wrapper
#'
#' @description DESeq2 wrapper
#' @param geneCounts A data frame containing 6 columns of raw gene counts.
#' @param orthoList A data frame containing two columns (Gene_ID and Human_Gene_ID)
#' @import DESeq2
#' @export
#' @examples
#' run_DESeq2(geneCounts, orthoList)

run_DESeq2 <- function(geneCounts, orthoList) {
  # Specify experimental design factors ---------------------------------------
    animal <- factor(c(rep(c("15979", "30760", "31151"), 2)))
    timepoint <- factor(c(rep("T-7", 3), rep("T15", 3)))
    sampleData <- data.frame(animal, timepoint)


  # Make a DESeqDataSet object ------------------------------------------------
    geneCounts <- as.matrix(geneCounts)
    geneCounts[, 1:6] <- as.integer(geneCounts)

    dds <- DESeqDataSetFromMatrix(
      countData = geneCounts,
      colData = sampleData,
      design = ~ animal + timepoint
    )

  # Wald test -----------------------------------------------------------------
    waldDDS <- DESeq(dds, quiet = TRUE)
    waldRes <- results(waldDDS, lfcThreshold = 1)
    waldSummary <- DESeq2::summary(waldRes, alpha = 0.05)
    waldRes <- data.frame(rownames(waldRes), waldRes[, c(2, 6)])
    colnames(waldRes) <- c("Gene_ID", "LFC", "padj")
    waldRes <- convert_IDs(waldRes, orthoList)

  # Wald test with LFC shrinkage ----------------------------------------------
    shrinkRes <- lfcShrink(
      waldDDS, type = "apeglm", lfcThreshold = 1, coef = (4), quiet = TRUE
    )
    shrinkSummary <- DESeq2::summary(shrinkRes, alpha = 0.05)
    shrinkRes <- data.frame(rownames(shrinkRes), shrinkRes[, c(2, 4)])
    colnames(shrinkRes) <- c("Gene_ID", "LFC", "padj")
    shrinkRes <- convert_IDs(shrinkRes, orthoList)


  # Return a list object with all of the data ---------------------------------
    list(
      Wald = waldRes, Wald_summary = waldSummary,
      LFCshrinkage = shrinkRes, LFCshrinkage_summary = shrinkSummary
    )
}
