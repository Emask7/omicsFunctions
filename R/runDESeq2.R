#' @title DESeq2 wrapper
#'
#' @description
#' @param geneCounts A data frame containing 6 columns of raw gene counts.
#' @param orthoList A data frame containing two columns (Gene_ID and Human_Gene_ID)
#' @export
#' @examples
#' runDESeq2()

runDESeq2 <- function(geneCounts, orthoList) {
  # Specify experimental design factors ---------------------------------------
    animal <- factor(c(rep(c("15979", "30760", "31151"), 2)))
    timepoint <- factor(c(rep("T-7", 3), rep("T15", 3)))
    sampleData <- data.frame(animal, timepoint)


  # Make a DESeqDataSet object ------------------------------------------------
    geneCounts <- as.matrix(geneCounts)
    geneCounts[, 1:6] <- as.integer(geneCounts)

    dds <- DESeqDataSetFromMatrix(
      countData = geneCounts, colData = sampleData, design = ~animal+timepoint
    )

  # Wald test -----------------------------------------------------------------
    waldDDS <- DESeq(dds)
    waldRes <- results(waldDDS, lfcThreshold = 1)
    # waldSummary <- summary(waldRes, alpha = 0.05)
    waldRes <- data.frame(rownames(waldRes), waldRes[, c(2, 6)])
    colnames(waldRes) <- c("Gene_ID", "LFC", "padj")
    waldRes <- convertIDs(waldRes, orthoList)

  # Wald test with LFC shrinkage ----------------------------------------------
    shrinkRes <- lfcShrink(
      waldDDS, type = "apeglm", lfcThreshold = 1, coef = (4)
    )
    # shrinkSummary <- summary(shrinkRes, alpha = 0.05)
    shrinkRes <- data.frame(rownames(shrinkRes), shrinkRes[, c(2, 4)])
    colnames(shrinkRes) <- c("Gene_ID", "LFC", "padj")
    shrinkRes <- convertIDs(shrinkRes, orthoList)


  # Return a list object with all of the data ---------------------------------
    filt <- function(x) subset(x, abs(x$LFC) >= 1 & x$padj <= 0.05)

    list(
      # Wald_summary = waldSummary, LFCshrinkage_summary = shrinkSummary,
      Wald = waldRes, Wald_FDRfiltered = filt(waldRes),
      LFCshrinkage = shrinkRes, LFCshrinkage_FDRfiltered = filt(shrinkRes)
    )

}
