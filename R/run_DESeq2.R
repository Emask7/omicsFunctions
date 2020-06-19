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
      countData = geneCounts, colData = sampleData,
      design = ~ animal + timepoint
    )

  # Wald test -----------------------------------------------------------------
    waldDDS <- DESeq(dds, quiet = TRUE)
    wald <- DESeq2::results(waldDDS, lfcThreshold = 1)

    wald_res <- data.frame(rownames(wald), wald[, c(2, 6)])
    colnames(wald_res) <- c("Gene_ID", "LFC", "padj")
    wald_res <- convert_IDs(wald_res, orthoList)

    print("DESeq2 Wald test DEG summary:", quote = FALSE)
    DESeq2::summary(wald, 0.05)
    print(" ", quote = FALSE)

  # Wald test with LFC shrinkage ----------------------------------------------
    shrink <- lfcShrink(
      waldDDS, type = "apeglm", lfcThreshold = 1, coef = (4), quiet = TRUE
    )

    shrink_res <- data.frame(rownames(shrink), shrink[, c(2, 4)])
    colnames(shrink_res) <- c("Gene_ID", "LFC", "padj")
    shrink_res <- convert_IDs(shrink_res, orthoList)

    print("DESeq2 Wald test with LFC shrinkage DEG summary:", quote = FALSE)
    DESeq2::summary(shrink, 0.05)

  # Return a list object with all of the data ---------------------------------
    list(Wald = wald_res, LFCshrinkage = shrink_res)
}
