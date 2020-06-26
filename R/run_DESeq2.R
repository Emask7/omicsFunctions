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
