#' @title DESeq2 wrapper
#'
#' @description DESeq2 wrapper
#' @param counts_data A data frame containing 6 columns of raw gene or transcript counts.
#' @param count_type string. Either "genes" or "transcripts"
#' @import DESeq2
#' @export
#' @examples
#' run_DESeq2_humanOnly(counts_data)

run_DESeq2_humanOnly <- function(counts_data, count_type) {
  # Specify experimental design factors ---------------------------------------
    animal <- factor(c(rep(c("15979", "30760", "31151"), 2)))
    timepoint <- factor(c(rep("T-7", 3), rep("T15", 3)))
    sampleData <- data.frame(animal, timepoint)

  # Make a DESeqDataSet object ------------------------------------------------
    counts_data <- as.matrix(counts_data)
    counts_data[, 1:6] <- as.integer(counts_data)

    dds <- DESeqDataSetFromMatrix(
      countData = counts_data, colData = sampleData,
      design = ~ animal + timepoint
    )

  # Wald test with LFC shrinkage ----------------------------------------------
    waldDDS <- DESeq(dds, quiet = TRUE)
    shrink <- lfcShrink(
      waldDDS, type = "apeglm", lfcThreshold = 1, coef = (4), quiet = TRUE
    )

    res <- data.frame(rownames(shrink), shrink[, c(2, 4)])
    if (count_type == "genes") colnames(res) <- c("Gene_ID", "LFC", "padj")
    else if (count_type == "transcripts") {
      colnames(res) <- c("Transcript_ID", "LFC", "padj")
    }

    res_filtered <- subset(res, res$padj <= 0.05 & abs(res$LFC) >= 1)

    print("DESeq2 Wald test with LFC shrinkage DEG summary:", quote = FALSE)
    DESeq2::summary(shrink, 0.05)

  # Return a list object with all of the data ---------------------------------
    list(results = res, filtered_results = res_filtered)
}
