#' @title edgeR wrapper
#'
#' @description
#' @param geneCounts A data frame containing 6 columns of raw gene counts.
#' @param minCount numberic. The minimum number of reads required for at least some samples. By default, minCount is set to 5.
#' @param normMethod The normalization method to be used for analysis. Options are "TMM", "TMMwsp", "RLE", "UQ", or "none". By default, this is set to "TMM".
#' @param showPlots logical. If TRUE, plots are shown. If FALSE, plots are not shown. By default, this is set to FALSE.
#' @export
#' @examples
#' runEdgeR()

runEdgeR <- function(geneCounts, minCount = 5, normMethod = "TMM", showPlots = FALSE) {
  # Specify experimental design factors ---------------------------------------
    animal <- factor(c(rep(c("15979", "30760", "31151"), 2)))
    timepoint <- factor(c(rep("T-7", 3), rep("T15", 3)))
    design <- model.matrix(~animal+timepoint)

  # Make a DGEList object from the raw count data and filter for low reads ----
    dgeList <- DGEList(counts = geneCounts, group = timepoint)
    keep <- filterByExpr(dgeList, min.count = minCount)
    dgeList <- dgeList[keep, , keep.lib.sizes = FALSE]

  # Calculate normalization factors -------------------------------------------
    dgeList <- calcNormFactors(dgeList, method = normMethod)

  # Data exploration: multidimensional scaling plot ---------------------------
    if(showPlots) plotMDS(dgeList, main = "MDS Plot")

  # Estimate dispersion -------------------------------------------------------
    dgeList <- estimateDisp(dgeList, design, robust = TRUE)

  # Data exploration: biological coefficient of variation plot ----------------
    if(showPlots) plotBCV(dgeList, main = "BCV Plot")

  # DE testing: GLM approach --------------------------------------------------
    lrtFit <- glmFit(dgeList, design)
    lrt <- glmLRT(lrtFit)
    lrtRes <- data.frame(
      rownames(lrt$table), lrt$table$logFC,
      p.adjust(lrt$table$PValue, method = "fdr")
    )
    colnames(lrtRes) <- c("Gene_ID", "LFC", "FDR")
    lrtSummary <- summary(decideTests(lrt, adjust.method = "fdr", lfc = 1))

  # # DE testing: QL GLM approach -----------------------------------------------
  #   qlFit <- glmQLFit(dgeList, design, robust = TRUE)
  #   qlf <- glmQLFTest(qlFit)
  #   qlfRes <- data.frame(
  #     rownames(qlf$table), qlf$table$logFC,
  #     p.adjust(qlf$table$PValue, method = "fdr")
  #   )
  #   colnames(qlfRes) <- c("Gene_ID", "LFC", "FDR")
  #   qlfSummary <- summary(decideTests(qlf, adjust.method = "fdr", lfc = 1))
  #
  #   if(showPlots) plotQLDisp(qlFit, main = "QL Dispersion")

  # DE above a FC threshold: GLM approach -------------------------------------
    lrtTreat <- glmTreat(lrtFit, lfc = 1, null = "interval")
    lrtTreatRes <- data.frame(
      rownames(lrtTreat$table), lrtTreat$table$logFC,
      p.adjust(lrtTreat$table$PValue, method = "fdr")
    )
    colnames(lrtTreatRes) <- c("Gene_ID", "LFC", "FDR")
    lrtTreatSummary <- summary(
      decideTests(lrtTreat, adjust.method = "fdr", lfc = 1)
    )

  # # DE above a FC threshold: QL GLM approach ----------------------------------
  #   qlfTreat <- glmTreat(qlFit, lfc = 1, null = "interval")
  #   qlfTreatRes <- data.frame(
  #     rownames(qlfTreat$table), qlfTreat$table$logFC,
  #     p.adjust(qlfTreat$table$PValue, method = "fdr")
  #   )
  #   colnames(qlfTreatRes) <- c("Gene_ID", "LFC", "FDR")
  #   qlfTreatSummary <-  summary(
  #     decideTests(qlfTreat, adjust.method = "fdr", lfc = 1)
  #   )

  # Mean Difference Plots -----------------------------------------------------
    if(showPlots) {
      plotMD(lrt, main = "Mean-Difference Plot - LRT")
      plotMD(lrtTreat, main = "Mean-Difference Plot - LRT with TREAT")
      # plotMD(qlf, main = "Difference Plot - QLF")
      # plotMD(qlfTreat, main = "Difference Plot - QLF with TREAT")
    }

  # Combine all summaries -----------------------------------------------------
    summary_all <- data.frame(
      lrtSummary[, 1], lrtTreatSummary[, 1],
      # qlfSummary[, 1], qlfTreatSummary[, 1],
      row.names = c("Down", "NotSig", "Up")
    )
    colnames(summary_all) <- c("LRT", "LRT + TREAT")
    # colnames(summary_all) <- c("LRT", "LRT + TREAT", "QLF", "QLF + TREAT")

  # Return a list object with all of the data ---------------------------------
    filt <- function(x) subset(x, abs(x$LFC) >= 1 & x$FDR <= 0.05)

    list(
      summary = summary_all,
      LRT = lrtRes, LRT_FDRfiltered = filt(lrtRes),
      LRT_TREAT = lrtTreatRes, LRT_TREAT_FDRfiltered = filt(lrtTreatRes)
      # results = list(
      #   LRT = lrtRes, LRT_TREAT = lrtTreatRes
      #   # QLF = qlfRes, QLF_TREAT = qlfTreatRes
      # ),
      # FDRfilteredResults = list(
      #   LRT = filt(lrtRes), LRT_TREAT = filt(lrtTreatRes)
      #   # QLF = filt(qlfRes), QLF_TREAT = filt(qlfTreatRes)
      # ),
      # summaries = list(
      #   all = summary_all,
      #   LRT = lrtSummary, LRT_TREAT = lrtTreatSummary,
      #   QLF = qlfSummary, QLF_TREAT = qlfTreatSummary
      # )
    )
}
