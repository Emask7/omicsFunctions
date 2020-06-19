#' @title edgeR wrapper
#'
#' @description edgeR wrapper
#' @param gene_counts A data frame containing 6 columns of raw gene counts.
#' @param ortho_list A data frame containing two columns (Gene_ID and Human_Gene_ID)
#' @param min_count numberic. The minimum number of reads required for at least some samples. By default, min_count is set to 5.
#' @param norm_method The normalization method to be used for analysis. Options are "TMM", "TMMwsp", "RLE", "UQ", or "none". By default, this is set to "TMM".
#' @param show_plots logical. If TRUE, plots are shown. If FALSE, plots are not shown. By default, this is set to FALSE.
#' @import edgeR
#' @export
#' @examples
#' run_edgeR(gene_counts, ortho_list)

run_edgeR <- function(gene_counts,
                     ortho_list,
                     min_count = 5,
                     norm_method = "TMM",
                     show_plots) {

  # Specify experimental design factors ---------------------------------------
    animal <- factor(c(rep(c("15979", "30760", "31151"), 2)))
    timepoint <- factor(c(rep("T-7", 3), rep("T15", 3)))
    design <- model.matrix(~ animal + timepoint)

  # Make a DGEList object from the raw count data and filter for low reads ----
    dgeList <- DGEList(counts = gene_counts, group = timepoint)
    keep <- filterByExpr(dgeList, min.count = min_count)
    dgeList <- dgeList[keep, , keep.lib.sizes = FALSE]

  # Calculate normalization factors -------------------------------------------
    dgeList <- calcNormFactors(dgeList, method = norm_method)

  # Data exploration: multidimensional scaling plot ---------------------------
    if (show_plots) limma::plotMDS(dgeList, main = "MDS Plot")

  # Estimate dispersion -------------------------------------------------------
    dgeList <- estimateDisp(dgeList, design, robust = TRUE)

  # Data exploration: biological coefficient of variation plot ----------------
    if (show_plots) plotBCV(dgeList, main = "BCV Plot")

  # DE testing: GLM approach --------------------------------------------------
    fit <- glmFit(dgeList, design)
    lrt <- glmLRT(fit)
    lrtRes <- data.frame(
      rownames(lrt$table), lrt$table$logFC,
      p.adjust(lrt$table$PValue, method = "fdr")
    )
    colnames(lrtRes) <- c("Gene_ID", "LFC", "padj")
    lrtRes <- convert_IDs(lrtRes, ortho_list)
    lrtSummary <- summary(
      limma::decideTests(lrt, adjust.method = "fdr", lfc = 1)
    )

  # DE above a FC threshold: GLM approach -------------------------------------
    lrtTreat <- glmTreat(fit, lfc = 1, null = "interval")
    lrtTreatRes <- data.frame(
      rownames(lrtTreat$table), lrtTreat$table$logFC,
      p.adjust(lrtTreat$table$PValue, method = "fdr")
    )
    colnames(lrtTreatRes) <- c("Gene_ID", "LFC", "padj")
    lrtTreatRes <- convert_IDs(lrtTreatRes, ortho_list)
    lrtTreatSummary <- summary(
      limma::decideTests(lrtTreat, adjust.method = "fdr", lfc = 1)
    )

  # Mean Difference Plots -----------------------------------------------------
    if (show_plots) {
      limma::plotMD(lrt, main = "Mean-Difference Plot - LRT")
      limma::plotMD(lrtTreat, main = "Mean-Difference Plot - LRT with TREAT")
    }

  # Combine all summaries -----------------------------------------------------
    summary_all <- data.frame(
      lrtSummary[, 1], lrtTreatSummary[, 1],
      row.names = c("Down", "NotSig", "Up")
    )
    colnames(summary_all) <- c("LRT", "LRT + TREAT")

  # Return a list object with all of the data ---------------------------------
    filt <- function(x) subset(x, abs(x$LFC) >= 1 & x$padj <= 0.05)

    list(
      summary = summary_all,
      LRT = lrtRes, LRT_FDRfiltered = filt(lrtRes),
      LRT_TREAT = lrtTreatRes, LRT_TREAT_FDRfiltered = filt(lrtTreatRes)
    )
}
