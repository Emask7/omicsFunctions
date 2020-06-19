#' @title edgeR wrapper
#'
#' @description edgeR wrapper
#' @param gene_counts A data frame containing 6 columns of raw gene counts.
#' @param ortho_list A data frame containing two columns (Gene_ID and Human_Gene_ID)
#' @param show_plots logical. If TRUE, plots are shown. If FALSE, plots are not shown.
#' @import edgeR
#' @export
#' @examples
#' run_edgeR(gene_counts, ortho_list)

run_edgeR <- function(gene_counts, ortho_list, show_plots) {

  # Specify experimental design factors ---------------------------------------
    animal <- factor(c(rep(c("15979", "30760", "31151"), 2)))
    timepoint <- factor(c(rep("T-7", 3), rep("T15", 3)))
    design <- model.matrix(~ animal + timepoint)

  # Make a DGEList object from the raw count data and filter for low reads ----
    dgeList <- DGEList(counts = gene_counts, group = timepoint)
    keep <- filterByExpr(dgeList, min.count = 5)
    dgeList <- dgeList[keep, , keep.lib.sizes = FALSE]

  # Calculate normalization factors -------------------------------------------
    dgeList <- calcNormFactors(dgeList, method = "TMM")

  # Data exploration: multidimensional scaling plot ---------------------------
    if (show_plots) limma::plotMDS(dgeList, main = "MDS Plot")

  # Estimate dispersion -------------------------------------------------------
    dgeList <- estimateDisp(dgeList, design, robust = TRUE)

  # Data exploration: biological coefficient of variation plot ----------------
    if (show_plots) plotBCV(dgeList, main = "BCV Plot")

  # DE testing ----------------------------------------------------------------
    fit <- glmFit(dgeList, design)
    lrtTreat <- glmTreat(fit, lfc = 1, null = "interval")

    res <- data.frame(
      rownames(lrtTreat$table), lrtTreat$table$logFC,
      p.adjust(lrtTreat$table$PValue, method = "fdr")
    )
    colnames(res) <- c("Gene_ID", "LFC", "padj")
    res <- filter_DEG_table(res, ortho_list)

  # Mean Difference Plots -----------------------------------------------------
    if (show_plots) limma::plotMD(lrtTreat, main = "Mean-Difference Plot")

  # Return a list object with all of the data ---------------------------------
    DE_summary <- base::summary(
      limma::decideTests(lrtTreat, adjust.method = "fdr", p.value = 0.05)
    )
    list(results = res, summary = DE_summary)
}
