#' @title Filter differential expression data
#'
#' @description Filters differential expression data to remove all genes with an FDR value of > 0.05, log-fold-change values between -1 and 1, and/or genes with only ENSPANG IDs
#' @param data_table data.fame of DEGs
#' @param human_IDs data.frame with human gene names
#' @param LFC_filter logical. If TRUE, DEG lists are filtered by LFC values in addition to FDR values
#' @export
#' @examples
#' filter_DEG_table(data_table, human_IDs, LFC_filter)

filter_DEG_table <- function(data_table, human_IDs, LFC_filter) {
  if (LFC_filter) {
    filtered <- subset(
      data_table, data_table$padj <= 0.05 & abs(data_table$LFC) >= 1
    )
  } else filtered <- subset(data_table, data_table$padj <= 0.05)

  if (nrow(filtered) < 1) {
    print("No significant DEGs", quote = FALSE)
    return(data.frame())
  }
  else {
    filtered <- dplyr::left_join(filtered, human_IDs, by = "Gene_ID")
    for (n in 1:nrow(filtered)) {
      if (!is.na(filtered[n, 4]) & grepl("ENSPANG", filtered[n, 1])) {
        filtered[n, 1] <- filtered[n, 4]
      }
    }
    return(filtered[, 1:3])
  }
}

