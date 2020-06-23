#' @title Write DEG list to Excel file
#'
#' @description Write DEG list to Excel file
#' @param CD4T_res data.frame with columns 'Gene_ID', 'LFC', and 'padj'
#' @param CD8T_res data.frame with columns 'Gene_ID', 'LFC', and 'padj'
#' @param NK_res data.frame with columns 'Gene_ID', 'LFC', and 'padj'
#' @param file_name string
#' @export
#' @examples
#' write_DEGs_to_Excel(CD4T_res, CD8T_res, NK_res, file_name)

write_DEGs_to_Excel <- function(CD4T_res, CD8T_res, NK_res, file_name) {
  temp_CD4T <- CD4T_res
  colnames(temp_CD4T) <- c("Gene_ID", "LFC_CD4T", "padj_CD4T")

  temp_CD8T <- CD8T_res
  colnames(temp_CD8T) <- c("Gene_ID", "LFC_CD8T", "padj_CD8T")

  all_cell_res <- dplyr::full_join(temp_CD4T, temp_CD8T)

  if (nrow(NK_res) > 0) {
    temp_NK <- NK_res
    colnames(temp_NK) <- c("Gene_ID", "LFC_NK", "padj_NK")

    all_cell_res <- dplyr::full_join(all_cell_res, temp_NK)
  }

  wb <- createWorkbook(file_name)

  addWorksheet(wb, "All Cell Types")
  writeData(wb, "All Cell Types", all_cell_res)

  addWorksheet(wb, "CD4T")
  writeData(wb, "CD4T", CD4T_res)

  addWorksheet(wb, "CD8T")
  writeData(wb, "CD8T", CD8T_res)

  addWorksheet(wb, "NK")
  writeData(wb, "NK", NK_res)

  saveWorkbook(wb, file_name, overwrite = TRUE)
}



#' @title Write list of DEGs mapped to a given GO term to an Excel file
#'
#' @description Write list of DEGs mapped to a given GO term to an Excel file
#' @param GO_term a string containing the Gene Ontology term of interest
#' @param GO_dat data.frame with Gene Ontology analysis results to be searched for the term specified by GO_term.
#' Must include columns for at least 'term' and 'genes'.
#' @param deg_dat data.frame with columns 'Gene_ID', 'LFC', and 'padj'
#' @param file_name string
#' @param sheet_name string
#' @param overwrite_file logical
#' @export
#' @examples
#' get_term_DEGs(GO_term, GO_dat, deg_dat, file_name, sheet_name, overwrite_file)

get_term_DEGs <- function(GO_term, GO_dat, deg_dat, file_name, sheet_name, overwrite_file) {
  res <- subset(GO_dat, GO_dat$term == GO_term)
  res <- data.frame(res$term, res$genes)
  colnames(res) <- c("Term", "Gene_ID")
  res <- dplyr::left_join(res, deg_dat)

  wb <- createWorkbook(file_name)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, res)
  saveWorkbook(wb, file_name, overwrite = overwrite_file)

  res
}
