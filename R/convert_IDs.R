#' @title ENSPANG gene ID to human gene name converter
#'
#' @description converts ENSPANG gene IDs to human gene names
#' @param data_table data.fame of DEGs
#' @param human_IDs data.frame with human gene names
#' @import dplyr
#' @export
#' @examples
#' convert_IDs()

convert_IDs <- function(data_table, human_IDs) {
  if (nrow(data_table) <= 0) print("Error: the input table has zero rows")
  else {
    res <- dplyr::left_join(data_table, human_IDs)
    for (n in 1:nrow(res)) {
      if (!is.na(res[n, 4]) & grepl("ENSPANG", res[n, 1])) {
        res[n, 1] <- res[n, 4]
      }
    }
    res[, 1:3]
  }
}
