#' @title ENSPANG gene ID to human gene name converter
#'
#' @description
#' @param dataTable
#' @param humanIDs
#' @export
#' @examples
#' convertIDs()

convertIDs <- function(dataTable, humanIDs) {
  if(nrow(dataTable) > 0) {
    res <- left_join(dataTable, humanIDs)
    for(n in 1:nrow(res)) {
      if(!is.na(res[n, 4]) & grepl("ENSPANG", res[n, 1])) {
        res[n, 1] <- res[n, 4]
      }
    }
    res[, 1:3]
  } else print("Error: the input table has zero rows")
}
