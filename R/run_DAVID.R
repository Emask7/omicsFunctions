#' @title DAVID Gene Ontology functional annotation wrapper
#'
#' @description DAVID Gene Ontology functional annotation wrapper
#' @param davidWS A DAVIDWebService class object.
#' @param annotation_cat a string or vector of strings indicating which set of GO annotation terms to use
#' @param split_at a character at which strings will be split. If NULL (the default), strings will not be split.
#' @param DE_data A data frame containing 3 columns ('Gene_ID' or 'Transcript_ID', 'LFC', and 'padj').
#' @import limma
#' @import stringi
#' @import RDAVIDWebService
#' @examples
#' run_FA(davidWS, annotation_cat, split_at, )

run_FA <- function(davidWS, annotation_cat, split_at, DE_data) {
  setAnnotationCategories(davidWS, c(annotation_cat))
  FA_chart <- getFunctionalAnnotationChart(davidWS)
  FA_sub <- subset(FA_chart, FA_chart$FDR <= 0.05)

  if (nrow(FA_sub) < 1) return(data.frame())

  gene_list <- data.frame(FA_sub$Genes)
  for (r in 1:nrow(gene_list)) {
    gene_names <- limma::strsplit2(gene_list[r, 1], ", ")[1, ]
    gene_names <- list(
      AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db, gene_names,
        keytype = "ENSEMBL", column = "SYMBOL"
      )
    )
    gene_names <- stri_join_list(gene_names, sep = ", ")
    gene_list[r, 1] <- gene_names
  }

  if (is.null(split_at)) {
    FA_sub <- data.frame(
      FA_sub$Category, FA_sub$Term, FA_sub$Count, FA_sub$X.,
      gene_list, FA_sub$Fold.Enrichment, FA_sub$FDR
    )
    colnames(FA_sub) <- c(
      "Category", "Term", "Count", "Percent_of_input_list",
      "Genes", "Fold.Enrichment", "FDR"
    )
    return(FA_sub)
  } else {
    split_IDs <- limma::strsplit2(FA_sub$Term, split_at)
    FA_sub <- data.frame(
      FA_sub$Category, split_IDs, FA_sub$Count, FA_sub$X.,
      gene_list, FA_sub$Fold.Enrichment, FA_sub$FDR
    )
    colnames(FA_sub) <- c(
      "Category", "ID", "Term", "Count", "Percent_of_input_list",
      "Genes", "Fold.Enrichment", "FDR"
    )

    terms_data <- data.frame(FA_sub[, c(1:3, 8, 6)])
    colnames(terms_data) <- c("category", "ID", "term", "adj_pval", "genes")

    genes_data <- data.frame(DE_data$Gene_ID, DE_data$LFC)
    colnames(genes_data) <- c("ID", "logFC")

    cd <- GOplot::circle_dat(terms_data, genes_data)

    z_scores <- cd[, c(2, 8)]
    z_scores <- z_scores[!duplicated(z_scores), ]
    z_scores <- dplyr::left_join(FA_sub, z_scores, by = "ID")
    z_scores <- subset(z_scores, abs(z_scores$zscore) >= 2)

    return(list(FA_results = z_scores, GOplot_results = cd))
  }
}



#' @title DAVID Gene Ontology Wrapper
#'
#' @description DAVID Gene Ontology Wrapper
#' @param davidWS A DAVIDWebService class object.
#' @param DE_data A data frame containing 3 columns ('Gene_ID' or 'Transcript_ID', 'LFC', and 'padj').
#' @param list_name A string that will be used to identify the submitted list.
#' @param list_type string. Either 'gene_symbol', 'gene_ID', or 'transcript_ID'
#' @param overwrite_file logical
#' @import limma
#' @import stringi
#' @import RDAVIDWebService
#' @import openxlsx
#' @import GOplot
#' @export
#' @examples
#' run_DAVID(davidWS, DE_data, list_name, list_type, overwrite_file)

run_DAVID <- function(davidWS, DE_data, list_name, list_type, overwrite_file) {
  # Make sure input is valid --------------------------------------------------
    # DE_data <- subset.data.frame(DE_data, !grepl("ENSPANG", DE_data$Gene_ID))

    if (nrow(DE_data) < 1) {
      print("Error: input has 0 rows")
      return(NULL)
    }
    if (colnames(DE_data)[2] != "LFC" | colnames(DE_data)[3] != "padj") {
      print("Error: input must have 3 columns (Gene_ID or Transcript_ID, LFC, and padj)")
      return(NULL)
    }
    if (!is.connected(david)) {
      print("Error: not connected to RDAVIDWebService")
      return(NULL)
    }

  # Check to see if gene list has been added ----------------------------------
    gene_lists <- getGeneListNames(davidWS)
    list_position <- 0

    if (length(getGeneListNames(davidWS)) > 0) {
      for (x in 1:length(gene_lists)) {
        if (gene_lists[x] == list_name) list_position <- x
      }
    }

    if (list_position == 0) {
      if (list_type == "gene_symbol") {
        input_IDs <- AnnotationDbi::mapIds(
          org.Hs.eg.db::org.Hs.eg.db, DE_data[, 1],
          keytype = "SYMBOL", column = "ENSEMBL"
        )
        addList(
          davidWS, input_IDs, idType = "ENSEMBL_GENE_ID",
          listName = list_name, listType = "Gene"
        )
      } else if (list_type == "gene_ID") {
        input_IDs <- DE_data[, 1]
        addList(
          davidWS, input_IDs, idType = "ENSEMBL_GENE_ID",
          listName = list_name, listType = "Gene"
        )
      } else if (list_type == "transcript_ID") {
        input_IDs <- DE_data[, 1]
        addList(
          davidWS, input_IDs, idType = "ENSEMBL_TRANSCRIPT_ID",
          listName = list_name, listType = "Gene"
        )
      }
      gene_lists <- getGeneListNames(davidWS)
    } else setCurrentGeneListPosition(davidWS, list_position)

    print("Running GO analysis on gene list:")
    print(gene_lists[getCurrentGeneListPosition(davidWS)])

  # Run GO analysis -----------------------------------------------------------
    # setAnnotationCategories(davidWS, c("GOTERM_BP_ALL"))
    bp_res <- run_FA(davidWS, c("GOTERM_BP_ALL"), "~", DE_data)

    if(nrow(bp_res$FA_results) > 0) {
      file_name <- stri_join(list_name, " GO Terms.xlsx")
      wb <- createWorkbook(file_name)
      addWorksheet(wb, "Bio_Process")
      writeData(wb, "Bio_Process", bp_res$FA_results)
      saveWorkbook(wb, file_name, overwrite = overwrite_file)

      return(
        list(
          DAVID_results = bp_res$FA_results,
          GOplot_results = bp_res$GOplot_results
        )
      )
    } else return(data.frame())
}
