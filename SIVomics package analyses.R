roxygen2::roxygenise()


library(SIVomics)
library(openxlsx)
library(RDAVIDWebService)



# Import data -----------------------------------------------------------------
  sample_list <- c(
    "CD4T_T-7_15979", "CD4T_T-7_30760", "CD4T_T-7_31151",
    "CD4T_T15_15979", "CD4T_T15_30760", "CD4T_T15_31151",
    "CD8T_T-7_15979", "CD8T_T-7_30760", "CD8T_T-7_31151",
    "CD8T_T15_15979", "CD8T_T15_30760", "CD8T_T15_31151",
    "NK_T-7_15979", "NK_T-7_30760", "NK_T-7_31151",
    "NK_T15_15979", "NK_T15_30760", "NK_T15_31151"
  )

  gCounts <- read.xlsx("copy - Partek_LG_RNA_Seq_20190304_raw_human_gene_counts.xlsx")
  rownames(gCounts) <- make.names(gCounts[, 5], unique = TRUE)
  gCounts <- gCounts[, 7:24]
  colnames(gCounts) <- c(sample_list)
  head(gCounts)

  tCounts <- read.xlsx("copy - Partek_LG_RNA_Seq_20190304_raw_human_transcript_counts.xlsx")
  rownames(tCounts) <- make.names(tCounts[, 16], unique = TRUE)
  tCounts <- tCounts[, c(5:6, 9, 20:37)]
  colnames(tCounts) <- c("Transcript", "Gene_ID", "Ensembl_Gene_ID", sample_list)
  head(tCounts)

  animal <- factor(c(rep(c("15979", "30760", "31151"), 2)))
  timepoint <- factor(c(rep("T-7", 3), rep("T15", 3)))
  sampleData <- data.frame(animal, timepoint)
  sampleData

  design <- ~ animal + timepoint

# Differential expression analyses --------------------------------------------
  DE_genes <- list(
    CD4T = run_DESeq2(gCounts[, 1:6], sampleData, design, c(4), "gene"),
    CD8T = run_DESeq2(gCounts[, 7:12], sampleData, design, c(4), "gene"),
    NK = run_DESeq2(gCounts[, 13:18], sampleData, design, c(4), "gene")
  )

  write_DEGs_to_Excel(
    DE_genes$CD4T$filtered_results,
    DE_genes$CD8T$filtered_results,
    DE_genes$NK$filtered_results,
    "DEGs - DESeq2 LFC Shrinkage Method - human alignment.xlsx"
  )


  DE_transcripts <- list(
    CD4T = run_DESeq2(tCounts[, 4:9], sampleData, design, c(4), "transcript"),
    CD8T = run_DESeq2(tCounts[, 10:15], sampleData, design, c(4), "transcript"),
    NK = run_DESeq2(tCounts[, 16:21], sampleData, design, c(4), "transcript")
  )

  transcript_info <- data.frame(rownames(tCounts), tCounts[, 1:3])
  colnames(transcript_info) <- c(
    "Ensembl_Transcript_ID" "Transcript_ID", "Gene_ID", "Ensembl_Gene_ID"
  )
  head(transcript_info)

  write_DETs_to_Excel(
    DE_transcripts$CD4T$filtered_results,
    DE_transcripts$CD8T$filtered_results,
    DE_transcripts$NK$filtered_results,
    transcript_info,
    "DE Transcripts - DESeq2 LFC Shrinkage Method - human alignment.xlsx"
  )



# DAVID Gene Ontology analyses ------------------------------------------------
  david <- DAVIDWebService(
    email = "emask@txbiomed.org",
    url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/"
  )

  is.connected(david)
  show(david)

  getIdTypes(david)
  #' run_DAVID(davidWS, DE_data, list_name, list_type)

  CD4T_GO <- list(
    genes = run_DAVID(
      david, DE_genes$CD4T$filtered_results,
      "human_gene_counts", "gene_symbol", TRUE
    ),
    transcripts = run_DAVID(
      david, DE_transcripts$CD4T$filtered_results,
      "human_transcript_counts", "transcript_ID", TRUE
    )
  )

  test <- run_DAVID(
    david, DE_genes$CD4T$filtered_results,
    "human_gene_counts", "gene_symbol", TRUE
  )

  run_DAVID(david, CD4T$human_DESeq2$LFCshrinkage, "CD4T_human_LFCshrink")
  CD8T_GO <- run_DAVID(david, CD8T$human_DESeq2$LFCshrinkage, "CD8T_human_LFCshrink")
  NK_GO <- run_DAVID(david, NK$human_DESeq2$LFCshrinkage, "NK_human_LFCshrink")


  # Troubleshooting the "Read timed out" error --------------------------------
    setTimeOut(david, 50000)
    getTimeOut(david)

    setHttpProtocolVersion(david, "HTTP/1.0")
    getHttpProtocolVersion(david)





#
#
#
#
#
# CD4T_temp <- GO_results$CD4T_human$LFCshrink[, c(3, 9)]
# colnames(CD4T_temp) <- c("Term", "Zscore_CD4T")
#
# CD8T_temp <- GO_results$CD8T_baboon$LFCshrink[, c(3, 9)]
# colnames(CD8T_temp) <- c("Term", "Zscore_CD8T")
#
# NK_temp <- GO_results$NK_baboon$LFCshrink[, c(3, 9)]
# colnames(NK_temp) <- c("Term", "Zscore_NK")
#
# LFCshrink_BP <- dplyr::full_join(CD4T_temp, CD8T_temp, by = "Term")
# LFCshrink_BP <- dplyr::full_join(LFCshrink_BP, NK_temp, by = "Term")
#
# for (x in 1:nrow(LFCshrink_BP)) {
#   if (grepl("negative regulation of", LFCshrink_BP[x, 1])) {
#     new_name <- limma::strsplit2(LFCshrink_BP[x, 1], " of ")[2]
#     new_name <- stringi::stri_join(new_name, ", negative regulation of")
#     LFCshrink_BP[x, 1] <- new_name
#   } else if (grepl("positive regulation of", LFCshrink_BP[x, 1])) {
#     new_name <- limma::strsplit2(LFCshrink_BP[x, 1], " of ")[2]
#     new_name <- stringi::stri_join(new_name, ", positive regulation of")
#     LFCshrink_BP[x, 1] <- new_name
#   } else if (grepl("regulation of", LFCshrink_BP[x, 1])) {
#     new_name <- limma::strsplit2(LFCshrink_BP[x, 1], " of ")[2]
#     new_name <- stringi::stri_join(new_name, ", regulation of")
#     LFCshrink_BP[x, 1] <- new_name
#   }
# }
#
#
# rownames(LFCshrink_BP) <- LFCshrink_BP[, 1]
# LFCshrink_BP <- LFCshrink_BP[, 2:4]
# head.DataTable(LFCshrink_BP, 10)
#
# wb <- openxlsx::createWorkbook("LFCshrink Biological Processes GO.xlsx")
#
# openxlsx::addWorksheet(wb, "sheet1")
# openxlsx::writeData(wb, "sheet1", LFCshrink_BP, rowNames = TRUE)
#
# openxlsx::saveWorkbook(wb, "LFCshrink Biological Processes GO.xlsx", overwrite = TRUE)
#
# for (n in 1:length(NK$human$LFCshrinkage$Gene_ID)) {
#   if (NK$human$LFCshrinkage$Gene_ID[n] == "ccl3") print("found ccl3")
# }
