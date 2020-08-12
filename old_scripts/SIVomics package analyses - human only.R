roxygen2::roxygenise()


library(omicsFunctions)
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
rownames(tCounts) <- make.names(tCounts[, 5], unique = TRUE)
tCounts <- tCounts[, c(6, 9, 16, 20:37)]
colnames(tCounts) <- c("Gene_ID", "Ensembl_Gene_ID", "Ensembl_Transcript_ID", sample_list)
head(tCounts)


colnames(tCounts)


scater::librarySizeFactors(tCounts[, 4:9])
scater::librarySizeFactors(tCounts[, 10:15])
scater::librarySizeFactors(tCounts[, 16:21])

# Differential expression analyses --------------------------------------------
DE_genes <- list(
  CD4T = run_DESeq2_humanOnly(gCounts[, 1:6], "genes"),
  CD8T = run_DESeq2_humanOnly(gCounts[, 7:12], "genes"),
  NK = run_DESeq2_humanOnly(gCounts[, 13:18], "genes")
)

DE_transcripts <- list(
  CD4T = run_DESeq2_humanOnly(tCounts[, 4:9], "transcripts"),
  CD8T = run_DESeq2_humanOnly(tCounts[, 10:15], "transcripts"),
  NK = run_DESeq2_humanOnly(tCounts[, 16:21], "transcripts")
)


write_DEGs_to_Excel(
  DE_genes$CD4T$filtered_results,
  DE_genes$CD8T$filtered_results,
  DE_genes$NK$filtered_results,
  "DEGs - DESeq2 LFC Shrinkage Method - human only.xlsx"
)


transcript_info <- data.frame(rownames(tCounts), tCounts[, 1:3])
colnames(transcript_info) <- c("Transcript_ID", "Transcript", "Gene", "Ensembl_Gene_ID")
head(transcript_info)

write_DETs_to_Excel(
  DE_transcripts$CD4T$filtered_results,
  DE_transcripts$CD8T$filtered_results,
  DE_transcripts$NK$filtered_results,
  transcript_info,
  "DE Transcripts - DESeq2 LFC Shrinkage Method - human only.xlsx"
)



# Summarize differential expression results -----------------------------------
get_DEG_number <- function(x) {
  up <- nrow(subset(x, x$LFC >= 1 & x$padj <= 0.05))
  down <- nrow(subset(x, x$LFC <= -1 & x$padj <= 0.05))
  c(up, down)
}

CD4T_summary = data.frame(
  get_DEG_number(CD4T$human_edgeR), get_DEG_number(CD4T$baboon_edgeR),
  get_DEG_number(CD4T$human_DESeq2$Wald),
  get_DEG_number(CD4T$baboon_DESeq2$Wald),
  get_DEG_number(CD4T$human_DESeq2$LFCshrinkage),
  get_DEG_number(CD4T$baboon_DESeq2$LFCshrinkage),
  row.names = c("Up (LFC >= 1)", "Down (LFC <= -1)")
)
colnames(CD4T_summary) <- c(
  "edgeR_human", "edgeR_baboon", "Wald_human", "Wald_baboon",
  "LFCshrinkage_human", "LFCshrinkage_baboon"
)

CD8T_summary = data.frame(
  get_DEG_number(CD8T$human_edgeR), get_DEG_number(CD8T$baboon_edgeR),
  get_DEG_number(CD8T$human_DESeq2$Wald),
  get_DEG_number(CD8T$baboon_DESeq2$Wald),
  get_DEG_number(CD8T$human_DESeq2$LFCshrinkage),
  get_DEG_number(CD8T$baboon_DESeq2$LFCshrinkage),
  row.names = c("Up (LFC >= 1)", "Down (LFC <= -1)")
)
colnames(CD8T_summary) <- c(
  "edgeR_human", "edgeR_baboon", "Wald_human", "Wald_baboon",
  "LFCshrinkage_human", "LFCshrinkage_baboon"
)

NK_summary = data.frame(
  get_DEG_number(NK$human_edgeR), get_DEG_number(NK$baboon_edgeR),
  get_DEG_number(NK$human_DESeq2$Wald),
  get_DEG_number(NK$baboon_DESeq2$Wald),
  get_DEG_number(NK$human_DESeq2$LFCshrinkage),
  get_DEG_number(NK$baboon_DESeq2$LFCshrinkage),
  row.names = c("Up (LFC >= 1)", "Down (LFC <= -1)")
)
colnames(NK_summary) <- c(
  "edgeR_human", "edgeR_baboon", "Wald_human", "Wald_baboon",
  "LFCshrinkage_human", "LFCshrinkage_baboon"
)

wb1 <- createWorkbook("DEG Summary.xlsx")
addWorksheet(wb1, "sheet1")
writeData(wb1, "sheet1", "CD4T", startCol = 1, startRow = 2)
writeData(wb1, "sheet1", "CD8T", startCol = 1, startRow = 5)
writeData(wb1, "sheet1", "NK", startCol = 1, startRow = 8)
writeData(
  wb1, "sheet1", CD4T_summary,
  startCol = 2, startRow = 1, colNames = TRUE, rowNames = TRUE
)
writeData(
  wb1, "sheet1", CD8T_summary,
  startCol = 2, startRow = 5, colNames = FALSE, rowNames = TRUE
)
writeData(
  wb1, "sheet1", NK_summary,
  startCol = 2, startRow = 8, colNames = FALSE, rowNames = TRUE
)
saveWorkbook(wb1, "DEG Summary.xlsx", overwrite = TRUE)
rm(wb1)

# DAVID Gene Ontology analyses ------------------------------------------------
david <- DAVIDWebService(
  email = "emask@txbiomed.org",
  url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/"
)

is.connected(david)
show(david)


GO_LFCshrink_human <- list(
  CD4T = run_DAVID(david, CD4T$human_DESeq2$LFCshrinkage, "CD4T_human_LFCshrink"),
  CD8T = run_DAVID(david, CD8T$human_DESeq2$LFCshrinkage, "CD8T_human_LFCshrink"),
  NK = run_DAVID(david, NK$human_DESeq2$LFCshrinkage, "NK_human_LFCshrink")
)

CD4T_GO <- run_DAVID(david, CD4T$human_DESeq2$LFCshrinkage, "CD4T_human_LFCshrink")
CD8T_GO <- run_DAVID(david, CD8T$human_DESeq2$LFCshrinkage, "CD8T_human_LFCshrink")
NK_GO <- run_DAVID(david, NK$human_DESeq2$LFCshrinkage, "NK_human_LFCshrink")

NK_GO_backup <- NK_GO

# Troubleshooting the "Read timed out" error --------------------------------
setTimeOut(david, 50000)
getTimeOut(david)

setHttpProtocolVersion(david, "HTTP/1.0")
getHttpProtocolVersion(david)





# GO_results <- list(
#   CD4T_human = list(
#     Wald = run_DAVID(david, CD4T$human$Wald, "CD4T_human_Wald"),
#     LFCshrink = run_DAVID(david, CD4T$human$LFCshrinkage, "CD4T_human_LFCshrink")
#   ),
#   CD4T_baboon = list(
#     Wald = run_DAVID(david, CD4T$baboon$Wald, "CD4T_baboon_Wald"),
#     LFCshrink = run_DAVID(david, CD4T$baboon$LFCshrinkage, "CD4T_baboon_LFCshrink")
#   ),
#   CD8T_human = list(
#     Wald = run_DAVID(david, CD8T$human$Wald, "CD8T_human_Wald"),
#     LFCshrink = run_DAVID(david, CD8T$human$LFCshrinkage, "CD8T_human_LFCshrink")
#   ),
#   CD8T_baboon = list(
#     Wald = run_DAVID(david, CD8T$baboon$Wald, "CD8T_baboon_Wald"),
#     LFCshrink = run_DAVID(david, CD8T$baboon$LFCshrinkage, "CD8T_baboon_LFCshrink")
#   ),
#   NK_human = list(
#     Wald = run_DAVID(david, NK$human$Wald, "NK_human_Wald"),
#     LFCshrink = run_DAVID(david, NK$human$LFCshrinkage, "NK_human_LFCshrink")
#   ),
#   NK_baboon = list(
#     Wald = run_DAVID(david, NK$baboon$Wald, "NK_baboon_Wald"),
#     LFCshrink = run_DAVID(david, NK$baboon$LFCshrinkage, "NK_baboon_LFCshrink")
#   )
# )







# nrow(GO_results$CD4T_human$Wald)
# nrow(GO_results$CD4T_human$LFCshrink)
#
# nrow(GO_results$CD4T_baboon$Wald)
# nrow(GO_results$CD4T_baboon$LFCshrink)
#
# nrow(GO_results$CD8T_human$Wald)
# nrow(GO_results$CD8T_human$LFCshrink)
#
# nrow(GO_results$CD8T_baboon$Wald)
# nrow(GO_results$CD8T_baboon$LFCshrink)
#
# nrow(GO_results$NK_human$Wald)
# nrow(GO_results$NK_human$LFCshrink)
#
# nrow(GO_results$NK_baboon$Wald)
# nrow(GO_results$NK_baboon$LFCshrink)
#
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
r
