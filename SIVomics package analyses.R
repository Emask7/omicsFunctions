roxygen2::roxygenise()


library(SIVomics)
library(openxlsx)
library(RDAVIDWebService)




sample_list <- c(
  "CD4T_T-7_15979", "CD4T_T-7_30760", "CD4T_T-7_31151",
  "CD4T_T15_15979", "CD4T_T15_30760", "CD4T_T15_31151",
  "CD8T_T-7_15979", "CD8T_T-7_30760", "CD8T_T-7_31151",
  "CD8T_T15_15979", "CD8T_T15_30760", "CD8T_T15_31151",
  "NK_T-7_15979", "NK_T-7_30760", "NK_T-7_31151",
  "NK_T15_15979", "NK_T15_30760", "NK_T15_31151"
)

humanCounts <- read.xlsx("copy - Partek_LG_RNA_Seq_20190304_raw_human_gene_counts.xlsx")
rownames(humanCounts) <- make.names(humanCounts[, 5], unique = TRUE)
humanCounts <- humanCounts[, 7:24]
colnames(humanCounts) <- c(sample_list)
head(humanCounts)

baboonCounts <- read.xlsx("raw_baboon_genecounts.xlsx")
rownames(baboonCounts) <- make.names(baboonCounts[, 5], unique = TRUE)
baboonCounts <- baboonCounts[, 9:26]
colnames(baboonCounts) <- c(sample_list)
head(baboonCounts)

humanHomologs <- read.xlsx("biomart_export_human homologs_no duplicates.xlsx")
head(humanHomologs)
colnames(humanHomologs)
humanHomologs <- humanHomologs[,c(1, 3)]
colnames(humanHomologs) <- c("Gene_ID", "Human_Gene_ID")
head(humanHomologs)



CD4T <- list(
  human = run_DESeq2(humanCounts[, 1:6], humanHomologs),
  baboon = run_DESeq2(baboonCounts[, 1:6], humanHomologs)
)

CD8T <- list(
  human = run_DESeq2(humanCounts[, 7:12], humanHomologs),
  baboon = run_DESeq2(baboonCounts[, 7:12], humanHomologs)
)

NK <- list(
  human = run_DESeq2(humanCounts[, 13:18], humanHomologs),
  baboon = run_DESeq2(baboonCounts[, 13:18], humanHomologs)
)

nrow(NK$human$Wald)
nrow(NK$human$LFCshrinkage)

nrow(NK$baboon$Wald)
nrow(NK$baboon$LFCshrinkage)



david <- DAVIDWebService(
  email = "emask@txbiomed.org",
  url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/"
)

# is.connected(david)
show(david)

GO_results <- list(
  CD4T_human = list(
    Wald = run_DAVID(david, CD4T$human$Wald, "CD4T_human_Wald"),
    LFCshrink = run_DAVID(david, CD4T$human$LFCshrinkage, "CD4T_human_LFCshrink")
  ),
  CD4T_baboon = list(
    Wald = run_DAVID(david, CD4T$baboon$Wald, "CD4T_baboon_Wald"),
    LFCshrink = run_DAVID(david, CD4T$baboon$LFCshrinkage, "CD4T_baboon_LFCshrink")
  ),
  CD8T_human = list(
    Wald = run_DAVID(david, CD8T$human$Wald, "CD8T_human_Wald"),
    LFCshrink = run_DAVID(david, CD8T$human$LFCshrinkage, "CD8T_human_LFCshrink")
  ),
  CD8T_baboon = list(
    Wald = run_DAVID(david, CD8T$baboon$Wald, "CD8T_baboon_Wald"),
    LFCshrink = run_DAVID(david, CD8T$baboon$LFCshrinkage, "CD8T_baboon_LFCshrink")
  ),
  NK_human = list(
    Wald = run_DAVID(david, NK$human$Wald, "NK_human_Wald"),
    LFCshrink = run_DAVID(david, NK$human$LFCshrinkage, "NK_human_LFCshrink")
  ),
  NK_baboon = list(
    Wald = run_DAVID(david, NK$baboon$Wald, "NK_baboon_Wald"),
    LFCshrink = run_DAVID(david, NK$baboon$LFCshrinkage, "NK_baboon_LFCshrink")
  )
)

nrow(GO_results$CD4T_human$Wald)
nrow(GO_results$CD4T_human$LFCshrink)

nrow(GO_results$CD4T_baboon$Wald)
nrow(GO_results$CD4T_baboon$LFCshrink)

nrow(GO_results$CD8T_human$Wald)
nrow(GO_results$CD8T_human$LFCshrink)

nrow(GO_results$CD8T_baboon$Wald)
nrow(GO_results$CD8T_baboon$LFCshrink)

nrow(GO_results$NK_human$Wald)
nrow(GO_results$NK_human$LFCshrink)

nrow(GO_results$NK_baboon$Wald)
nrow(GO_results$NK_baboon$LFCshrink)






CD4T_temp <- GO_results$CD4T_human$LFCshrink[, c(3, 9)]
colnames(CD4T_temp) <- c("Term", "Zscore_CD4T")

CD8T_temp <- GO_results$CD8T_baboon$LFCshrink[, c(3, 9)]
colnames(CD8T_temp) <- c("Term", "Zscore_CD8T")

NK_temp <- GO_results$NK_baboon$LFCshrink[, c(3, 9)]
colnames(NK_temp) <- c("Term", "Zscore_NK")

LFCshrink_BP <- dplyr::full_join(CD4T_temp, CD8T_temp, by = "Term")
LFCshrink_BP <- dplyr::full_join(LFCshrink_BP, NK_temp, by = "Term")

for (x in 1:nrow(LFCshrink_BP)) {
  if (grepl("negative regulation of", LFCshrink_BP[x, 1])) {
    new_name <- limma::strsplit2(LFCshrink_BP[x, 1], " of ")[2]
    new_name <- stringi::stri_join(new_name, ", negative regulation of")
    LFCshrink_BP[x, 1] <- new_name
  } else if (grepl("positive regulation of", LFCshrink_BP[x, 1])) {
    new_name <- limma::strsplit2(LFCshrink_BP[x, 1], " of ")[2]
    new_name <- stringi::stri_join(new_name, ", positive regulation of")
    LFCshrink_BP[x, 1] <- new_name
  } else if (grepl("regulation of", LFCshrink_BP[x, 1])) {
    new_name <- limma::strsplit2(LFCshrink_BP[x, 1], " of ")[2]
    new_name <- stringi::stri_join(new_name, ", regulation of")
    LFCshrink_BP[x, 1] <- new_name
  }
}


rownames(LFCshrink_BP) <- LFCshrink_BP[, 1]
LFCshrink_BP <- LFCshrink_BP[, 2:4]
head.DataTable(LFCshrink_BP, 10)

wb <- openxlsx::createWorkbook("LFCshrink Biological Processes GO.xlsx")

openxlsx::addWorksheet(wb, "sheet1")
openxlsx::writeData(wb, "sheet1", LFCshrink_BP, rowNames = TRUE)

openxlsx::saveWorkbook(wb, "LFCshrink Biological Processes GO.xlsx", overwrite = TRUE)
