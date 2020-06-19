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
  human = list(
    edgeR = run_edgeR(humanCounts[, 1:6], humanHomologs, show_plots = FALSE),
    DESeq2 = run_DESeq2(humanCounts[, 1:6], humanHomologs)
  ),
  baboon = list(
    edgeR = run_edgeR(baboonCounts[, 1:6], humanHomologs, show_plots = FALSE),
    DESeq2 = run_DESeq2(baboonCounts[, 1:6], humanHomologs)
  )
)

CD8T <- list(
  human = list(
    edgeR = run_edgeR(humanCounts[, 7:12], humanHomologs, show_plots = FALSE),
    DESeq2 = run_DESeq2(humanCounts[, 7:12], humanHomologs)
  ),
  baboon = list(
    edgeR = run_edgeR(baboonCounts[, 7:12], humanHomologs, show_plots = FALSE),
    DESeq2 = run_DESeq2(baboonCounts[, 7:12], humanHomologs)
  )
)

NK <- list(
  human = list(
    edgeR = run_edgeR(humanCounts[, 13:18], humanHomologs, show_plots = FALSE),
    DESeq2 = run_DESeq2(humanCounts[, 13:18], humanHomologs)
  ),
  baboon = list(
    edgeR = run_edgeR(baboonCounts[, 13:18], humanHomologs, show_plots = FALSE),
    DESeq2 = run_DESeq2(baboonCounts[, 13:18], humanHomologs)
  )
)

nrow(NK$human$edgeR$results)
nrow(NK$human$DESeq2$Wald)
nrow(NK$human$DESeq2$LFCshrinkage)

nrow(NK$baboon$edgeR$results)
nrow(NK$baboon$DESeq2$Wald)
nrow(NK$baboon$DESeq2$LFCshrinkage)



david <- DAVIDWebService(
  email = "emask@txbiomed.org",
  url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/"
)

is.connected(david)
show(david)

GO_results <- list(
  CD4T_human = list(
    edgeR = run_DAVID(david, CD4T$human$edgeR$results, "CD4T_human_edgeR"),
    DESeq2_Wald = run_DAVID(david, CD4T$human$DESeq2$Wald, "CD4T_human_DESeq2_Wald"),
    DESeq2_LFCshrink = run_DAVID(david, CD4T$human$DESeq2$LFCshrinkage, "CD4T_human_DESeq2_LFCshrink")
  ),
  CD4T_baboon = list(
    edgeR = run_DAVID(david, CD4T$baboon$edgeR$results, "CD4T_baboon_edgeR"),
    DESeq2_Wald = run_DAVID(david, CD4T$baboon$DESeq2$Wald, "CD4T_baboon_DESeq2_Wald"),
    DESeq2_LFCshrink = run_DAVID(david, CD4T$baboon$DESeq2$LFCshrinkage, "CD4T_baboon_DESeq2_LFCshrink")
  ),
  CD8T_human = list(
    edgeR = run_DAVID(david, CD8T$human$edgeR$results, "CD8T_human_edgeR"),
    DESeq2_Wald = run_DAVID(david, CD8T$human$DESeq2$Wald, "CD8T_human_DESeq2_Wald"),
    DESeq2_LFCshrink = run_DAVID(david, CD8T$human$DESeq2$LFCshrinkage, "CD8T_human_DESeq2_LFCshrink")
  ),
  CD8T_baboon = list(
    edgeR = run_DAVID(david, CD8T$baboon$edgeR$results, "CD8T_baboon_edgeR"),
    DESeq2_Wald = run_DAVID(david, CD8T$baboon$DESeq2$Wald, "CD8T_baboon_DESeq2_Wald"),
    DESeq2_LFCshrink = run_DAVID(david, CD8T$baboon$DESeq2$LFCshrinkage, "CD8T_baboon_DESeq2_LFCshrink")
  ),
  NK_human = list(
    edgeR = run_DAVID(david, NK$human$edgeR$results, "NK_human_edgeR"),
    DESeq2_Wald = run_DAVID(david, NK$human$DESeq2$Wald, "NK_human_DESeq2_Wald"),
    DESeq2_LFCshrink = run_DAVID(david, NK$human$DESeq2$LFCshrinkage, "NK_human_DESeq2_LFCshrink")
  ),
  NK_baboon = list(
    # edgeR = run_DAVID(david, NK$baboon$edgeR$results, "NK_baboon_edgeR"),
    # DESeq2_Wald = run_DAVID(david, NK$baboon$DESeq2$Wald, "NK_baboon_DESeq2_Wald"),
    DESeq2_LFCshrink = run_DAVID(david, NK$baboon$DESeq2$LFCshrinkage, "NK_baboon_DESeq2_LFCshrink")
  )
)







