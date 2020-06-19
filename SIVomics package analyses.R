library(SIVomics)
library(openxlsx)
library(RDAVIDWebService)

humanHomologs <- read.xlsx("biomart_export_human homologs_no duplicates.xlsx")
head(humanHomologs)
colnames(humanHomologs)
humanHomologs <- humanHomologs[,c(1, 3)]
colnames(humanHomologs) <- c("Gene_ID", "Human_Gene_ID")
head(humanHomologs)


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
colnames(humanCounts)

CD4T_hu_edgeR <- run_edgeR(humanCounts[, 1:6], humanHomologs, show_plots = FALSE)
CD8T_hu_edgeR <- run_edgeR(humanCounts[, 7:12], humanHomologs, show_plots = FALSE)
NK_hu_edgeR <- run_edgeR(humanCounts[, 13:18], humanHomologs, show_plots = FALSE)

CD4T_hu_DESeq2 <- run_DESeq2(humanCounts[, 1:6], humanHomologs)
CD8T_hu_DESeq2 <- run_DESeq2(humanCounts[, 7:12], humanHomologs)
NK_hu_DESeq2 <- run_DESeq2(humanCounts[, 13:18], humanHomologs)



baboonCounts <- read.xlsx("raw_baboon_genecounts.xlsx")
rownames(baboonCounts) <- make.names(baboonCounts[, 5], unique = TRUE)
baboonCounts <- baboonCounts[, 9:26]
colnames(baboonCounts) <- c(sample_list)
head(baboonCounts)
colnames(baboonCounts)

CD4T_bn_edgeR <- run_edgeR(baboonCounts[, 1:6], humanHomologs, show_plots = FALSE)
CD8T_bn_edgeR <- run_edgeR(baboonCounts[, 7:12], humanHomologs, show_plots = FALSE)
NK_bn_edgeR <- run_edgeR(baboonCounts[, 13:18], humanHomologs, show_plots = FALSE)

CD4T_bn_DESeq2 <- run_DESeq2(baboonCounts[, 1:6], humanHomologs)
CD8T_bn_DESeq2 <- run_DESeq2(baboonCounts[, 7:12], humanHomologs)
NK_bn_DESeq2 <- run_DESeq2(baboonCounts[, 13:18], humanHomologs)




CD4T <- list(
  human = list(
    
  )
)



david <- DAVIDWebService(
  email = "emask@txbiomed.org", 
  url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/"
)

is.connected(david)
show(david)

add_gene_list(david, CD4T_hu_edgeR$LRT, "CD4T_human_edgeR_LRT")
add_gene_list(david, CD4T_hu_edgeR$LRT_TREAT, "CD4T_human_edgeR_LRT_TREAT")
add_gene_list(david, CD4T_hu_DESeq2$Wald, "CD4T_human_DESeq2_Wald")
add_gene_list(david, CD4T_hu_DESeq2$LFCshrinkage, "CD4T_human_DESeq2_LFCshrink")

add_gene_list(david, CD8T_hu_edgeR$LRT, "CD8T_human_edgeR_LRT")
add_gene_list(david, CD8T_hu_edgeR$LRT_TREAT, "CD8T_human_edgeR_LRT_TREAT")
add_gene_list(david, CD8T_hu_DESeq2$Wald, "CD8T_human_DESeq2_Wald")
add_gene_list(david, CD8T_hu_DESeq2$LFCshrinkage, "CD8T_human_DESeq2_LFCshrink")

add_gene_list(david, NK_hu_edgeR$LRT, "NK_human_edgeR_LRT")
add_gene_list(david, NK_hu_edgeR$LRT_TREAT, "NK_human_edgeR_LRT_TREAT")
add_gene_list(david, NK_hu_DESeq2$Wald, "NK_human_DESeq2_Wald")
add_gene_list(david, NK_hu_DESeq2$LFCshrinkage, "NK_human_DESeq2_LFCshrink")



add_gene_list(david, CD4T_bn_edgeR$LRT, "CD4T_baboon_edgeR_LRT")
add_gene_list(david, CD4T_bn_edgeR$LRT_TREAT, "CD4T_baboon_edgeR_LRT_TREAT")
add_gene_list(david, CD4T_bn_DESeq2$Wald, "CD4T_baboon_DESeq2_Wald")
add_gene_list(david, CD4T_bn_DESeq2$LFCshrinkage, "CD4T_baboon_DESeq2_LFCshrink")

add_gene_list(david, CD8T_bn_edgeR$LRT, "CD8T_baboon_edgeR_LRT")
add_gene_list(david, CD8T_bn_edgeR$LRT_TREAT, "CD8T_baboon_edgeR_LRT_TREAT")
add_gene_list(david, CD8T_bn_DESeq2$Wald, "CD8T_baboon_DESeq2_Wald")
add_gene_list(david, CD8T_bn_DESeq2$LFCshrinkage, "CD8T_baboon_DESeq2_LFCshrink")

add_gene_list(david, NK_bn_edgeR$LRT, "NK_baboon_edgeR_LRT")
# add_gene_list(david, NK_bn_edgeR$LRT_TREAT, "NK_baboon_edgeR_LRT_TREAT")
# add_gene_list(david, NK_bn_DESeq2$Wald, "NK_baboon_DESeq2_Wald")
add_gene_list(david, NK_bn_DESeq2$LFCshrinkage, "NK_baboon_DESeq2_LFCshrink")

show(david)





CD4T_hu_edgeR_LRT_GO <- run_DAVID(david, "CD4T_human_edgeR_LRT")
CD4T_hu_edgeR_LRT_TREAT_GO <- run_DAVID(david, "CD4T_human_edgeR_LRT_TREAT")
CD4T_hu_DESeq2_Wald_GO <- run_DAVID(david, "CD4T_human_DESeq2_Wald")
CD4T_hu_DESeq2_LFCshrink_GO <- run_DAVID(david, "CD4T_human_DESeq2_LFCshrink")

CD8T_hu_edgeR_LRT_GO <- run_DAVID(david, "CD8T_human_edgeR_LRT")
CD8T_hu_edgeR_LRT_TREAT_GO <- run_DAVID(david, "CD8T_human_edgeR_LRT_TREAT")
CD8T_hu_DESeq2_Wald_GO <- run_DAVID(david, "CD8T_human_DESeq2_Wald")
CD8T_hu_DESeq2_LFCshrink_GO <- run_DAVID(david, "CD8T_human_DESeq2_LFCshrink")

NK_hu_edgeR_LRT_GO <- run_DAVID(david, "NK_human_edgeR_LRT")
NK_hu_edgeR_LRT_TREAT_GO <- run_DAVID(david, "NK_human_edgeR_LRT_TREAT")
NK_hu_DESeq2_Wald_GO <- run_DAVID(david, "NK_human_DESeq2_Wald")
NK_hu_DESeq2_LFCshrink_GO <- run_DAVID(david, "NK_human_DESeq2_LFCshrink")


CD4T_bn_edgeR_LRT_GO <- run_DAVID(david, "CD4T_baboon_edgeR_LRT")
CD4T_bn_edgeR_LRT_TREAT_GO <- run_DAVID(david, "CD4T_baboon_edgeR_LRT_TREAT")
CD4T_bn_DESeq2_Wald_GO <- run_DAVID(david, "CD4T_baboon_DESeq2_Wald")
CD4T_bn_DESeq2_LFCshrink_GO <- run_DAVID(david, "CD4T_baboon_DESeq2_LFCshrink")

CD8T_bn_edgeR_LRT_GO <- run_DAVID(david, "CD8T_baboon_edgeR_LRT")
CD8T_bn_edgeR_LRT_TREAT_GO <- run_DAVID(david, "CD8T_baboon_edgeR_LRT_TREAT")
CD8T_bn_DESeq2_Wald_GO <- run_DAVID(david, "CD8T_baboon_DESeq2_Wald")
CD8T_bn_DESeq2_LFCshrink_GO <- run_DAVID(david, "CD8T_baboon_DESeq2_LFCshrink")

NK_bn_edgeR_LRT_GO <- run_DAVID(david, "NK_baboon_edgeR_LRT")
# NK_bn_edgeR_LRT_TREAT_GO <- run_DAVID(david, "NK_baboon_edgeR_LRT_TREAT")
# NK_bn_DESeq2_Wald_GO <- run_DAVID(david, "NK_baboon_DESeq2_Wald")
NK_bn_DESeq2_LFCshrink_GO <- run_DAVID(david, "NK_baboon_DESeq2_LFCshrink")










GO_CD4T_h_Wald <- run_DAVID(david, "DESeq2_Wald_human")

head(GO_CD4T_h_Wald$BP)
head(GO_CD4T_h_Wald$HIV_interaction)
head(GO_CD4T_h_Wald$KEGG_pathway)
