library(omicsFunctions)
library(openxlsx)


# Import data -----------------------------------------------------------------
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

  humanHomologs <- read.xlsx("biomart_export_human homologs_no duplicates.xlsx")
  head(humanHomologs)
  colnames(humanHomologs)
  humanHomologs <- humanHomologs[,c(1, 3)]
  colnames(humanHomologs) <- c("Gene_ID", "Human_Gene_ID")
  head(humanHomologs)


# Differential expression analyses --------------------------------------------
  CD4T_unfiltered <- run_DESeq2(humanCounts[, 1:6], humanHomologs, FALSE)
  CD8T_unfiltered <- run_DESeq2(humanCounts[, 7:12], humanHomologs, FALSE)
  NK_unfiltered <- run_DESeq2(humanCounts[, 13:18], humanHomologs, FALSE)

  write_DEGs_to_Excel(
    CD4T_unfiltered$LFCshrinkage,
    CD8T_unfiltered$LFCshrinkage,
    NK_unfiltered$LFCshrinkage,
    "DEGs - DESeq2 LFC Shrinkage Method - human align - no LFC filter.xlsx"
  )

