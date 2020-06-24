library(openxlsx)
library(SIVomics)


immune_process <- list(
  CD4T = get_term_DEGs(
    c("immune system process"),
    GO_LFCshrink_human$CD4T, CD4T$human_DESeq2$LFCshrinkage,
    "Immune system process DEGs.xlsx", "CD4T", FALSE
  ),
  CD8T = get_term_DEGs(
    c("immune system process"),
    GO_LFCshrink_human$CD8T, CD8T$human_DESeq2$LFCshrinkage,
    "Immune system process DEGs.xlsx", "CD8T", FALSE
  ),
  NK = get_term_DEGs(
    c("immune system process"),
    GO_LFCshrink_human$NK, NK$human_DESeq2$LFCshrinkage,
    "Immune system process DEGs.xlsx", "NK", FALSE
  ),
  regulation_of = list(
    CD4T = get_term_DEGs(
      c("regulation of immune system process"),
      GO_LFCshrink_human$CD4T, CD4T$human_DESeq2$LFCshrinkage,
      "Immune System Processes DEGs.xlsx", "CD4T - regulation of ISP", FALSE
    ),
    CD8T = get_term_DEGs(
      c("regulation of immune system process"),
      GO_LFCshrink_human$CD8T, CD8T$human_DESeq2$LFCshrinkage,
      "Immune System Processes DEGs.xlsx", "CD8T - regulation of ISP", FALSE
    ),
    NK = get_term_DEGs(
      c("regulation of immune system process"),
      GO_LFCshrink_human$NK, NK$human_DESeq2$LFCshrinkage,
      "Immune System Processes DEGs.xlsx", "NK - regulation of ISP", FALSE
    )
  ),
  positive_regulation_of = list(
    CD4T = get_term_DEGs(
      c("positive regulation of immune system process"),
      GO_LFCshrink_human$CD4T, CD4T$human_DESeq2$LFCshrinkage,
      "Immune System Processes DEGs.xlsx",
      "CD4T - positive reg.", FALSE
    ),
    CD8T = get_term_DEGs(
      c("positive regulation of immune system process"),
      GO_LFCshrink_human$CD8T, CD8T$human_DESeq2$LFCshrinkage,
      "Immune System Processes DEGs.xlsx",
      "CD8T - positive reg.", FALSE
    ),
    NK = get_term_DEGs(
      c("positive regulation of immune system process"),
      GO_LFCshrink_human$NK, NK$human_DESeq2$LFCshrinkage,
      "Immune System Processes DEGs.xlsx",
      "NK - positive reg.", FALSE
    )
  )
)



inflammatory_response <- list(
  CD4T = get_term_DEGs(
    c("inflammatory response"),
    GO_LFCshrink_human$CD4T, CD4T$human_DESeq2$LFCshrinkage,
    "Inflammatory Response DEGs.xlsx", "CD4T", FALSE
  ),
  CD8T = get_term_DEGs(
    c("inflammatory response"),
    GO_LFCshrink_human$CD8T, CD8T$human_DESeq2$LFCshrinkage,
    "Inflammatory Response DEGs.xlsx", "CD8T", FALSE
  ),
  NK = get_term_DEGs(
    c("inflammatory response"),
    GO_LFCshrink_human$NK, NK$human_DESeq2$LFCshrinkage,
    "Inflammatory Response DEGs.xlsx", "NK", FALSE
  )
)
