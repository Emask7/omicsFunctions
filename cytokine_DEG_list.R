library(openxlsx)
library(SIVomics)


cytokine_production <- list(
  CD4T = get_term_DEGs(
    c("cytokine production"),
    CD4T_GO, CD4T$human_DESeq2$LFCshrinkage,
    "Cytokine Production DEGs.xlsx", "CD4T", FALSE
  ),
  CD8T = get_term_DEGs(
    c("cytokine production"),
    CD8T_GO, CD8T$human_DESeq2$LFCshrinkage,
    "Cytokine Production DEGs.xlsx", "CD8T", FALSE
  ),
  NK = get_term_DEGs(
    c("cytokine production"),
    GO_LFCshrink_human$NK, NK$human_DESeq2$LFCshrinkage,
    "Cytokine Production DEGs.xlsx", "NK", FALSE
  )
)

regulation_of_cytokine_production <- list(
  CD4T = get_term_DEGs(
    c("regulation of cytokine production"),
    CD4T_GO, CD4T$human_DESeq2$LFCshrinkage,
    "Regulation of Cytokine Production DEGs.xlsx", "CD4T", FALSE
  ),
  CD8T = get_term_DEGs(
    c("regulation of cytokine production"),
    CD8T_GO, CD8T$human_DESeq2$LFCshrinkage,
    "Regulation of Cytokine Production DEGs.xlsx", "CD8T", FALSE
  ),
  NK = get_term_DEGs(
    c("regulation of cytokine production"),
    GO_LFCshrink_human$NK, NK$human_DESeq2$LFCshrinkage,
    "Regulation of Cytokine Production DEGs.xlsx", "NK", FALSE
  ),
  positive = list(
    CD4T = get_term_DEGs(
      c("positive regulation of cytokine production"),
      CD4T_GO, CD4T$human_DESeq2$LFCshrinkage,
      "Regulation of Cytokine Production DEGs.xlsx",
      "CD4T - positive reg.", FALSE
    ),
    CD8T = get_term_DEGs(
      c("positive regulation of cytokine production"),
      CD8T_GO, CD8T$human_DESeq2$LFCshrinkage,
      "Regulation of Cytokine Production DEGs.xlsx",
      "CD8T - positive reg.", FALSE
    ),
    NK = get_term_DEGs(
      c("positive regulation of cytokine production"),
      GO_LFCshrink_human$NK, NK$human_DESeq2$LFCshrinkage,
      "Regulation of Cytokine Production DEGs.xlsx",
      "NK - positive reg.", FALSE
    )
  ),
  negative = list(
    CD4T = get_term_DEGs(
      c("negative regulation of cytokine production"),
      CD4T_GO, CD4T$human_DESeq2$LFCshrinkage,
      "Regulation of Cytokine Production DEGs.xlsx",
      "CD4T - negative reg.", FALSE
    ),
    CD8T = get_term_DEGs(
      c("negative regulation of cytokine production"),
      CD8T_GO, CD8T$human_DESeq2$LFCshrinkage,
      "Regulation of Cytokine Production DEGs.xlsx",
      "CD8T - negative reg.", FALSE
    ),
    NK = get_term_DEGs(
      c("negative regulation of cytokine production"),
      GO_LFCshrink_human$NK, NK$human_DESeq2$LFCshrinkage,
      "Regulation of Cytokine Production DEGs.xlsx",
      "NK - negative reg.", FALSE
    )
  )
)





# cytokine_production <- list(
#   cytokine_production = get_term_DEGs(
#     c("cytokine production"),
#     GO_results$CD4T_human$LFCshrink$GOplot_results,
#     CD4T$human$LFCshrinkage,
#     "equal_to"
#   ),
#   regulation_of = get_term_DEGs(
#     c("regulation of cytokine production"),
#     GO_results$CD4T_human$LFCshrink$GOplot_results,
#     CD4T$human$LFCshrinkage,
#     "equal_to"
#   ),
#   neg_regulation_of = get_term_DEGs(
#     c("negative regulation of cytokine production"),
#     GO_results$CD4T_human$LFCshrink$GOplot_results,
#     CD4T$human$LFCshrinkage,
#     "equal_to"
#   ),
#   pos_regulation_of = get_term_DEGs(
#     c("positive regulation of cytokine production"),
#     GO_results$CD4T_human$LFCshrink$GOplot_results,
#     CD4T$human$LFCshrinkage,
#     "equal_to"
#   )
# )
#
# response_to_cytokine_CD4T <- get_term_DEGs(
#   c("response to cytokine"),
#   GO_results$CD4T_human$LFCshrink$GOplot_results,
#   CD4T$human$LFCshrinkage,
#   "equal_to"
# )
# colnames(response_to_cytokine_CD4T) <- c("Term", "Gene_ID", "LFC (CD4T)", "padj (CD4T)")
#
# response_to_cytokine_CD8T <- get_term_DEGs(
#   c("response to cytokine"),
#   GO_results$CD8T_human$LFCshrink$GOplot_results,
#   CD8T$human$LFCshrinkage,
#   "equal_to"
# )
# colnames(response_to_cytokine_CD8T) <- c("Term", "Gene_ID", "LFC (CD8T)", "padj (CD8T)")
#
#
# response_to_cytokine <- dplyr::full_join(
#   response_to_cytokine_CD4T, response_to_cytokine_CD8T
# )
#
#
#
#
#
# immune_system_CD4T <- get_term_DEGs(
#   c("immune", "inflamm"),
#   GO_results$CD4T_human$LFCshrink$GOplot_results,
#   CD4T$human$LFCshrinkage,
#   "contains"
# )
# colnames(immune_system_CD4T) <- c("Term", "Gene_ID", "LFC (CD4T)", "padj (CD4T)")
#
# immune_system_CD8T <- get_term_DEGs(
#   c("immune", "inflamm"),
#   GO_results$CD8T_human$LFCshrink$GOplot_results,
#   CD8T$human$LFCshrinkage,
#   "contains"
# )
# colnames(immune_system_CD8T) <- c("Term", "Gene_ID", "LFC (CD8T)", "padj (CD8T)")
#
# immune_system <- dplyr::full_join(
#   immune_system_CD4T, immune_system_CD8T
# )
#
#
#
#
#
#
#
# wb <- createWorkbook("GO Term DEG Lists.xlsx")
#
# addWorksheet(wb, "Cytokine Production")
# writeData(wb, "Cytokine Production", cytokine_production$cytokine_production)
#
# addWorksheet(wb, "Reg. of Cytokine Production")
# writeData(wb, "Reg. of Cytokine Production", cytokine_production$regulation_of)
#
# addWorksheet(wb, "- Reg. of Cytokine Production")
# writeData(
#   wb, "- Reg. of Cytokine Production", cytokine_production$neg_regulation_of
# )
#
# addWorksheet(wb, "+ Reg. of Cytokine Production")
# writeData(
#   wb, "+ Reg. of Cytokine Production", cytokine_production$pos_regulation_of
# )
#
# addWorksheet(wb, "Response to Cytokine")
# writeData(wb, "Response to Cytokine", response_to_cytokine)
#
# addWorksheet(wb, "Immune Response")
# writeData(wb, "Immune Response", immune_system)
#
# saveWorkbook(wb, "GO Term DEG Lists.xlsx", overwrite = TRUE)
