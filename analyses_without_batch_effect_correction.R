CD4T_noBatch <- run_DESeq2_noBatchCorrection(humanCounts[, 1:6], humanHomologs, TRUE)
CD8T_noBatch <- run_DESeq2_noBatchCorrection(humanCounts[, 7:12], humanHomologs, TRUE)
NK_noBatch <- run_DESeq2_noBatchCorrection(humanCounts[, 13:18], humanHomologs, TRUE)

CD4T_noBatch$dds
