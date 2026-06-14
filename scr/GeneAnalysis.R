# Raw CEL files
# ↓
# Quality Control
# ↓
# RMA normalization
# ↓
# PCA & clustering
# ↓
# limma differential expression
# ↓
# Volcano plot
# ↓
# ComplexHeatmap of top DEGs
# ↓
# fgsea Hallmark pathways
# ↓
# clusterProfiler GO analysis
# ↓
# ReactomePA pathway analysis
# ↓
# Publication figures with ggplot2 + patchwork


# 1. DOWNLOAD DATA ------------------

library(GEOquery)

getGEOSuppFiles("GSE106571")


# 1. BACKGROUND CORRECTION AND NORMALIZATION ---------------
library(oligo)
eset <- rma(rawData)



