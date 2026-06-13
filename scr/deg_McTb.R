
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142756

# INSTALL PACKAGES ======
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "GEOquery",
  "limma",
  "oligo",
  "affy",
  "pheatmap",
  "EnhancedVolcano",
  "Biobase",
  "STRINGdb",
  "enrichplot"
))


# DOWNLOAD DATASET ======
library(GEOquery)

gse <- getGEO("GSE142756", GSEMatrix = TRUE)

expr <- exprs(gse[[1]])
pheno <- pData(gse[[1]])

dim(expr)
head(pheno)
boxplot(expr[,1:10])
boxplot(log2(expr[,1:30]))

# CHANGE ID OF GENES =====================
# Download the platform
library(GEOquery)
library(stringr)

gpl <- getGEO("GPL27959")
tab <- Table(gpl)
colnames(tab)
head(tab)
View(tab)

# extract gene IDs
tab$locus_tag <- str_extract(
  tab$SPOT_Id,
  "(?<=locus_tag=)[^ ]+"
)

# exxtract locus Tags
tab$GeneID <- str_extract(
  tab$SPOT_Id,
  "(?<=GeneID:)\\d+"
)

# extract gene symbol
tab$gene_symbol <- str_extract(
  tab$SPOT_Id,
  "(?<=Name=)[^ ]+"
)

# extract gene biotype
tab$gene_biotype <- str_extract(
  tab$SPOT_Id,
  "(?<=biotype=)[^ ]+"
)

# check
head(tab[, c("ID","Indentifier Source","gene_symbol","locus_tag","gene_biotype")])

# Conting table of genes by biotype
table(tab$gene_biotype)

# There are NA in the data
summary(tab)

# Check NA
tmp <- tab[is.na(tab$locus_tag),]
head(tmp[,c("ID","gene_symbol","locus_tag","gene_biotype")])
head(tab[, c("ID","gene_symbol","locus_tag","gene_biotype")])

# remove no-protein coding
table(tab$gene_biotype)
tab2 <- tab[tab$gene_biotype=="protein_coding",]
dim(tab2)
table(tab2$gene_biotype)
summary(tab2$locus_tag)
tab2 <- tab2[!is.na(tab2$locus_tag),]

# preserve only genes in my expression matrix
idx <- match(tab2$ID,rownames(expr))
expr2 <- expr[idx,]

# rename expression matrix
idx <- match(rownames(expr), tab$ID)  # index of gene ID with rowname of expression matrix
rownames(expr) <- tab$locus_tag[idx]  # change the names

head(expr)

# SELECT SAMPLE TO ANALYSE =====
sample_names <- colnames(expr)
sample_names

no.trt <- c(
  "GSM4239580",
  "GSM4239588",
  "GSM4239596",
  "GSM4239604"
)

sptm <- c(
  "GSM4239617",
  "GSM4239618",
  "GSM4239619",
  "GSM4239620",
  "GSM4239621",
  "GSM4239622",
  "GSM4239623",
  "GSM4239624",
  "GSM4239625"
)


selected <- c(no.trt, sptm)

expr_sub <- expr[, selected]


# CREATE DESIGN MATRIX ====
group <- factor(c(
  rep("no.trt", 4),
  rep("sptm", 9)
))

design <- model.matrix(~0 + group)

colnames(design) <- levels(group)

design

# DIFFERENTIAL EXPRESSION USING LIMMA  ====
library(limma)

fit <- lmFit(expr_sub, design)

contrast.matrix <- makeContrasts(
  sptm - no.trt,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

results <- topTable(
  fit2,
  number = Inf,
  adjust.method = "BH"
)

head(results)

# SAVE RESULTS =====
write.csv(
  results,
  "degs_MTb.csv"
)

# SIGNIFICANCE GENES =====
deg <- subset(
  results,
  adj.P.Val < 0.05 &
    abs(logFC) > 1  
)

nrow(deg)
head(deg)

# VOLCANO PLOT =====
library(EnhancedVolcano)

EnhancedVolcano(
  results,
  lab = results$ID,
  x = "logFC",
  y = "adj.P.Val",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "BEFORE vs AFTER Treatment"
)

# HEAT MAP TOP DEGs =======
library(pheatmap)

topgenes <- deg$ID[!is.na(deg$ID)]

mat <- expr_sub[topgenes, ]

mat <- t(scale(t(mat)))

annotation <- data.frame(
  Group = group
)

rownames(annotation) <- colnames(mat)

pheatmap(
  mat,
  annotation_col = annotation,
  show_rownames = TRUE
)

# PCA ANALYSIS ============
pca <- prcomp(t(expr_sub))

plot(
  pca$x[,1],
  pca$x[,2],
  col = c(rep("blue",3),rep("red",3)),
  pch = 19,
  xlab = "PC1",
  ylab = "PC2"
)

legend(
  "topright",
  legend = levels(group),
  col = c("blue","red"),
  pch = 19
)

# UMAP =====================
library(uwot)
library(ggplot2)

expr_t <- t(expr_sub)

# Scale data
expr_scaled <- scale(expr_t)

umap_coords <- umap(
  expr_t,
  n_neighbors = 3,
  min_dist = 0.05,
  n_components = 2,
  metric = "correlation"
)

umap_df <- data.frame(
  UMAP1 = umap_coords[,1],
  UMAP2 = umap_coords[,2]
)

# color
umap_df$Group <- ifelse(
  rownames(umap_df) %in% no.trt,
  "No Treatment",
  ifelse(
    rownames(umap_df) %in% sptm,
    "SPTM",
    "Other"
  )
)


ggplot(umap_df, aes(x= UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 4) +
  theme_linedraw()


# GSEA ====================================
library(KEGGREST)

mtb_genes <- keggLink("pathway", "mtu")
head(mtb_genes)

# term to gene
term2gene <- data.frame(
  pathway = sub("path:", "", mtb_genes),
  gene = sub("mtu:", "", names(mtb_genes))
)
head(term2gene)

# Run enrichment analysis
library(clusterProfiler)
genes <- deg$ID[!is.na(deg$ID)]
kk <- enricher(
  gene = genes,
  TERM2GENE = term2gene,
  pvalueCutoff = 0.05
)

# plot
dotplot(kk, showCategory = 20)
barplot(kk, showCategory = 15)
cnetplot(kk)

# GENE SET ENRICHMENT ANALYSIS ===========
geneList <- results$logFC
names(geneList) <- results$ID
geneList <- sort(geneList, decreasing = TRUE)

gsea <- GSEA(
  geneList,
  TERM2GENE = term2gene
)
# plot
ridgeplot(gsea)
dotplot(gsea)
library(enrichplot)
gseaplot2(gsea, 1)

# MYCOBROWSER CATEGORIES ============
myco_t2g <- data.frame(
  Category = annotation$Functional_Category,
  Gene = annotation$Locus
)

enricher(
  gene = genes,
  TERM2GENE = myco_t2g
)


# FOR PROTEIN INTERACTION ==============
library(STRINGdb)

string_db <- STRINGdb$new(
  version = "12",
  species = 1773
)

