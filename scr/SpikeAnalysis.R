# INSTALL PACKAGES =====================
# from bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c(
  "rentrez",
  "Biostrings",
  "DECIPHER",
  "phangorn",
  "ape",
  "msa"
))

# from CRAN
install.packages(c(
  "ggplot2",
  "dplyr",
  "viridis",
  "progress",
  "tidyverse",
  "hrbrthemes"
))

# LOAD LIBRARIES =======================
library(rentrez)
library(Biostrings)
library(DECIPHER)
library(phangorn)
library(ape)
library(ggplot2)
library(dplyr)
library(viridis)

# DOWNLOAD SPIKE sequence ============
# Search sequences
search <- entrez_search(
  db = "protein",
  term = 'Severe acute respiratory syndrome coronavirus 2[Organism] AND Spike[Name]',
  retmax = 2000,
  use_history = TRUE
)
length(search$ids)

# Download spike sequences
fasta <- entrez_fetch(
  db = "protein",
  id = search$ids,
  rettype = "fasta",
  retmode = "text",
)
# NOTE: Get eroor because there are many sequences, so 
# we will split the download sequences in batchs

# Download sequences in batches
ids <- search$ids         # define Ids vector
web_history_obj <- search$web_history # define "history"

# Create progress bar
library(progress)
pb <- progress_bar$new(
  format = "Batch :current/:total [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
  total = length(seq(0, 2000, by = 50)),
  clear = FALSE,
  width = 90
)

# Create address of output directory
outfile <- file.path(getwd(), 
                     "data", 
                     "Spike2000.fasta")


for(start in seq(0, 2000, by = 50)) {
  # Download fasta files (batch)
  fasta <- entrez_fetch(
    db = "protein",
    web_history = search$web_history,
    rettype = "fasta",
    retmode = "text",
    retstart = start,
    retmax = 50
  )
  
  # Save it
  cat(
    fasta,
    file = outfile,
    append = TRUE
  )
  
  # update progress bar
  pb$tick()
  
  # to avoid exced the NCBI's rate-request limit
  Sys.sleep(1)
}


# EDIT SEQUENCES =====================
# upload sequences like readFasta format
seqs <- readAAStringSet(
  outfile
)
length(seqs)
typeof(seqs)
class(seqs)
summary(seqs)

# control quality
ambiguous <- sapply(seqs,function(x){sum(strsplit(as.character(x),"")[[1]]%in%c("X","B","Z","J","*"))})
seqs <- seqs[ambiguous < 1]
length(seqs)

# remove duplicates
seqs <- seqs[!duplicated(as.character(seqs))]
length(seqs)

# remove incomplete sequences
len <- data.frame(x=width(seqs))

# Libraries
library(tidyverse)
library(hrbrthemes)
library(ggplot2)
# plot
ggplot( len, aes(x=x)) +
  geom_histogram( binwidth=100, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("SPIKE") +
  labs(subtitle = "(Sequence legnths)") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  )
# Remove outliers
seqs <- seqs[
  len$x > 1250 &
    len$x < 1290
]
length(seqs)

p <- boxplot(len$x)
p$out
seqs <- seqs[-p$out]
length(seqs)

len <- data.frame(x=width(seqs))
summary(len$x)
ggplot( len, aes(x=x)) +
  geom_histogram( binwidth=2, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("SPIKE") +
  labs(subtitle = "(Sequence legnths)") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  )


# SEQUENCE ALIGNMENT ===================
library(msa)
alignment <- msa(
  seqs,
  method = "ClustalOmega", # Omega more efficient for >1000 seq
  verbose = TRUE
)

library(Biostrings)
aligned <- msaConvert(
  alignment,
  type = "seqinr::alignment"
)

# Evaluate alignment quality
library(DECIPHER)
alignment2 <- AAStringSet(alignment)
BrowseSeqs(alignment2)

# save
writeXStringSet(
  alignment2,
  filepath = file.path(getwd(),"result","Spike_alignment.fasta")
)


# Calculated gaps
mat <- as.matrix(alignment2)
gap_fraction <- apply(mat,1,function(x) mean(x=="-"))
hist(gap_fraction)
summary(gap_fraction)

# remove sequence with many gaps
keep <- gap_fraction < 0.0731
alignment2 <- alignment2[keep,]
BrowseSeqs(alignment2)


# CONCESUS SEQUENCE =============
consensus <- ConsensusSequence(
  alignment2
)

# PHYLOGENY ==============================

# transform data
mat <- as.matrix(alignment2)

phy <- phyDat(
               mat,
                type = "AA"
                )

# calculate distance matrix
dm <- dist.ml(phy)

# Neightbor join tree
tree <- NJ(dm)
tree <- midpoint(tree) # root the midpoint

# Optimize tree
fit <- pml(tree,phy)
fit <- optim.pml(
  fit,
  model="JTT",
  optGamma=TRUE,
  optInv=TRUE
)

# best tree
tree.ml <- fit$tree

# plot
library(ggtree)
ggtree(tree)
ggtree(tree.ml)


# Bootstrap tree (take hours)
# bs <- bootstrap.pml(fit,
#                     bs=100,
#                     optNni=TRUE
#                     )


write.tree(            # save the tree
  tree,
  file = file.path(getwd(),"result","NJ_Spike_tree.nwk")
)

write.tree(            # save the tree
  tree.ml,
  file = file.path(getwd(),"result","LH_Spike_tree.nwk")
)

# plot phylogeny tree
pdf(
  file.path(getwd(),"plot","Spike_phylogeny.pdf"),
  width = 60,
  height = 60
)

ggtree(tree)
ggtree(tree.ml)

dev.off()

# IDENTIFY CLADES ==============================
class(tree.ml)
class(tree)
class(fit)
class(bs)

# convert distance
cophenetic_matrix <- cophenetic(tree.ml)

# hierarchical clustering
hc <- hclust(as.dist(cophenetic_matrix))

# Cluster identify
library(cluster)

for(k in 2:10){
  cl <- cutree(hc,k)
  sil <- silhouette(cl,dist(cophenetic_matrix))
  print(c(k,mean(sil[,3])))
}

# print tree
BiocManager::install("ggtree")
library(ggtree)

p <- ggtree(tree_ml)
p +
  geom_tiplab(size=2) +
  theme_tree2()
p +
  geom_tippoint(
    aes(color=factor(clades)),
    size=2
  )

# save
ggsave(
  file.path(getwd(),"plot","Spike_phylogeny.pdf"),
  width=12,
  height=10
)

# PLOT ALIGNMENT ================================
library(ggmsa)
ggmsa(
  alignment_filtered,
  start = 1,
  end = 200
)

# Conservation map
library(msa)
msaPrettyPrint(
  alignment,
  output="pdf",
  file= file.path(getwd(),"plot","alignment.pdf")
)

# LOGO plot
library(ggseqlogo)
mat <- as.matrix(alignment_filtered)

ggseqlogo(mat)

# VARIABLE SITE ANALYSIS =========================
n_sites <- ncol(mat)
variable <- logical(n_sites)

for(i in seq_len(n_sites)) {
  
  residues <- unique(
    mat[, i]
  )
  
  residues <- residues[
    residues != "-"
  ]
  
  variable[i] <- length(residues) > 1
}

# COUNT VARIATION FREQUENCY ==========
# Create the variation frequency
variation_frequency <- numeric(n_sites)
for(i in seq_len(n_sites)) {
  
  tab <- table(
    mat[, i]
  )
  
  variation_frequency[i] <-
    1 - max(tab) / sum(tab)
}

# plot 
variation_df <- data.frame(
  Position = 1:n_sites,
  Variation = variation_frequency
)

ggplot(
  variation_df,
  aes(Position, Variation)
) +
  geom_line() +
  theme_bw() +
  labs(
    title = "Spike Protein Variation Frequency",
    x = "Amino Acid Position",
    y = "Variation Frequency"
  )

ggsave(
  "Spike_variation_plot.pdf",
  width = 12,
  height = 4
)

# most variable residuos
top_sites <- variation_df %>%
  arrange(desc(Variation)) %>%
  head(50)

write.csv(
  top_sites,
  "Top_Variable_Sites.csv",
  row.names = FALSE
)