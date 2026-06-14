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
# cut names
names(alignment2) <- substr(names(alignment2), 1, 8)

# transform data
mat <- as.matrix(alignment2)
library(phangorn)
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
tree.ml <- midpoint(tree.ml)

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
# convert distance
cophenetic_matrix <- cophenetic(tree.ml)
heatmap(cophenetic_matrix)

# hierarchical clustering
hc <- hclust(as.dist(cophenetic_matrix))
plot(hc)

# Cluster identify (clade)
library(cluster)
res1 <- numeric(length = length(2:25))
names(res1) <- as.character(2:25)
for(k in 2:25){
  cl <- cutree(hc,k)
  sil <- silhouette(cl,dist(cophenetic_matrix))
  res1[k-1] <- mean(sil[,3])
}
plot(res1,type="l")


# mark the 4 clades
groups <- cutree(hc, k = 15)
table(groups)

#create a metadata
library(dplyr)
meta <- data.frame(
  label = names(groups),
  Clade = factor(groups)
)
head(meta)

# save meta
write.csv(meta,
          file.path(getwd(),"result","spike_clades.csv"),
          row.names = FALSE)

# save metadata for iTOL external tool
library(readr)
library(dplyr)

# print tree
library(ggtree)

p <- ggtree(tree.ml)
p <- p %<+% meta
p +
  geom_tiplab(size=2) +
  theme_tree2()
p +
  geom_tippoint(aes(color = Clade), size = 2) +
  theme_tree2()

ggtree(tree.ml, layout = "circular") %<+% meta +
  geom_tippoint(aes(color = Clade), size = 1.5)

# save
ggsave(
  file.path(getwd(),"plot","Spike_phylogeny.pdf"),
  width=12,
  height=10
)

# PLOT ALIGNMENT ================================
library(ggmsa)
ggmsa(
  alignment2,
  start = 1,
  end = 20
)

# # Conservation map
# library(msa)
# msaPrettyPrint(
#   alignment,
#   output="pdf",
#   file= file.path(getwd(),"plot","alignment.pdf")
# )


# EXTRACT CLADES SEQUENCE =====================================

consensus_list <- lapply(sort(unique(groups)), function(k){
  seq_names <- names(groups[groups == k])
  seqs <- alignment2[seq_names]
  consensusString(seqs)
})

# remove characters "x", "?", or "-"
consensus_list <- lapply(
  consensus_list,
  function(x) gsub("[X?\\-]", "", x)
)

names(consensus_list) <- paste0("Clade_", sort(unique(groups)))

library(Biostrings)

consensusAA <- AAStringSet(
  lapply(consensus_list, AAString)
)

# IMPORT REFERENCE SEQUENCE ===========================
library(Biostrings)
library(rentrez)
# Download reference spike sequences
fasta <- entrez_fetch(
  db = "protein",
  id = "YP_009724390.1",  # "Wuham Spike"
  rettype = "fasta",
  retmode = "text",
)

# save refseq
outfile <- file.path(getwd(),"result","refSeq.fasta")
cat(
  fasta,
  file = outfile,
  append = TRUE
)
RefSeqWuham <- readAAStringSet(outfile)
names(RefSeqWuham) <- "WuhamREFSEQ"

# Join refseq
length(consensusAA)
consensusAA <- c(RefSeqWuham,consensusAA)
length(consensusAA)

# ALIGNM SEQUENCES ================================
library(msa)
alignment3 <- msa(
  consensusAA,
  method = "Muscle", # Omega more efficient for >1000 seq
  verbose = TRUE
)

# Evaluate alignment quality
library(DECIPHER)
alignment4 <- AAStringSet(alignment3)
BrowseSeqs(alignment4)


# VARIABLE SITE ANALYSIS =========================
# convert a matrix
aln_mat <- as.matrix(alignment4)

# compute entrophy
site_entropy <- apply(
  aln_mat,
  2,
  function(x){p <- table(x)/length(x);-sum(p*log2(p))}
)
plot(site_entropy,type="l")
plot(smooth(site_entropy),type="l")
site_entropy <- smooth(site_entropy)

# define variable position
var.pos <- which(site_entropy>0.5)
var.pos <- var.pos[var.pos>100&var.pos<1200] # remove extreme

# 247 248 249 412 413 414 415 483 484 515 516 521 522
# 523 524 525 526 719 720 721 722 723

var.pos <- list(c(240:260),
                c(405:425),
                c(475:495),
                c(505:525),
                c(515:535),
                c(710:730))

# Plot
# install.packages("ggseqlogo")
library(ggseqlogo)
library(ggplot2)

s1 <- apply(aln_mat[,var.pos[[1]]],
            1,
            paste0,
            collapse = "")
ggseqlogo(
  s1,
  method = "probability",
  seq_type = "aa"
)
px <- list()
for (i in 1:length(var.pos)) {
  tmp <- var.pos[[i]]
  px[[i]] <- ggseqlogo(
    apply(aln_mat[,tmp],1,paste0,collapse = ""),
    method = "probability",
    seq_type = "aa"
    ) +
    labs(title = paste("POSITION: ",
                       tmp[1],
                       ":",
                       tmp[length(tmp)],
                       sep=""),
         )
}

names(px) <- paste("p",1:6,sep = "")

# install.packages("patchwork")
library(patchwork)
pdf(file=file.path(getwd(),"plot","logo_Spike.pdf"),
    width = 40,height = 20)
(px$p1+px$p2)/(px$p3+px$p4)/(px$p5+px$p6)
dev.off()

# PHYLOGENY ==============================
# transform data
mat <- as.matrix(alignment4)
library(phangorn)
phy <- phyDat(
  mat,
  type = "AA"
)

# calculate distance matrix
dm <- dist.ml(phy)

# Neightbor join tree
tree2 <- NJ(dm)
tree2 <- midpoint(tree2) # root the midpoint

# Optimize tree
fit <- pml(tree2,phy)
fit <- optim.pml(
  fit,
  model="JTT",
  optGamma=TRUE,
  optInv=TRUE
)

# best tree
tree2.ml <- fit$tree
tree2.ml <- midpoint(tree2.ml)

# plot
library(ggtree)
ggtree(tree2)
ggtree(tree2.ml) +
  geom_tiplab(align = T) +
  #geom_text(aes(label=node), hjust=-.3)
  geom_cladelabel(node=23, label="C-21", 
                  color="red2", offset=.01) + 
  geom_cladelabel(node=27, label="C-22", 
                  color="blue", offset=.025) +
  geom_cladelabel(node=29, label="C-23", 
                  color="orange", offset=.01) 

# 23, 27, 29
# save
write.tree(            # save the tree
  tree2.ml,
  file = file.path(getwd(),"result","HL_Spike_tree2.nwk")
)
