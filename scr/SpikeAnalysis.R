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
  total = length(seq(0, 1900, by = 50)),
  clear = FALSE,
  width = 90
)

# Create address of output directory
outfile <- file.path(getwd(), 
                     "data", 
                     "Spike2000.fasta")


for(start in seq(0, 1900, by = 50)) {
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
p <- boxplot(len$x)
p$out

summary(len$x)

seqs <- seqs[
  len$x > 1220 &
    len$x < 1320
]

length(seqs)
len <- data.frame(x=width(seqs))
summary(len$x)
ggplot( len, aes(x=x)) +
  geom_histogram( binwidth=3, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("SPIKE") +
  labs(subtitle = "(Sequence legnths)") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  )

# SEQUENCE ALIGNMENT ===================
alignment <- AlignSeqs(
  seqs,
  processors = NULL
)
writeXStringSet(
  alignment,
  filepath = file.path(getwd(),"result","Spike_alignment.fasta")
)

# CONCESUS SEQUENCE =============
consensus <- ConsensusSequence(
  alignment
)

writeLines(
  as.character(consensus),
  file.path(getwd(),"result","Spike_consensus.fasta")
)

# PHYLOGENY ==============================

# transform data
mat <- as.matrix(alignment)

phy <- phyDat(
               mat,
                type = "AA"
                )

# calculate distance matrix
dm <- dist.ml(phy)

# Neightbor join tree
tree <- NJ(dm)
tree <- midpoint(tree) # root the midpoint
write.tree(            # save the tree
  tree,
  file = "Spike_tree.nwk"
)

# plot phylogeny tree
pdf(
  "Spike_phylogeny.pdf",
  width = 14,
  height = 14
)

plot(
  tree,
  cex = 0.15,
  no.margin = TRUE
)

dev.off()

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