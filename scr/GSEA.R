# installl libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
BiocManager::install("hugene10sttranscriptcluster.db")
BiocManager::install("GEOquery")
BiocManager::install("GOstats")
BiocManager::install("topGO")

# I. Import database -------------------

# download GPL
GPL <- read.delim("~/Project/Bioinformatic_2026/data/GPL6244-17930.txt", comment.char="#")

# upload gene list
eb <- read.delim("~/Project/Bioinformatic_2026/data/EB_deg.txt")
eb$smb

# II. Change gen codes by same database -------

# loading BiomaRt
library(biomaRt)
# exploring command
biomaRt::listEnsembl()
biomaRt::listMarts()

# connect to database
mart <- biomaRt::useEnsembl(biomart = "genes")

# check list dataset
biomaRt::listDatasets(mart = mart)

# search dataset 
biomaRt::searchDatasets(mart = mart, pattern = "hsapiens")

# select dataset
mart <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# convert gene code to useful one
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") 

# search attributes
biomaRt::searchAttributes(mart = mart,pattern = "entrez")     # entrezgene_id 
biomaRt::searchAttributes(mart = mart,pattern = "hgnc_symbol")# symbol of the gene
biomaRt::searchAttributes(mart = mart,pattern = "name_1006") # GO term name 
biomaRt::searchAttributes(mart = mart,pattern = "go_id")     # GO id acc num

# create table with ensembl code of the gene of interest
eb.ensm <-getBM(c("entrezgene_id","hgnc_symbol"), 
               filters="hgnc_symbol", 
               values=eb$smb, 
               mart=mart)
head(eb.ensm)

# create table with ensble code of the gene in the microarrays platform  
library(hugene10sttranscriptcluster.db)
library(AnnotationDbi)
library(GEOquery)

gse <- getGEO("GSE106571")   # call the Gen set project GSE106571
eset <- gse[[1]]             # unlist the object
expr <- exprs(eset)          # extract the expression matrix
head(colnames(expr))         # Id samples are columnnames
head(rownames(expr))         # Id genes are rownames 

AnnotationDbi::columns(hugene10sttranscriptcluster.db) # available columns to match

# match id_probe with gene symbol
gpl.ensm <- select(
              hugene10sttranscriptcluster.db,   # Dbi
              keys = rownames(expr),            # Id_probe
              columns = c("SYMBOL","ENTREZID"), # Output columns
              keytype = "PROBEID"               # Input column
              )

head(head(gpl.ensm))                                  # seems empty, but...
head(sort(table(gpl.ensm$SYMBOL),decreasing = T),20)  # there are many probe to the same gene

# III. Create a GO annotated tables ---------------

head(sort(table(eb.ensm$hgnc_symbol), decreasing = T))  # check if there are duplicate genes
# we have 3 duplicate genes
# TPTEP1 ZNF117 ZNF846   AARD  AASDH  ACAA2 
# 2      2      2      1      1      1 


goannot<-getBM(c("hgnc_symbol", "entrezgene_id", "go_id", "name_1006"),   # output columns
               filters="hgnc_symbol",                      # match column
               values=eb.ensm$hgnc_symbol,                 # vector index
               mart=mart)                                  # mart object
goannot[1:15,1:3]  # check output

# second option using the Dbi 
library(org.Hs.eg.db)            # load Dbi human library
AnnotationDbi::columns(org.Hs.eg.db) # check available columns
goannot2<-select(org.Hs.eg.db,   # load reference Dbi
                 columns=c("SYMBOL","ENTREZID","GO"),   # output column
                 keys = eb.ensm$hgnc_symbol, # input column
                 keytype="SYMBOL")           # match column
head(goannot2)

# KEGG pathways
AnnotationDbi::columns(org.Hs.eg.db) # check available columns
keggannot<-select(org.Hs.eg.db,   # load reference Dbi
                 columns=c("SYMBOL","ENTREZID","PATH"),   # output column
                 keys = eb.ensm$hgnc_symbol, # input column
                 keytype="SYMBOL")           # match column
head(keggannot)

# Chromosome location
AnnotationDbi::columns(org.Hs.eg.db) # check available columns
chrmannot<-select(org.Hs.eg.db,   # load reference Dbi
                  columns=c("SYMBOL","ENTREZID","MAP"),   # output column
                  keys = eb.ensm$hgnc_symbol, # input column
                  keytype="SYMBOL")           # match column
head(chrmannot)

# create GO term table from the microarray platform  
goannot<-getBM(c("hgnc_symbol", "entrezgene_id", "go_id", "name_1006"),   # output columns
               filters="hgnc_symbol",                      # match column
               values=gpl.ensm$SYMBOL,                 # vector index
               mart=mart)                                  # mart object
goannot[1:15,1:3]  # check output



# IV. Hypergeometric test using GOstats =======
library(GOstats)
library(org.Hs.eg.db)


# Over represented GO terms: MF
# STEP 1 - prepare data
paramsGO <- new("GOHyperGParams", 
                geneIds = as.character(eb.ensm$entrezgene_id),        # target gene identifiers
                universeGeneIds = as.character(gpl.ensm$ENTREZID),    # universe of genes
                annotation = "org.Hs.eg.db",            # name of the annotation data package for the chip used to generate the dat 
                ontology = "MF",                        # BP, CC, MF
                pvalueCutoff = 0.01, 
                conditional = FALSE, 
                testDirection = "over")                 # over or under representate genes
# STEP 2 - test
hgt.Over.GO <- hyperGTest(paramsGO) 
head(summary(hgt.Over.GO))

# under represented GO terms: BP
# STEP 1 - prepare data
paramsGO <- new("GOHyperGParams", 
                geneIds = as.character(eb.ensm$entrezgene_id),        # target gene identifiers
                universeGeneIds = as.character(gpl.ensm$ENTREZID),    # universe of genes
                annotation = "org.Hs.eg.db",            # name of the annotation data package for the chip used to generate the dat 
                ontology = "BP",                        # BP, CC, MF
                pvalueCutoff = 0.01, 
                conditional = TRUE, 
                testDirection = "under")                 # over or under representate genes
# STEP 2 - test
hgt.Under.GO <- hyperGTest(paramsGO) 
head(summary(hgt.Under.GO))


# over represented KEGG pathways
# STEP 1 - prepare data
paramsKEGG <- new("KEGGHyperGParams", 
                geneIds = as.character(eb.ensm$entrezgene_id),        # target gene identifiers
                universeGeneIds = as.character(gpl.ensm$ENTREZID),    # universe of genes
                annotation = "org.Hs.eg.db",            # name of the annotation data package for the chip used to generate the dat 
                pvalueCutoff = 0.05, 
                testDirection = "over")                 # over or under representate genes
# STEP 2 - test
hgt.Under.Kegg <- hyperGTest(paramsKEGG) 
head(summary(hgt.Under.Kegg))

# V. ORA analysis using topGO =================
library(topGO)

# STEP 1: prepare the data
n <- length(unique(goannot$entrezgene_id))  # number of genes in the universe
# split GO terms by gene
gene2GO <- split(
  goannot$go_id,
  goannot$entrezgene_id
)
gene2GO$`343066`   # Check

# STEP 2: create a index of genes of interest
geneList<-rep(0,n) 
geneList[names(gene2GO) %in% eb.ensm$entrezgene_id] <- 1 
names(geneList)<-names(gene2GO) 
geneList<-factor(geneList)
head(geneList)

# STEP 3: build the 
GOdata.MF <- new("topGOdata", 
                 ontology = "MF", 
                 description= "MF on Epidermolisis Bullosa", 
                 allGenes = geneList, 
                 annot =annFUN.gene2GO, 
                 gene2GO = gene2GO)

# STEP 4: Statistical analysis
resultFisher.MF.classic <- runTest(GOdata.MF, 
                                   algorithm = "classic", 
                                   statistic = "fisher") 
(allRes.MF.classic <- GenTable(GOdata.MF, 
                              classicFisher = resultFisher.MF.classic,
                              topNodes=20) )

resultFisher.MF.weight <- runTest(GOdata.MF, 
                                  algorithm = "weight", 
                                  statistic = "fisher") 
(allRes.MF.weight <- GenTable(GOdata.MF, 
                             classicFisher = resultFisher.MF.weight,
                             topNodes=20)) 

resultFisher.MF.parentchild <- runTest(GOdata.MF, 
                                       algorithm = "parentchild", 
                                       statistic = "fisher")
(allRes.MF.parentchild <- GenTable(GOdata.MF, 
                                  classicFisher = resultFisher.MF.parentchild,
                                  topNodes=20))




