
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
organism = "org.Mm.eg.db"
library(organism, character.only = TRUE)

Bulk_De <- read.csv('DE_VP16_EX8_T_Cell.csv')
df = Bulk_De

# we want the log2 fold change 
original_gene_list <- df$avg_logFC

# name the vector
names(original_gene_list) <- df$X
df$X
# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 <-  df[df$X %in% dedup_ids$SYMBOL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y <-  dedup_ids$ENTREZID

# Create a vector of the gene universe
kegg_gene_list <- df2$avg_logFC

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)


kegg_organism = "mmu"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               minGSSize    = 1,
               maxGSSize = 5000,
               pvalueCutoff = .8,
               pAdjustMethod = "BH",
               keyType       = "ncbi-geneid")

ego3 <- gseGO(geneList     = kegg_gene_list,
              OrgDb        = org.Mm.eg.db,
              ont          = "CC",
              nPerm        = 1000,
              minGSSize    = 1,
              maxGSSize    = 5000,
              pvalueCutoff = 1,
              verbose      = FALSE)

dotplot(kk2, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

emapplot(kk2)
upsetplot(kk2)

dotplot(ego3, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

library(msigdbr)
msigdbr_show_species()

m_t2g <- msigdbr(species = "Mus musculus", category = "C2")
head(m_t2g, 2) %>% as.data.frame
egmt2 <- GSEA(gene_list, TERM2GENE = m_t2g)

# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(kk2, foldChange=gene_list, circular = TRUE, colorEdge = TRUE)
cnetplot(kk2, categorySize="pvalue", foldChange=gene_list)

cnetplot(kk2, categorySize="pvalue", foldChange=gene_list)
library(pathview)

# Produce the native KEGG plot (PNG)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="mmu04310", species = kegg_organism)