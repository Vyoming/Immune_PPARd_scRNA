
# take raw data and normalise it
count_raw <- Immune_VP16@assays$RNA@counts
count_norm <- apply(count_raw, 2, function(x)
  (x/sum(x))*10000)
write.table(matrix1, 'cellphonedb_count.txt', sep='\t', quote=F)
# generating meta file
meta_data <- cbind(rownames(Immune_VP16@meta.data), Immune_VP16@meta.data[,'Cell_Type', drop=F]) 
# cluster is the user's corresponding cluster column
write.table(meta_data, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)


rownames(count_norm)



require(EWCE)
require(tibble)
require(biomaRt)
require(tidyr)
require(dplyr)
# Basic function to convert mouse to human gene names
alldata <- Immune_VP16
alldata <- NormalizeData(alldata, normalization.method = "LogNormalize", scale.factor = 10000)
allgenes <- rownames(alldata)
matrix1 <- as.data.frame(alldata@assays$RNA@data)
matrix1 <- matrix1[rowSums(matrix1[,2:dim(matrix1)[2]])!=0,]

### If you are using a mouse data, then its needed to convert the gene names to human orthologs
human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(alldata@assays$RNA@data) , mart = mouse, attributesL = c("hgnc_symbol","hgnc_id",'ensembl_gene_id'), martL = human, uniqueRows=T)
matrix1 <- matrix1[match(genesV2$MGI.symbol,rownames(alldata),nomatch=F),]
genesV2 <- genesV2[match(genesV2$MGI.symbol,rownames(alldata),nomatch=F),]

matrix1$gene <- genesV2$Gene.stable.ID

#rownames(matrix1) <-  genesV2$Gene.stable.ID

### Subseting the matrix
s1 <- grepl('VP16PPARd',alldata@meta.data$Type)
s2 <- grepl('Control',alldata@meta.data$Type)
s1[match('gene',colnames(matrix1))] <- TRUE
s2[match('gene',colnames(matrix1))] <- TRUE

## Checking the dimensions
print(dim(matrix1[,s1]))
print(dim(matrix1[,s2]))


## If the cluster names are categorical, you will need to convert it to numerical
alldata@meta.data$Cell_Type <- as.factor(alldata@meta.data$Cell_Type)
print(levels(alldata@meta.data$Cell_Type))
levels(alldata@meta.data$Cell_Type) <- 1:length(levels(alldata@meta.data$Cell_Type))
print(1:length(levels(alldata@meta.data$Cell_Type)))
alldata@meta.data$Cell_Type <- as.numeric(alldata@meta.data$Cell_Type)

write.table(matrix1, 'cellphonedb_count.csv',row.names=T,sep=',')
write.table(matrix1[,s1], 'Vp16_filtered_hcount.csv',row.names=T,sep=',')
write.table(matrix1[,s2], 'Contr_filtered_hcount.csv',row.names=T,sep=',')
metadata <- data.frame(cells=rownames(alldata@meta.data[grepl('state1',alldata@meta.data$stim),]),cluster=alldata@meta.data$Cell_Type[grepl('state1',alldata@meta.data$stim)])
metadata_s2 <- data.frame(cells=rownames(alldata@meta.data[!grepl('state1',alldata@meta.data$stim),]),cluster=alldata@meta.data$Cell_Type[!grepl('state1',alldata@meta.data$stim)]) ## Just negate grepl('state1',alldata@meta.data$stim),]
print('Writing Metadata')
write.csv(metadata, 's1_filtered_meta.csv', row.names=FALSE)
write.csv(metadata_tac, 's2_filtered_meta.csv', row.names=FALSE)```