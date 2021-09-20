# Tumorigenic Immune SCRNA analysis
# Vyom Shah - Beyaz Laboratory
library(CIPR)
library(Rmagic)
library(nichenetr)
library(Seurat)
library(tidyverse)
library(cowplot)
library(ggplot2)
library(MAST)
library(DESeq2)
library(EnhancedVolcano)
library(limma)
library(scales)
library(metR)
library(ggpubr)
library(rstatix)
library(svglite)
library(viridis)


# Load in the Data
Control <- Read10X(data.dir = "/Users/vyom/data/tumor_infiltrating_lymphocytes/IE01_control/outs/filtered_feature_bc_matrix/")
CPT1aKO <- Read10X(data.dir = "/Users/vyom/data/tumor_infiltrating_lymphocytes/IE01_CPT1aKO/outs/filtered_feature_bc_matrix/")
VP16PPARd <- Read10X(data.dir = "/Users/vyom/data/tumor_infiltrating_lymphocytes/IE01_VP16PPARd/outs/filtered_feature_bc_matrix/")

Control <- CreateSeuratObject(Control, project = "Control")
CPT1aKO <- CreateSeuratObject(CPT1aKO, project = "CPT1aKO")
VP16PPARd <- CreateSeuratObject(VP16PPARd, project = "VP16PPARd")

d <- merge(Control, y = c(CPT1aKO,VP16PPARd), add.cell.ids = c("Control", "CPT1aKO","VP16PPARd"), project = "Immune")
d
d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "mt-")
VlnPlot(d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0 ,ncol = 3)
d <- subset(d, subset = nCount_RNA > 100 & nCount_RNA < 20000 & nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 8)
d
d$orig.ident
Data.list <- SplitObject(d, split.by = "ident")
Data.list <- Data.list[c("Control", "CPT1aKO","VP16PPARd")]
for (i in 1:length(Data.list)) {
  
  Data.list[[i]] <- SCTransform(Data.list[[i]], verbose = FALSE)
}

# Normilization
#select highly variable genes 
Data.features <- SelectIntegrationFeatures(object.list = Data.list, nfeatures = 5000)
options (future.globals.maxSize = 4000 * 1024^5)
Data.list <- PrepSCTIntegration(object.list = Data.list, anchor.features = Data.features, 
                                verbose = FALSE)
Data.anchors <- FindIntegrationAnchors(object.list = Data.list, normalization.method = "SCT", 
                                       anchor.features = Data.features, verbose = FALSE)
rm(Data.list, d)
Immune_obj <- IntegrateData(anchorset = Data.anchors, normalization.method = "SCT", 
                          verbose = TRUE)
rm(Data.anchors)
# Visulization and Clustering
Immune_obj <- RunPCA(Immune_obj)
VizDimLoadings(Immune_obj, dims = 1:2, reduction = "pca")

DimPlot(Immune_obj, reduction = "pca")
ElbowPlot(Immune_obj, ndims = 50, reduction = "pca")

Immune_obj <- RunUMAP(Immune_obj, dims = 1:10)

Idents(Immune_obj) <- levels(factor(Immune_obj@meta.data$orig.ident))
Immune_obj[["Type"]] <- Immune_obj$orig.ident
Immune_obj$Type
Immune_obj@meta.data$Type <- factor(Immune_obj@meta.data$Type, levels = c("Control", "CPT1aKO", "VP16PPARd")) 

Immune_obj <- FindNeighbors(Immune_obj, dims = 1:10)
Immune_obj <- FindClusters(Immune_obj, resolution = 1)
DimPlot(Immune_obj, reduction = "umap", label = TRUE)


markers <- FindAllMarkers(Immune_obj, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 1)
markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(markers,'Immutumor_Sigs_Per_Clust_anno.csv')


allmarkers <- FindAllMarkers(Immune_obj, min.pct = 0.5)
# cluster genes using CIPR
CIPR(input_dat = allmarkers,
     comp_method = "logfc_dot_product", 
     reference = "immgen", 
     plot_ind = T,
     plot_top = F)


FeaturePlot(object = Immune_obj, features = c('Plac8'),split.by = ,pt.size = 1)  + scale_color_viridis(option = 'B')
DotPlot(Immune_obj, features = c('C1qa','Foxp3'), group.by = 'integrated_snn_res.1')

VlnPlot(Immune_obj, group.by = 'Cell_Type' , ncol = 1, split.by = "Type", features = c('Gzmb'), pt.size = 0) + theme(legend.position = 'none') + geom_boxplot(width=0.5,outlier.shape = NA, coef = 0) + xlab('')
macrophage_markers <- c('Slc7a2',	'Serpinb2',	'Ppap2a',	'AA467197',	'Slc7a11',	'Al504432',	'Il4i1',	'Cd40',	'AI504432',	'F10',	'Rasgrp1',	'Jak2',	'Malt1',	'Mmp13',	'Upp1',	'Itgb8',	'Syt7',	'Cst7',	'H2-Ab1',	'Rcl1',	'Ch25h',	'Ccl22',	'Socs1',	'Aqp9',	'Slamf6',	'E430024C06Rik',	'Fmnl2',	'Adhfe1',	'Adora2b',	'Pde4b',	'St3gal1',	'Plekhg1',	'Ptgir',	'Ier3',	'Cav1',	'Trib3')
Immune_obj <- AddModuleScore(object = Immune_obj, features = list(macrophage_markers), name = 'macrophage_markers')



Prop_table<- prop.table(x = table(Immune_obj$integrated_snn_res.1, Immune_obj$Type), margin = 2)

#______________________________Clustering________________________________________#
Immune_obj <- FindClusters(Immune_obj, resolution = 1)
new.cluster.ids <- c('Macrophage','Macrophage','Macrophage','Macrophage','Macrophage','ILC','Neutrophil','Macrophage','Macrophage','Macrophage','Macrophage','T Cell','Macrophage','NK Cell','NK Cell','T Cell','DC','Macrophage','T Cell','B Cell','T Cell','Mast Cell','T Cell')
Immune_obj[["Cell_type"]] <- Idents(Immune_obj)
names(new.cluster.ids) <- levels(Immune_obj)
Immune_obj <- RenameIdents(Immune_obj, new.cluster.ids)
Immune_obj[["Cell_type"]] <- Idents(Immune_obj)
DimPlot(Immune_obj, reduction = "umap", group.by= 'Cell_type')
DimPlot(Immune_obj, reduction = "umap", group.by= 'Type')

#proportions Control vs All
Prop_table<- prop.table(x = table(Immune_obj$Cell_type, Immune_obj$Type), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table1 <- Prop_Table
my_levels <- c('T Cell','B Cell','NK Cell','Neutrophil','Macrophage','Monocyte','DC','Mast Cell','ILC')
Prop_Table1$Var1 <- factor(Prop_Table1$Var1,levels = my_levels)
Prop_Table1$Var2 <- factor(Prop_Table1$Var2,levels = c('Control','CPT1aKO','VP16PPARd'))
plot <- ggplot(data = Prop_Table1, aes(Var1, Freq, fill=Var2)) + geom_bar(position="dodge", stat="identity", na.rm = TRUE) +  scale_fill_manual(values = c('#1b9e77' ,'#d95f02', '#7570b3')) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot

#proportions Control vs VP16
Idents(Immune_obj) <- Immune_obj$Type
Immune_VP16 <- subset(Immune_obj,  idents =c('Control', 'VP16PPARd'))
my_levels <- c('Control' ,'VP16PPARd')
Immune_VP16$Type <- factor(x = Immune_VP16$Type, levels = my_levels)

Prop_table<- prop.table(x = table(Immune_VP16$Cell_type, Immune_VP16$Type), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table1 <- Prop_Table
my_levels <- c('T Cell','B Cell','NK Cell','Neutrophil','Macrophage','Monocyte','DC','Mast Cell','ILC')
Prop_Table1$Var1 <- factor(Prop_Table1$Var1,levels = my_levels)
Prop_Table1$Var2 <- factor(Prop_Table1$Var2,levels = c('Control','VP16PPARd'))
plot <- ggplot(data = Prop_Table1, aes(Var1, Freq, fill=Var2)) + geom_bar(position="dodge", stat="identity", na.rm = TRUE) +  scale_fill_manual(values = c('#1b9e77' ,'#7570b3')) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot

#proportions Control vs CPT1aKO
Idents(Immune_obj) <- Immune_obj$Type
Immune_CPT1a <- subset(Immune_obj,  idents =c('Control', 'CPT1aKO'))
my_levels <- c('Control' ,'CPT1aKO')
Immune_CPT1a$Type <- factor(x = Immune_CPT1a$Type, levels = my_levels)

Prop_table<- prop.table(x = table(Immune_CPT1a$Cell_type, Immune_CPT1a$Type), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table1 <- Prop_Table
my_levels <- c('T Cell','B Cell','NK Cell','Neutrophil','Macrophage','Monocyte','DC','Mast Cell','ILC')
Prop_Table1$Var1 <- factor(Prop_Table1$Var1,levels = my_levels)
Prop_Table1$Var2 <- factor(Prop_Table1$Var2,levels = c('Control','CPT1aKO'))
plot <- ggplot(data = Prop_Table1, aes(Var1, Freq, fill=Var2)) + geom_bar(position="dodge", stat="identity", na.rm = TRUE) +  scale_fill_manual(values = c('#1b9e77' ,'#7570b3')) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot


# Further Subcluster Macrophages and T Cells

Idents(Immune_obj) <- Immune_obj$Cell_type
Immune_Macrophage <- subset(Immune_obj, idents = 'Macrophage')

Immune_Macrophage <- RunPCA(Immune_Macrophage, verbose = FALSE)

Immune_Macrophage <- RunUMAP(Immune_Macrophage, dims = 1:10)

Immune_Macrophage <- FindNeighbors(Immune_Macrophage, dims = 1:10)
Immune_Macrophage <- FindClusters(Immune_Macrophage, resolution = 2 )
DimPlot(Immune_Macrophage, reduction = "umap", label = TRUE)

markers <- FindAllMarkers(Immune_Macrophage, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 1)
markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

new.cluster.ids <- c('M1 Macrophage','Monocyte','Monocyte','M1 Macrophage','Monocyte','Monocyte','M2 Macrophage','M2 Macrophage','M2 Macrophage','Monocyte','M1 Macrophage','M2 Macrophage','M1 Macrophage','Monocyte','M1 Macrophage','Monocyte','Monocyte','M1 Macrophage','Monocyte','M2 Macrophage','M2 Macrophage','M1 Macrophage','M1 Macrophage','M1 Macrophage','M1 Macrophage','M1 Macrophage','Monocyte','M2 Macrophage','M1 Macrophage')
Immune_Macrophage[["Cell_Type"]] <- Idents(Immune_Macrophage)
names(new.cluster.ids) <- levels(Immune_Macrophage)
Immune_Macrophage <- RenameIdents(Immune_Macrophage, new.cluster.ids)
Immune_Macrophage[["Cell_Type"]] <- Idents(Immune_Macrophage)
DimPlot(Immune_Macrophage, reduction = "umap", label = TRUE)


Prop_table<- prop.table(x = table(Immune_Macrophage$Cell_Type, Immune_Macrophage$Type), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table1 <- Prop_Table
Prop_Table1$Var1 <- factor(Prop_Table1$Var1)
Prop_Table1$Var2 <- factor(Prop_Table1$Var2,levels = c('Control','CPT1aKO','VP16PPARd'))
plot <- ggplot(data = Prop_Table1, aes(Var1, Freq, fill=Var2)) + geom_bar(position="dodge", stat="identity", na.rm = TRUE) +  scale_fill_manual(values = c('#1b9e77' ,'#d95f02', '#7570b3')) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot


#further subscluster T Cells
Idents(Immune_obj) <- Immune_obj$Cell_type
Immune_T_cell <- subset(Immune_obj, idents = 'T Cell')

Immune_T_cell <- RunPCA(Immune_T_cell, verbose = FALSE)

Immune_T_cell <- RunUMAP(Immune_T_cell, dims = 1:10)

Immune_T_cell <- FindNeighbors(Immune_T_cell, dims = 1:10)
Immune_T_cell <- FindClusters(Immune_T_cell, resolution = 2)
DimPlot(Immune_T_cell, reduction = "umap", label = TRUE)

markers <- FindAllMarkers(Immune_T_cell, only.pos = TRUE, min.pct = 0.5, logfc.threshold = .5)
library(readxl)
write.csv(markers,'Gene_de_Sigs_Per_Clust_immune.csv')

DotPlot(Immune_obj, group.by = 'integrated_snn_res.1', features = c('Cd3g'))
Idents(Immune_obj) <-  Immune_obj$integrated_snn_res.1
markers <- FindAllMarkers(Immune_obj, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 1)

Prop_table<- prop.table(x = table(Immune_T_cell$integrated_snn_res.2, Immune_T_cell$Type), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table1 <- Prop_Table
Prop_Table1$Var1 <- factor(Prop_Table1$Var1)
Prop_Table1$Var2 <- factor(Prop_Table1$Var2,levels = c('Control','CPT1aKO','VP16PPARd'))
plot <- ggplot(data = Prop_Table1, aes(Var1, Freq, fill=Var2)) + geom_bar(position="dodge", stat="identity", na.rm = TRUE) +  scale_fill_manual(values = c('#1b9e77' ,'#d95f02', '#7570b3')) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot

new.cluster.ids <- c('TReg','CD8+ T Cell','GD T lymphocytes','CD8+ T Cell','Ex CD8+ T Cell','NKT Cell','Ex CD8+ T Cell','TReg','Ex CD8+ T Cell','TReg','Ex CD8+ T Cell','CD4+ T Cell(Th2)','CD4+ T Cell(Th2)','CD8+ T Cell','CD8+ T Cell','CD8+ T Cell','CD8+ T Cell','GD T lymphocytes','Proliferating TReg','Stromal Cell','CD8+ T Cell','CD4+ T Cell(Th17)')

Immune_T_cell[["Cell_Type"]] <- Idents(Immune_T_cell)
names(new.cluster.ids) <- levels(Immune_T_cell)
Immune_T_cell <- RenameIdents(Immune_T_cell, new.cluster.ids)
Immune_T_cell[["Cell_Type"]] <- Idents(Immune_T_cell)
DimPlot(Immune_T_cell, reduction = "umap", label = TRUE)

Prop_table<- prop.table(x = table(Immune_T_cell$integrated_snn_res.2, Immune_T_cell$Type), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table1 <- Prop_Table
Prop_Table1$Var1 <- factor(Prop_Table1$Var1)
Prop_Table1$Var2 <- factor(Prop_Table1$Var2,levels = c('Control','CPT1aKO','VP16PPARd'))
plot <- ggplot(data = Prop_Table1, aes(Var1, Freq, fill=Var2)) + geom_bar(position="dodge", stat="identity", na.rm = TRUE) +  scale_fill_manual(values = c('#1b9e77' ,'#d95f02', '#7570b3')) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot

#recluster based on Macrophages and T cells
DefaultAssay(Immune_obj) <-  "integrated"
Immune_obj <- FindClusters(Immune_obj, resolution = 1)
new.cluster.ids <- c('','','','','','MDSC','Neutrophil','','','','','','','NK Cell','NK Cell','','DC','','','B Cell','','M1 Macrophage','')
Immune_obj[["Cell_type"]] <- Idents(Immune_obj)
names(new.cluster.ids) <- levels(Immune_obj)
Immune_obj <- RenameIdents(Immune_obj, new.cluster.ids)
Immune_obj[["Cell_type"]] <- Idents(Immune_obj)

DimPlot(Immune_obj, reduction = "umap", group.by= 'Cell_type')
DimPlot(Immune_obj, reduction = "umap", group.by= 'Type')
Immune_obj$Cell_type <- as.character(Idents(Immune_obj))

Immune_obj$Cell_type[WhichCells(Immune_Macrophage)] <- paste(Idents(Immune_Macrophage))
Immune_obj$Cell_type[WhichCells(Immune_T_cell)] <- paste(Idents(Immune_T_cell))
Idents(Immune_obj) <- Immune_obj$Cell_type
DimPlot(Immune_obj, reduction = "umap", group.by= 'Cell_type')


VlnPlot(Immune_obj, group.by = 'Cell_Type' , ncol = 4, split.by = "Type", features = c('Cpt1a', 'Havcr2', 'Mrc1','Tox', 'Cx3cr1','Pdcd1','Lag3','Batf'), pt.size = 0, assay = "RNA", cols = c('#1b9e77' ,'#d95f02', 'blue')) + theme(legend.position = 'none') + geom_boxplot(width=0.5,outlier.shape = NA, coef = 0) + xlab('')
DotPlot(Immune_obj ,features = c('Cpt1a', 'Havcr2', 'Mrc1','Tox', 'Cx3cr1','Pdcd1','Lag3','Batf'), assay = 'RNA', split.by = 'Type', cols = c('#1b9e77','#d95f02','blue')) 
VlnPlot(Immune_obj, group.by = 'Cell_Type' , ncol = 1, split.by = "Type", features = c('Cpt1a'), pt.size = 0, assay = "RNA", cols = c('#1b9e77' ,'#d95f02', 'blue'))  + geom_boxplot(width=0.5,outlier.shape = NA, coef = 0) + xlab('')
DotPlot(Immune_T_cell ,features = c('Ldha'), assay = 'RNA', split.by = 'Type', cols = c('#1b9e77','#d95f02','blue')) 
FeaturePlot(Immune_obj, features = c('Cd3g', 'Cd4','Cd8a', 'Foxp3', 'Tox','Cxcr3'), cols = c("#313695",'#ffffbf' ,"#a50026"), pt.size = 1.5)
DotPlot(Immune_T_cell, features = c('Tcra','Tcrb'), cols = c("#313695",'#ffffbf' ,"#a50026"), assay = 'RNA')
FeaturePlot(Immune_T_cell, features = c('Cx3cr1', 'Havcr2','Cd3e','Trac','Batf', 'Lag3'), cols = c("#313695",'#ffffbf' ,"#a50026"), pt.size = 1.5)

DimPlot(Immune_obj, group.by = 'integrated_snn_res.0.5')
DimPlot(Immune_obj, group.by = 'Cell_Type')

Prop_table<- prop.table(x = table(Immune_obj$Cell_Type, Immune_obj$Type), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table1 <- Prop_Table
Prop_Table1$Var1 <- factor(Prop_Table1$Var1)
Prop_Table1$Var2 <- factor(Prop_Table1$Var2,levels = c('Control','CPT1aKO','VP16PPARd'))
plot <- ggplot(data = Prop_Table1, aes(Var1, Freq, fill=Var2)) + geom_bar(position="dodge", stat="identity", na.rm = TRUE) +  scale_fill_manual(values = c('#1b9e77' ,'#d95f02', '#7570b3')) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot

#cluster sanity check
Idents(Immune_obj) <- Immune_obj$Cell_type
Immune_obj$Cell_Type <- Immune_obj$Cell_type
my_levels <- c('CD8+ T Cell','Ex CD8+ T Cell','GD T lymphocytes','CD4+ T Cell(Th17)','CD4+ T Cell(Th2)','TReg','Proliferating TReg','NKT Cell','B Cell','NK Cell','M1 Macrophage','M2 Macrophage','Monocyte','MDSC','DC','Neutrophil','Stromal Cell')
Immune_obj$Cell_Type <- factor(x = Immune_obj$Cell_type, levels = my_levels)
DimPlot(Immune_obj, group.by = "Cell_Type", label = FALSE, pt.size=1, label.size = 1)


DotPlot_Sig <- unique(c('Cd3g','Cd3e','Cd8a','Cd4','Trac','Tcrg-C1','Lag3','Pdcd1','Havcr2','Tox','Tcf7','Gzmb','Tbx21','Ifng','Gata3','Il5','Il13','Rorc','Il17a','Il17f','Foxp3','Il10','Il2rb','Il2ra','Klrd1','Cd19','Ighm','Ighg1','Cd74','Ciita','Nrc1','Klre1','Itgam','Itgax','H2-Eb1','H2-Ab1','Arg1','Mrc1','Tgfbi','Ccr2','Vegfa','Prdx1','Clec4d','Ccl5','Cd83','Ccr7','Fcn1','Msrb1','Ly6g','Col3a1','Sparc'))
DotPlot(Immune_obj, features = DotPlot_Sig, assay = 'RNA', group.by = 'Cell_Type') + labs(y= "Cell type", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 5)) +
  theme(text = element_text(size=5), axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust= .01), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())

cluster.markers <- FindAllMarkers(Immune_obj, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.01, assay = 'SCT')
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
check_heatmap <- DoHeatmap(Immune_obj, features = top10$gene, group.by = 'Cell_Type', assay = 'SCT',size = 2 , angle = 60) + NoLegend() + scale_fill_viridis(option="inferno")#top10$Gene
ggsave(file = 'Immune_cluster_Heatmap.pdf', plot=check_heatmap, width=10, height=10)
write.csv(cluster.markers,'Immune_cluster_markers.csv')

#Gene level and DE testing
DefaultAssay(Immune_obj) <- "RNA"
Immune_obj <- NormalizeData(object = Immune_obj, normalization.method = "LogNormalize", assay = "RNA")
Immune_obj <- magic(Immune_obj)
DimPlot(Immune_obj)

#differential expression anaylsis
Idents(Immune_obj) <- Immune_obj$Type
DE_VP16_Control <- FindMarkers(Immune_obj, ident.1 = "VP16PPARd", ident.2 = "Control", test.use = "MAST", logfc.threshold = .00, min.pct = 0, assay = 'RNA')
EnhancedVolcano(DE_VP16_Control, lab = rownames(DE_VP16_Control), x = 'avg_logFC', y = 'p_val_adj', title = 'VP16 vs Control', pCutoff = 10e-15, FCcutoff = 0.25, xlim = c(-1.5,1.5), ylim = c(0,320), subtitle = 'All Cells' )
write.csv(DE_VP16_Control,'DE_VP16.csv')

DimPlot(Immune_obj, group.by = 'Cell_Type')

Idents(Immune_obj) <-  Immune_obj$Cell_Type
Immune_subset <- subset(Immune_obj,  idents =c('CD8+ T Cell', 'Ex CD8+ T Cell', 'CD4+ T Cell(Th2)', 'CD4+ T Cell(Th17)','TReg','Proliferating TReg', 'NKT Cells','GD T lymphocytes'))
Idents(Immune_subset) <-  Immune_subset$Type
DE_VP16_all_T_Cell <- FindMarkers(Immune_subset, ident.1 = "VP16PPARd", ident.2 = "Control", test.use = "MAST", logfc.threshold = .001, min.pct = 0.0, assay = 'RNA')
write.csv(DE_VP16_all_T_Cell,'DE_VP16_all_T_Cell.csv')

Idents(Immune_obj) <-  Immune_obj$Cell_Type
Immune_subset <- subset(Immune_obj,  idents =c('CD8+ T Cell', 'Ex CD8+ T Cell'))
Idents(Immune_subset) <-  Immune_subset$Type
DE_VP16_all_T_Cell <- FindMarkers(Immune_subset, ident.1 = "VP16PPARd", ident.2 = "Control", test.use = "MAST", logfc.threshold = .001, min.pct = 0.0, assay = 'RNA')
write.csv(DE_VP16_all_8T_Cell,'DE_VP16_all_8T_Cell.csv')

Idents(Immune_obj) <-  Immune_obj$Cell_Type
Immune_subset <- subset(Immune_obj,  idents =c('Ex CD8+ T Cell'))
Idents(Immune_subset) <-  Immune_subset$Type
DE_VP16_EX8_T_Cell <- FindMarkers(Immune_subset, ident.1 = "VP16PPARd", ident.2 = "Control", test.use = "MAST", logfc.threshold = .001, min.pct = 0, assay = 'RNA')
write.csv(DE_VP16_EX8_T_Cell,'DE_VP16_EX8_T_Cell.csv')

Idents(Immune_obj) <-  Immune_obj$Cell_Type
Immune_subset <- subset(Immune_obj,  idents =c('CD8+ T Cell'))
Idents(Immune_subset) <-  Immune_subset$Type
DE_VP16_Cd8_T_Cell <- FindMarkers(Immune_subset, ident.1 = "VP16PPARd", ident.2 = "Control", test.use = "MAST", logfc.threshold = .001, min.pct = 0, assay = 'RNA')
write.csv(DE_VP16_Cd8_T_Cell,'DE_VP16_Cd8_T_Cell.csv')

Idents(Immune_obj) <-  Immune_obj$Cell_Type
Immune_subset <- subset(Immune_obj,  idents =c('CD4+ T Cell(Th2)', 'CD4+ T Cell(Th17)','TReg','Proliferating TReg'))
Idents(Immune_subset) <-  Immune_subset$Type
DE_VP16_cd4_T_Cell <- FindMarkers(Immune_subset, ident.1 = "VP16PPARd", ident.2 = "Control", test.use = "MAST", logfc.threshold = .001, min.pct = 0, assay = 'RNA')
write.csv(DE_VP16_cd4_T_Cell,'DE_VP16_cd4_T_Cell.csv')

Idents(Immune_obj) = Immune_obj$Cell_Type
Immune_subset <- subset(Immune_obj,  idents =c('M2 Macrophage'))
Idents(Immune_subset) = Immune_subset$Type
DE_VP16_M2 <- FindMarkers(Immune_subset, ident.1 = "VP16PPARd", ident.2 = "Control", test.use = "MAST", logfc.threshold = .001, min.pct = 0, assay = 'RNA')
write.csv(DE_VP16_M2,'DE_VP16_M2.csv')

DE_VP16_all_T_Cell <- read.csv('DE_VP16_all_T_Cell.csv')
rownames(DE_VP16_all_T_Cell) <-  DE_VP16_all_T_Cell$X
EnhancedVolcano(DE_VP16_all_T_Cell, lab = rownames(DE_VP16_all_T_Cell), x = 'avg_logFC', y = 'p_val_adj', title = 'VP16 vs Control', pCutoff = 10e-15, FCcutoff = 0.25, xlim = c(-2.5, 2.5), ylim = c(0,320), subtitle = 'All Cells' )

DE_VP16_EX8_T_Cell <- read.csv('DE_VP16_EX8_T_Cell.csv') 
rownames(DE_VP16_EX8_T_Cell) <-  DE_VP16_EX8_T_Cell$X
EnhancedVolcano(DE_VP16_EX8_T_Cell, lab = rownames(DE_VP16_EX8_T_Cell), x = 'avg_logFC', y = 'p_val', title = 'VP16 vs Control', pCutoff = 10e-5, FCcutoff = 0.25, xlim = c(-1, 1), ylim = c(0,20), subtitle = 'Ex Cd8+ T Cells' )


EnhancedVolcano(DE_VP16_mem_T_Cell, lab = rownames(DE_VP16_mem_T_Cell), x = 'avg_logFC', y = 'p_val_adj', title = 'VP16 vs Control', pCutoff = 10e-15, FCcutoff = 0.25, xlim = c(-2.5, 2.5), ylim = c(0,320), subtitle = 'All Cells' )
EnhancedVolcano(DE_VP16_M2, lab = rownames(DE_VP16_M2), x = 'avg_logFC', y = 'p_val_adj', title = 'VP16 vs Control', pCutoff = 10e-15, FCcutoff = 0.25, xlim = c(-2.5, 2.5), ylim = c(0,320), subtitle = 'All Cells' )

#prep vp16 only set
Idents(Immune_obj) <- Immune_obj$Type
Immune_VP16 <- subset(Immune_obj,  idents =c('Control', 'VP16PPARd'))
my_levels <- c('Control' ,'VP16PPARd')
Immune_VP16$Type <- factor(x = Immune_VP16$Type, levels = my_levels)

DefaultAssay(Immune_VP16) <- "RNA"
Immune_VP16 <- NormalizeData(object = Immune_VP16, normalization.method = "LogNormalize", assay = "RNA")
Immune_VP16 <- magic(Immune_VP16)
DimPlot(Immune_VP16, group.by = 'Cell_Type')

levels(Immune_VP16$Type)
Gene.list <- c('Ifi208',	'Hspa1a',	'Hspa1b',	'Jun',	'Klf2',	'Gzmk',	'H2-D1',	'Ifit3',	'Ifi209',	'Btg2',	'Ifit1bl1',	'Bst2',	'Ifit1',	'Fos',	'Isg15',	'Cxcr3',	'Ccl9',	'Irf7',	'Apobec1',	'Slfn1',	'Rps27',	'Bcl2',	'Ms4a4c',	'Ier5',	'Ly6c2',	'Itgb7',	'Ifi206')
Gene.list.filtered <- c('Tcf7','Cx3cr1','Cxcr1','Pdcd1','Havcr2','Tox')
Gene.list.Ifn <- c('Ifi208','Ifit3','Ifi209','Ifit1bl1','Ifit1','Isg15','Irf7','Itgb7','Ifi206')

Gene.list.filtered <- c('Cxcr3','Cxcl9','Cxcl10','Cxcl11','Gzmk','Pdcd1','Havcr2','Tox')
All <- VlnPlot(Immune_VP16, group.by = 'Cell_Type' ,ncol = 2,  split.by = "Type", features = Gene.list.filtered, pt.size = 0, assay = "MAGIC_RNA",  cols = c('#1b9e77' ,'#d95f02'), log = FALSE, split.plot = TRUE, combine = TRUE, stack=T, flip=T) + 
  theme(legend.position = 'none') + 
  geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd=.25) + xlab('') + 
  scale_y_continuous(expand = c(0.2,0), breaks = scales::breaks_extended(n = 3)) + 
  theme(text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))
All$layers[[1]]$aes_params$size = .15
All
ggsave(file = paste0('Immune_Gene_filtered.pdf'), plot=All, width=12, height=6, units="in")

VlnPlot(Immune_VP16, group.by = 'Cell_Type' , ncol = 2, split.by = "Type", features = c('Cxcr6','Cxcl10', 'Cxcl9','Cxcr3'), pt.size = 0, assay = "MAGIC_RNA", cols = c('#1b9e77' ,'#d95f02'), split.plot = TRUE) +  theme(legend.position = 'none', axis.title.x = element_blank()) + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0)
ggsave(file="Highlighted_genes.pdf", plot=plot, width=20, height=15)

DotPlot(Immune_VP16, features = Gene.list, assay = 'MAGIC_RNA', group.by = 'Cell_Type') + labs(y= "Cell type", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 5)) +
  theme(text = element_text(size=5), axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust= .01), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())

FeaturePlot(object = Immune_VP16, features = c('Cxcl9','Cxcl10','Ifna','Gzmb','Prf1'), pt.size = .001) + scale_color_viridis(option = 'B') +
  theme(plot.title = element_blank(), text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())


#do differential expression for each cell type

Idents(Immune_obj) <- Immune_obj$Type
DE_VP16_Control <- FindMarkers(Immune_obj, ident.1 = "VP16PPARd", ident.2 = "Control", test.use = "MAST", logfc.threshold = .01, min.pct = 0, assay = 'RNA')
EnhancedVolcano(DE_VP16_Control, lab = rownames(DE_VP16_Control), x = 'avg_logFC', y = 'p_val_adj', title = 'VP16 vs Control', pCutoff = 10e-15, FCcutoff = 0.25, xlim = c(-1.5,1.5), ylim = c(0,320), subtitle = 'All Cells' )

Cell_Types <- levels(Immune_obj$Cell_Type)
Gene.list.Ifn <- c('Ifi208','Ifit3','Ifi209','Ifit1bl1','Ifit1','Isg15','Irf7','Itgb7','Ifi206')

DE_table_all = Idents(Immune_VP16) %>% levels() %>% intersect(Cell_Types) %>% lapply(get_lfc_celltype, seurat_obj = Immune_VP16, condition_colname = "Type", condition_oi = 'VP16PPARd', condition_reference = 'Control', celltype_col = NULL, expression_pct = 0.10) %>% reduce(full_join)
DE_table_all[is.na(DE_table_all)] = 0

DE_table_all_frame <- data.frame(DE_table_all)
write.csv(DE_table_all_frame,'DE_VP16_Celltype.csv')
rownames(DE_table_all_frame) <- DE_table_all_frame$gene
DE_table_all_filtered <- DE_table_all_frame[Gene.list.Ifn,]

# make LFC heatmap
lfc_matrix = DE_table_all_filtered  %>% select(-gene) %>% as.matrix() %>% magrittr::set_rownames(DE_table_all_filtered$gene)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

vis_lfc = lfc_matrix

colnames(vis_lfc) = vis_lfc %>% colnames() %>% make.names()
plot_lfc = vis_lfc %>% make_threecolor_heatmap_ggplot("Genes","Cell Type", low_color = "midnightblue",mid_color = "white", mid = median(vis_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"), legend.text = element_text(size = 9)) + coord_flip()
plot_lfc

#Create a module score for interferon for umaps
Immune_VP16 <- AddModuleScore(object = Immune_VP16, features = list(Gene.list.Ifn), name = 'Interferon_Score', assay = 'MAGIC_RNA')
FeaturePlot(object = Immune_VP16, features = 'Interferon_Score1', pt.size = .001)  + scale_color_viridis(option = 'B') +
  theme( text = element_text(size=6),  axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())


# Check if we see the opposite thing in CPT1aKO

Idents(Immune_obj) <- Immune_obj$Type
Immune_CPT1a <- subset(Immune_obj,  idents =c('Control', 'CPT1aKO'))
my_levels <- c('Control' ,'CPT1aKO')
Immune_CPT1a$Type <- factor(x = Immune_CPT1a$Type, levels = my_levels)

DefaultAssay(Immune_CPT1a) <- "RNA"
Immune_CPT1a <- NormalizeData(object = Immune_CPT1a, normalization.method = "LogNormalize", assay = "RNA")
Immune_CPT1a <- magic(Immune_CPT1a)
DimPlot(Immune_CPT1a, group.by = 'Cell_Type')

#Vln plots
Gene.list <- c("Cxcr3",'Cxcl10','Cxcl9','Tcf7','Cx3cr1','Pdcd1','Havcr2','Tox','Ifit3','Ifi209','Ifit1bl1','Ifit1','Isg15','Irf7','Itgb7','Ifi206')

All <- VlnPlot(Immune_CPT1a, group.by = 'Cell_Type' ,  split.by = "Type", features = Gene.list, pt.size = 0, assay = "MAGIC_RNA",  cols = c('#1b9e77' ,'#d95f02'), log = FALSE, split.plot = TRUE, combine = TRUE, stack=T, flip=T) + 
  theme(legend.position = 'none') + 
  geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd=.25) + xlab('') + 
  scale_y_continuous(expand = c(0.2,0), breaks = scales::breaks_extended(n = 3)) + 
  theme(text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))
All$layers[[1]]$aes_params$size = .15
All
ggsave(file = paste0('Immune_Gene_filtered.pdf'), plot=All, width=5, height=5, units="in")

#do DE for CPT1aKO
Cell_Types <- levels(Immune_obj$Cell_Type)
Idents(Immune_CPT1a) <-  Immune_CPT1a$Cell_Type
DE_table_all = Idents(Immune_CPT1a) %>% levels() %>% intersect(Cell_Types) %>% lapply(get_lfc_celltype, seurat_obj = Immune_CPT1a, condition_colname = "Type", condition_oi = 'CPT1aKO', condition_reference = 'Control', celltype_col = NULL, expression_pct = 0.10) %>% reduce(full_join)
DE_table_all[is.na(DE_table_all)] = 0

DE_table_all_frame <- data.frame(DE_table_all)
write.csv(DE_table_all_frame,'DE_CPT1aKO_Celltype.csv')
rownames(DE_table_all_frame) <- DE_table_all_frame$gene
DE_table_all_filtered <- DE_table_all_frame[Gene.list,]

# make LFC heatmap
lfc_matrix = DE_table_all_filtered  %>% select(-gene) %>% as.matrix() %>% magrittr::set_rownames(DE_table_all_filtered$gene)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

vis_lfc = lfc_matrix

colnames(vis_lfc) = vis_lfc %>% colnames() %>% make.names()
plot_lfc = vis_lfc %>% make_threecolor_heatmap_ggplot("Genes","Cell Type", low_color = "midnightblue",mid_color = "white", mid = median(vis_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"), legend.text = element_text(size = 9)) + coord_flip()
plot_lfc

#differential expression anaylsis
Idents(Immune_obj) <- Immune_obj$Type
DE_CPT1aKO_Control <- FindMarkers(Immune_obj, ident.1 = "CPT1aKO", ident.2 = "Control", test.use = "MAST", logfc.threshold = .00, min.pct = 0, assay = 'RNA')
EnhancedVolcano(DE_CPT1aKO_Control, lab = rownames(DE_CPT1aKO_Control), x = 'avg_logFC', y = 'p_val_adj', title = 'CPT1aKO vs Control', pCutoff = 10e-15, FCcutoff = 0.25, xlim = c(-1.5,1.5), ylim = c(0,320), subtitle = 'All Cells' )
write.csv(DE_CPT1aKO_Control,'DE_CPT1aKO.csv')

DimPlot(Immune_obj, group.by = 'Cell_Type')

Idents(Immune_obj) <-  Immune_obj$Cell_Type
Immune_subset <- subset(Immune_obj,  idents =c('CD8+ T Cell', 'Ex CD8+ T Cell', 'CD4+ T Cell(Th2)', 'CD4+ T Cell(Th17)','TReg','Proliferating TReg', 'NKT Cell','GD T lymphocytes'))
Idents(Immune_subset) <-  Immune_subset$Type
DE_CPT1aKO_all_T_Cell <- FindMarkers(Immune_subset, ident.1 = "CPT1aKO", ident.2 = "Control", test.use = "MAST", logfc.threshold = .001, min.pct = 0.0, assay = 'RNA')
write.csv(DE_CPT1aKO_all_T_Cell,'DE_CPT1aKO_all_T_Cell.csv')

Idents(Immune_obj) <-  Immune_obj$Cell_Type
Immune_subset <- subset(Immune_obj,  idents =c('CD8+ T Cell', 'Ex CD8+ T Cell'))
Idents(Immune_subset) <-  Immune_subset$Type
DE_CPT1aKO_all_T_Cell <- FindMarkers(Immune_subset, ident.1 = "CPT1aKO", ident.2 = "Control", test.use = "MAST", logfc.threshold = .001, min.pct = 0.0, assay = 'RNA')
write.csv(DE_CPT1aKO_all_8T_Cell,'DE_CPT1aKO_all_8T_Cell.csv')

Idents(Immune_obj) <-  Immune_obj$Cell_Type
Immune_subset <- subset(Immune_obj,  idents =c('Ex CD8+ T Cell'))
Idents(Immune_subset) <-  Immune_subset$Type
DE_CPT1aKO_EX8_T_Cell <- FindMarkers(Immune_subset, ident.1 = "CPT1aKO", ident.2 = "Control", test.use = "MAST", logfc.threshold = .001, min.pct = 0, assay = 'RNA')
write.csv(DE_CPT1aKO_EX8_T_Cell,'DE_CPT1aKO_EX8_T_Cell.csv')

Idents(Immune_obj) <-  Immune_obj$Cell_Type
Immune_subset <- subset(Immune_obj,  idents =c('CD8+ T Cell'))
Idents(Immune_subset) <-  Immune_subset$Type
DE_CPT1aKO_Cd8_T_Cell <- FindMarkers(Immune_subset, ident.1 = "CPT1aKO", ident.2 = "Control", test.use = "MAST", logfc.threshold = .001, min.pct = 0, assay = 'RNA')
write.csv(DE_CPT1aKO_Cd8_T_Cell,'DE_CPT1aKO_Cd8_T_Cell.csv')

Idents(Immune_obj) <-  Immune_obj$Cell_Type
Immune_subset <- subset(Immune_obj,  idents =c('CD4+ T Cell(Th2)', 'CD4+ T Cell(Th17)','TReg','Proliferating TReg'))
Idents(Immune_subset) <-  Immune_subset$Type
DE_CPT1aKO_cd4_T_Cell <- FindMarkers(Immune_subset, ident.1 = "CPT1aKO", ident.2 = "Control", test.use = "MAST", logfc.threshold = .001, min.pct = 0, assay = 'RNA')
write.csv(DE_CPT1aKO_cd4_T_Cell,'DE_CPT1aKO_cd4_T_Cell.csv')

Idents(Immune_obj) <-  Immune_obj$Cell_Type
Immune_subset <- subset(Immune_obj,  idents =c('M2 Macrophage'))
Idents(Immune_subset) <-  Immune_subset$Type
DE_CPT1aKO_M2_Cell <- FindMarkers(Immune_subset, ident.1 = "CPT1aKO", ident.2 = "Control", test.use = "MAST", logfc.threshold = .001, min.pct = 0, assay = 'RNA')
write.csv(DE_CPT1aKO_M2_Cell,'DE_CPT1aKO_M2_Cell.csv')

Idents(Immune_obj) <-  Immune_obj$Cell_Type
Immune_subset <- subset(Immune_obj,  idents =c('M1 Macrophage'))
Idents(Immune_subset) <-  Immune_subset$Type
DE_CPT1aKO_M1_Cell <- FindMarkers(Immune_subset, ident.1 = "CPT1aKO", ident.2 = "Control", test.use = "MAST", logfc.threshold = .001, min.pct = 0, assay = 'RNA')
write.csv(DE_CPT1aKO_M1_Cell,'DE_CPT1aKO_M1_Cell.csv')

#make fancy smancy plots for presentation

DimPlot(Immune_obj, reduction = "umap", label = FALSE, pt.size = 1)
ggsave(file = 'Immune_Umap.png', width=7.5, height=6, units="in")

Prop_table<- prop.table(x = table(Immune_VP16$Cell_Type, Immune_VP16$Type), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table1 <- Prop_Table
Prop_Table1$Var1 <- factor(Prop_Table1$Var1)
Prop_Table1$Var2 <- factor(Prop_Table1$Var2,levels = c('Control','VP16PPARd'))
plot <- ggplot(data = Prop_Table1, aes(Var1, Freq, fill=Var2)) + geom_bar(position="dodge", stat="identity", na.rm = TRUE) +  scale_fill_manual(values = c('#1b9e77' ,'#d95f02')) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot
ggsave(file="Immune_vp16_prop.png", plot=plot, width=10, height=8)

Idents(Immune_VP16) <- Immune_VP16$Cell_Type
Immune_VP16_Tcell <- subset(Immune_VP16,  idents = c('CD8+ T Cell', 'Ex CD8+ T Cell','GD T lymphocytes','CD4+ T Cell(Th17)', 'CD4+ T Cell(Th2)', 'TReg', 'Proliferating TReg', 'NKT Cell'))
my_levels <- c('CD8+ T Cell', 'Ex CD8+ T Cell','GD T lymphocytes','CD4+ T Cell(Th17)', 'CD4+ T Cell(Th2)', 'TReg', 'Proliferating TReg', 'NKT Cell')
Immune_VP16_Tcell$Cell_Type <- factor(x = Immune_VP16_Tcell$Cell_Type, levels = my_levels)

Prop_table<- prop.table(x = table(Immune_VP16_Tcell$Cell_Type, Immune_VP16_Tcell$Type), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table1 <- Prop_Table

my_levels <- c('CD8+ T Cell', 'Ex CD8+ T Cell','GD T lymphocytes','CD4+ T Cell(Th17)', 'CD4+ T Cell(Th2)', 'TReg', 'Proliferating TReg', 'NKT Cell')
Prop_Table1$Var1 <- factor(Prop_Table1$Var1,levels = my_levels)
Prop_Table1$Var2 <- factor(Prop_Table1$Var2,levels = c('Control','VP16PPARd'))
plot <- ggplot(data = Prop_Table1, aes(Var1, Freq, fill=Var2)) + geom_bar(position="dodge", stat="identity", na.rm = TRUE) +  scale_fill_manual(values = c('#1b9e77' ,'#d95f02')) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot
ggsave(file="Immune_vp16_prop.png", plot=plot, width=10, height=8)

Heatmap_data <- data.frame(DE_table_all_filtered)
library(reshape)
Heatmap_data <-melt(Heatmap_data)
pathwayColorsDiff = rev(brewer.pal(9, "RdBu"))
ggplot(Heatmap_data, aes(gene, variable, fill= value)) + geom_tile() + scale_fill_gradientn(name = "VP16PPARd vs Control difference", colors = pathwayColorsDiff)

Gene.list.filtered <- c('Cxcr3','Cxcl9','Cxcl10','Cxcl11')
All <- VlnPlot(Immune_CPT1a, group.by = 'Cell_Type' ,ncol = 1,  split.by = "Type", features = Gene.list.filtered, pt.size = 0, assay = "MAGIC_RNA",  cols = c('#1b9e77' ,'#d95f02'), log = FALSE, split.plot = TRUE, combine = TRUE, stack=T, flip=T) + 
  theme(legend.position = 'none') + 
  geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd=.25) + xlab('') + 
  scale_y_continuous(expand = c(0.2,0), breaks = scales::breaks_extended(n = 3)) + 
  theme(text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))
All$layers[[1]]$aes_params$size = .15
All
ggsave(file = paste0('Immune_Gene_filtered.png'), plot=All, width=3, height=3, units="in")

De_CPT1aKO <- read.csv('De_CPT1aKO.csv')

EnhancedVolcano(De_CPT1aKO, lab = De_CPT1aKO$X, x = 'avg_log2FC', y = 'p_val_adj', title = 'CPT1aKO vs Control', pCutoff = 10e-25, FCcutoff = 0.25, xlim = c(-1.5,1.5), ylim = c(0,175), subtitle = 'All Cells', selectLab = c('Gzmb','Gzmc','Ifng','Ctsb','Anxa2','Fos','Ccl7','Cd8b1','Gzma', 'Nkg7','Ctla','Slpi'), raster = TRUE, drawConnectors = TRUE )
ggsave(file = paste0('CPT1a_Volcano.png'), width=6, height=6, units="in")

Idents(Immune_CPT1a) <- Immune_CPT1a$Cell_Type
Immune_VP16_Tcell <- subset(Immune_CPT1a,  idents = c('CD8+ T Cell', 'Ex CD8+ T Cell','GD T lymphocytes','CD4+ T Cell(Th17)', 'CD4+ T Cell(Th2)', 'TReg', 'Proliferating TReg', 'NKT Cell'))
my_levels <- c('CD8+ T Cell', 'Ex CD8+ T Cell','GD T lymphocytes','CD4+ T Cell(Th17)', 'CD4+ T Cell(Th2)', 'TReg', 'Proliferating TReg', 'NKT Cell')
Immune_VP16_Tcell$Cell_Type <- factor(x = Immune_VP16_Tcell$Cell_Type, levels = my_levels)

Prop_table<- prop.table(x = table(Immune_VP16_Tcell$Cell_Type, Immune_VP16_Tcell$Type), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table1 <- Prop_Table

my_levels <- c('CD8+ T Cell', 'Ex CD8+ T Cell','GD T lymphocytes','CD4+ T Cell(Th17)', 'CD4+ T Cell(Th2)', 'TReg', 'Proliferating TReg', 'NKT Cell')
Prop_Table1$Var1 <- factor(Prop_Table1$Var1,levels = my_levels)
Prop_Table1$Var2 <- factor(Prop_Table1$Var2,levels = c('Control','CPT1aKO'))
plot <- ggplot(data = Prop_Table1, aes(Var1, Freq, fill=Var2)) + geom_bar(position="dodge", stat="identity", na.rm = TRUE) +  scale_fill_manual(values = c('#1b9e77' ,'#d95f02')) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot
ggsave(file="Immune_Cpt1a_prop.png", plot=plot, width=10, height=8)

Prop_table<- prop.table(x = table(Immune_obj$Cell_Type, Immune_obj$Type), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table1 <- Prop_Table
Prop_Table1$Var1 <- factor(Prop_Table1$Var1)
Prop_Table1$Var2 <- factor(Prop_Table1$Var2,levels = c('Control','VP16PPARd','CPT1aKO'))
plot <- ggplot(data = Prop_Table1, aes(Var1, Freq, fill=Var2)) + geom_bar(position="dodge", stat="identity", na.rm = TRUE) +  scale_fill_manual(values = c('#1b9e77' ,'#d95f02',"#145FB4")) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot
ggsave(file="Immune_Cpt1a_prop.png", plot=plot, width=10, height=8)

DimPlot(Immune_obj, group.by = 'Cell_Type', label = FALSE, repel = TRUE, label.size = 2)+ theme_cem +
  theme(plot.title = element_blank(),legend.key.size = unit(.1, "cm"), legend.text=element_text(size=3), legend.title =element_text(size=2),  text = element_text(size=5), axis.text.x = element_text(size =  6), axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7), title = )#+ geom_shadowtext(data = clusterMedian, aes(x = UMAP_1,y= UMAP_2, label = CellType, group = as.factor(CellType)), size = .87, bg.colour="black")

ggsave(file="Immune_UMAP_labels.pdf",  width=4.5, height=3.5, units="in")



Idents(Immune_obj) <- Immune_obj$Cell_Type
DE_prolif_treg <- FindMarkers(Immune_obj, ident.1 = "Proliferating TReg",  test.use = "MAST", logfc.threshold = .15, min.pct = .15, assay = 'RNA')
EnhancedVolcano(DE_prolif_treg, lab = rownames(DE_VP16_Control), x = 'avg_logFC', y = 'p_val_adj', title = 'VP16 vs Control', pCutoff = 10e-15, FCcutoff = 0.25, xlim = c(-1.5,1.5), ylim = c(0,320), subtitle = 'All Cells' )
write.csv(DE_prolif_treg,'DE_prolif_treg.csv')


Immune_obj <- readRDS("~/data/tumor_infiltrating_lymphocytes/Immune_obj_Sobj.rds")
Immune_VP16 <- readRDS("~/data/tumor_infiltrating_lymphocytes/Immune_VP16_Sobj.rds")

