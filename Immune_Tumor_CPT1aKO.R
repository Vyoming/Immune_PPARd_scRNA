library(hdf5r)
library(Seurat)
library(scales)
library(plyr)
library(dplyr)
library(future)
library(tidyverse)
library(magrittr)
library(fgsea)
library(viridis)
library(biomaRt)
library(mgcv)
library(shadowtext)
library(ggrepel)

source("~/analysis/scRNA_Intestine/Single_Cell_Functions/cem_helper.R")
source("~/analysis/scRNA_Intestine/Single_Cell_Functions/cem_polyOutline.R")
source("~/analysis/scRNA_Intestine/Single_Cell_Functions/cem_scRNAseq_vyom.R")

getPathways = function(genesOrth = NULL)
{
  pathwaysM = c(gmtPathways("~/analysis/scRNA_Intestine/Beyaz_AA_final.gmt"))
  return(pathwaysM)
}

# parameters
minGenes = 500
minCells = 100

geneColors = brewer.pal(9, "YlOrRd")[-c(1)] ## brewer.pal(9, "YlOrRd") ## rev(viridis::magma(10))
geneColorsDiff = brewer.pal(9, "PiYG")
dimredPlotWidth = 13
dimredPlotHeight = 10


colorsType = c(
  Control = "#1b9e77",
  CPT1aKO = "#d95f02",
  AA = "#d95f02",
  Pge2 = "#d95f02"
)

swatch(colorsType)


# Input files
curFilename = "CPT1aKO"
curTitle = "CPT1aKO vs Control"

seuratAll = Immune_CPT1a


seuratObj = seuratAll
dimred = data.frame(seuratObj@reductions$umap@cell.embeddings)
dimred$Cell = rownames(dimred)
all.equal(names(seuratObj$seurat_clusters), dimred$Cell)
all.equal(names(seuratObj$Type), dimred$Cell)

dimred$Cluster = seuratObj$Cell_Type # Overwrite clusters with cell types for current plots
dimred$Sample = paste0(splitGet(dimred$Cell, "_", 1), "_", splitGet(dimred$Cell, "_", 2))
dimred$Type = seuratObj$Type
dimred$CellType = dimred$Cluster




## Plot clusters on UMAP embedding
{
  library(dplyr)
  
  source("~/analysis/scRNA_Intestine/Single_Cell_Functions/cem_polyOutline.R")
  clusterMedian = dimred %>% group_by(Cluster) %>% summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
  
  hulls = concaveHull(dimred, "UMAP_1", "UMAP_2", "Cluster", alpha = 0.5, extend = T, minPoint = 20)
  
  numClust = length(unique(dimred$Cluster))
  numBatches = length(unique(dimred$Batch))
  numSamples = length(unique(dimred$Sample))
  numTypes = length(unique(dimred$Type))
  
  ggplot(dimred, aes(x = UMAP_1, y = UMAP_2, fill = as.factor(Cluster), color = as.factor(Cluster))) + geom_point(size=.1) + theme_cem + 
    geom_polygon(data = hulls, aes(x = x, y = y, fill = as.factor(Cluster), col = as.factor(Cluster)), alpha = 0.3) + 
    geom_text_repel(data = clusterMedian, aes(label = Cluster, group = as.factor(Cluster)), size = 2, box.padding = .1, force = 5, color = "Black") + 
    scale_fill_manual(values = iwanthue(numClust)) + scale_color_manual(values = iwanthue(numClust)) +
    theme(legend.position = "none", text = element_text(size=5), axis.text.x = element_text(size =  6), axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7))#+ geom_shadowtext(data = clusterMedian, aes(x = UMAP_1,y= UMAP_2, label = CellType, group = as.factor(CellType)), size = .87, bg.colour="black")
  
  ggsave(file= 'CPT1aKO_Hull.pdf', width=2.5, height=2.5, units="in")
  
  
  ggplot(dimred, aes(x = UMAP_1, y = UMAP_2, fill = as.factor(Cluster), color = as.factor(Cluster))) + geom_point() + theme_cem + geom_polygon(data = hulls, aes(x = x, y = y, fill = as.factor(Cluster), col = as.factor(Cluster)), alpha = 0.3) + geom_text(data = clusterMedian, aes(label = Cluster, group = as.factor(Cluster)), size = 5, color = "Black") + theme(legend.position = "none")  + scale_fill_manual(values = iwanthue(numClust)) + scale_color_manual(values = iwanthue(numClust)) + ggtitle(paste0(curTitle, ", clusters")) + facet_wrap(~Cluster)
  ggsave(paste0(curFilename, "_dimred_cluster_hull_facet.png"), width = dimredPlotWidth, height = dimredPlotHeight)
  ggsave(paste0(curFilename, "_dimred_cluster_hull_facet.pdf"), width = dimredPlotWidth, height = dimredPlotHeight)
  
  ggplot(dimred, aes(x = UMAP_1, y = UMAP_2, fill = as.factor(Cluster), color = as.factor(Cluster))) + geom_point(shape = 21, col = "black") + theme_cem + geom_text(data = clusterMedian, aes(label = Cluster, group = as.factor(Cluster)), size = 5, color = "Black") + theme(legend.position = "none")  + scale_fill_manual(values = iwanthue(numClust)) + scale_color_manual(values = iwanthue(numClust)) + ggtitle(paste0(curTitle, ", clusters"))
  ggsave(paste0(curFilename, "_dimred_cluster.png"), width = dimredPlotWidth, height = dimredPlotHeight)
  ggsave(paste0(curFilename, "_dimred_cluster.pdf"), width = dimredPlotWidth, height = dimredPlotHeight)
  
  ggplot(dimred, aes(x = UMAP_1, y = UMAP_2, fill = as.factor(Type), color = as.factor(Type))) + geom_point() + theme_cem  + scale_fill_manual(values = colorsType) + scale_color_manual(values = colorsType) + ggtitle(paste0(curTitle, ", by Type"))
  ggsave(paste0(curFilename, "_dimred_type.png"), width = dimredPlotWidth, height = dimredPlotHeight)
  ggsave(paste0(curFilename, "_dimred_type.pdf"), width = dimredPlotWidth, height = dimredPlotHeight)
  
  nc = if(numTypes <= 3) numTypes else ceiling(sqrt(numTypes))
  nr = ceiling(numTypes/nc)
  ggplot(dimred, aes(x = UMAP_1, y = UMAP_2, fill = as.factor(Type), color = as.factor(Type))) + geom_point() + theme_cem  + scale_fill_manual(values = colorsType) + scale_color_manual(values = colorsType) + ggtitle(paste0(curTitle, ", by Type")) + facet_wrap(~Type, ncol=nc)
  ggsave(paste0(curFilename, "_dimred_type_facet.png"), width = dimredPlotWidth*nc, height = dimredPlotHeight*nr) 
  ggsave(paste0(curFilename, "_dimred_type_facet.pdf"), width = dimredPlotWidth*nc, height = dimredPlotHeight*nr) 
  
  ggplot(dimred, aes(x = UMAP_1, y = UMAP_2, fill = as.factor(Sample), color = as.factor(Sample))) + geom_point() + theme_cem  + scale_fill_manual(values = iwanthue(numSamples)) + scale_color_manual(values = iwanthue(numSamples)) + ggtitle(paste0(curTitle, ", by Sample"))
  ggsave(paste0(curFilename, "_dimred_sample.png"), width = dimredPlotWidth, height = dimredPlotHeight)
  ggsave(paste0(curFilename, "_dimred_sample.pdf"), width = dimredPlotWidth, height = dimredPlotHeight)
  
  nc = if(numSamples <= 3) numSamples else ceiling(sqrt(numSamples))
  nr = ceiling(numSamples/nc)
  ggplot(dimred, aes(x = UMAP_1, y = UMAP_2, fill = as.factor(Sample), color = as.factor(Sample))) + geom_point() + theme_cem  + scale_fill_manual(values = iwanthue(numSamples)) + scale_color_manual(values = iwanthue(numSamples)) + ggtitle(paste0(curTitle, ", by Sample")) + facet_wrap(Sample~., ncol=nc)
  ggsave(paste0(curFilename, "_dimred_sample_facet.png"), width = dimredPlotWidth*nc, height = dimredPlotHeight*nr) 
  ggsave(paste0(curFilename, "_dimred_sample_facet.pdf"), width = dimredPlotWidth*nc, height = dimredPlotHeight*nr) 
  
  ggplot(dimred, aes(x = UMAP_1, y = UMAP_2, fill = as.factor(Batch), color = as.factor(Batch))) + geom_point() + theme_cem  + scale_fill_manual(values = iwanthue(numBatches)) + scale_color_manual(values = iwanthue(numBatches)) + ggtitle(paste0(curTitle, ", by Batch"))
  ggsave(paste0(curFilename, "_dimred_batch.png"), width = dimredPlotWidth, height = dimredPlotHeight)
  ggsave(paste0(curFilename, "_dimred_batch.pdf"), width = dimredPlotWidth, height = dimredPlotHeight)
  
  nc = if(numBatches <= 3) numBatches else ceiling(sqrt(numBatches))
  nr = ceiling(numBatches/nc)
  ggplot(dimred, aes(x = UMAP_1, y = UMAP_2, fill = as.factor(Batch), color = as.factor(Batch))) + geom_point() + theme_cem  + scale_fill_manual(values = iwanthue(numBatches)) + scale_color_manual(values = iwanthue(numBatches)) + ggtitle(paste0(curTitle, ", by Batch")) + facet_wrap(~Batch, ncol=nc)
  ggsave(paste0(curFilename, "_dimred_batch_facet.png"), width = dimredPlotWidth*nc, height = dimredPlotHeight*nr) 
  ggsave(paste0(curFilename, "_dimred_batch_facet.pdf"), width = dimredPlotWidth*nc, height = dimredPlotHeight*nr) 
  
  
}


## Cluster distribution
{
  clustCounts = dimred %>% group_by(Cluster, Type)  %>% summarise(Count=n())
  clustCounts = dcast(clustCounts, Cluster~Type, value.var = "Count")
  clustCounts$Cluster = factor(clustCounts$Cluster, levels=sort(unique(as.character(clustCounts$Cluster))))
  clustCounts[ is.na(clustCounts) ] = 0
  
  clustCountsNorm = clustCounts
  clustCountsNorm[2:ncol(clustCountsNorm)] = clustCountsNorm[2:ncol(clustCountsNorm)] / rep.row(colSums(clustCountsNorm[2:ncol(clustCountsNorm)]), nrow(clustCountsNorm))
  
  clustCountsNorm2 = melt(clustCountsNorm, by = "Cluster")
  clustCountsNorm2$Type = clustCountsNorm2$variable
  
  
  #clustCountsNorm2$Cluster = #factor(clustCountsNorm2$Cluster, levels = seq(0, max(as.numeric(clustCountsNorm2$Cluster))))
  ggplot(clustCountsNorm2, aes(x = Type, y = value, fill = Type)) + geom_boxplot() + geom_point(position = position_dodge(width = 0.75)) + theme_cem + scale_fill_manual(values = colorsType) + facet_wrap(~Cluster, scales = "free") + expand_limits(y = c(0)) + scale_y_continuous("% of cells in cluster", label=scales::percent) + ggtitle("Cluster assignment by Type")
  ggsave(paste0(curFilename, "_clusterByType_box.png"), width = 13, height = 11)
  ggsave(paste0(curFilename, "_clusterByType_box.pdf"), width = 13, height = 11)
  
  ggplot(clustCountsNorm2, aes(x = Type, y = value, fill = Type)) + geom_bar(stat="identity") + geom_point(position = position_dodge(width = 0.75)) + theme_cem + scale_fill_manual(values = colorsType) + facet_wrap(~Cluster, scales = "free") + expand_limits(y = c(0)) + scale_y_continuous("% of cells in cluster", label=scales::percent) + ggtitle("Cluster assignment by Type")
  ggsave(paste0(curFilename, "_clusterByType_bar2.png"), width = 13, height = 11)
  ggsave(paste0(curFilename, "_clusterByType_bar2.pdf"), width = 13, height = 11)
  
  clustCountsNormSum = summarySE(clustCountsNorm2, measurevar="value", groupvars=c("Cluster", "Type"))
  
  ggplot(clustCountsNormSum, aes(x=Cluster, y=value, fill=Type, group=Type)) + geom_bar(stat="identity", position=position_dodge(0.9)) + theme_cem  + scale_fill_manual(values=colorsType)  + scale_color_manual(values=colorsType) + scale_y_continuous("% of cells in cluster", label=scales::percent) + ggtitle("Cluster assignment by Type")# + geom_point(data=clustCountsNorm2, shape=21, position=position_dodge(0.9)) + geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1, position=position_dodge(0.9)) 
  ggsave(paste0(curFilename, "_clusterByType_bar.png"), width = 13, height = 11)
  ggsave(paste0(curFilename, "_clusterByType_bar.pdf"), width = 13, height = 11)
  
  
  
  clustCountsNormByCluster = clustCounts
  clustCountsNormByCluster[2:ncol(clustCountsNormByCluster)] = clustCountsNormByCluster[2:ncol(clustCountsNormByCluster)] / rep.col(rowSums(clustCountsNormByCluster[2:ncol(clustCountsNormByCluster)]), ncol(clustCountsNormByCluster))
  
  clustCountsNormByCluster2 = melt(clustCountsNormByCluster, by = "Cluster")
  clustCountsNormByCluster2$Type = clustCountsNorm2$variable
  
  clustCountsNormByCluster3 = clustCountsNormByCluster2 %>% group_by(Cluster, Type) %>% summarise(value=sum(value))
  
  ggplot(clustCountsNormByCluster2, aes(x = Type, y = value, fill = Type)) + geom_bar(stat="identity") + theme_cem + scale_fill_manual(values = colorsType) + facet_wrap(~Cluster, scales = "free") + expand_limits(y = c(0)) + scale_y_continuous("% of cells in cluster", label=scales::percent) + ggtitle("Cluster assignment by Type")
  ggsave(paste0(curFilename, "_clusterByTypePerCluster_bar.png"), width = 13, height = 11)
  ggsave(paste0(curFilename, "_clusterByTypePerCluster_bar.pdf"), width = 13, height = 11)
}


# Pathways
{
  pathways = getPathways()
  pathwayNames = gsub("(_UP|_DN|_downreg|_upreg|_down_|_up_)", "\\*", names(pathways))
  pathwaysSelected = pathwayNames # pathways[ pathwayNames %in% c( "SCHUHMACHER_MYC_TARGETS*", "HALLMARK_MYC_TARGETS_V1", "DECP_vs_DECN*", "GCB_GFP_MYC*", "DZgenes_microarray_Nussenzweig", "LZgenes_microarray_Nussenzweig", "CCvsCB*_Teater", "CC_vs_CB_Agirre*", "MEM_vs_CC_Agirre*", "Ezh2mut_vs_WT__CB*", "Ezh2mut_vs_WT__CC*", "PC_vs_GCB_Teater*", "termDiffGenes") ]
  
  pathwayScoresZscore = calculatePathwayScores(seuratObj, dimred, pathways, method="seurat")
  
  
  
  
  
  
  pathwayColors = brewer.pal(9, "RdYlGn")
  pathwayColorsDiff = rev(brewer.pal(9, "RdBu"))
  
  pathwayScores = pathwayScoresZscore
  
  pathwayScoresW = dcast(pathwayScores, Cell+Cluster+Sample+Type+CellType+UMAP_1+UMAP_2~Pathway, value.var="Score")
  colnames(pathwayScoresW)[9:ncol(pathwayScoresW)] = paste0("Score_", colnames(pathwayScoresW)[9:ncol(pathwayScoresW)])
  
  dimredP = pathwayScoresW
  
  cellColors = iwanthue(length(unique(pathwayScoresW$CellType)))
  names(cellColors) = unique(pathwayScoresW$CellType)
  annColors = list(Pseudotime = brewer.pal(11, "Spectral"), Type = colorsType)
  
  source("~/analysis/scRNA_Intestine/Single_Cell_Functions/cem_polyOutline.R")
  clusterMedian = dimredP %>%
    group_by(CellType) %>%
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
  
  hulls2 = concaveHull(dimredP, "UMAP_1", "UMAP_2", "CellType", alpha = 0.5, extend = T, minPoint = 20)
  
  curPath = getwd()
  plotPath = paste0(curPath, "/PathwayPlots_", curFilename)
  dir.create(plotPath, showWarnings = F)
  dir.create(paste0(plotPath, "/umap_pointAll"), showWarnings = F)
  dir.create(paste0(plotPath, "/umap_pointByType"), showWarnings = F)
  dir.create(paste0(plotPath, "/umap_densityAll"), showWarnings = F)
  dir.create(paste0(plotPath, "/umap_densityByType"), showWarnings = F)
  dir.create(paste0(plotPath, "/umap_densityDifference"), showWarnings = F)
  dir.create(paste0(plotPath, "/byCluster_density"), showWarnings = F)
  setwd(plotPath)
  
  for(curPathway in unique(pathwayScoresZscore$Pathway)){	
    curPathwayScore = pathwayScoresZscore[pathwayScoresZscore$Pathway %in% curPathway, ]
    if(nrow(curPathwayScore) == 0 | sum(is.nan(curPathwayScore$Score) > 0)) next
    
    rangeL = quantile(curPathwayScore$Score, 0.01, na.rm = T)
    rangeH = quantile(curPathwayScore$Score, 0.95, na.rm = T)
    if(rangeL == rangeH) { rangeL = min(curPathwayScore$Score); rangeH = max(curPathwayScore$Score); }
    
    #ggplot(curPathwayScore, aes(x=UMAP_1, y=UMAP_2, fill=Score, color=Score)) + geom_point(alpha=0.7) + theme_cem + scale_color_gradientn(colors=pathwayColors, limits = c(rangeL, rangeH), oob = squish ) + scale_fill_gradientn(colors=pathwayColors, limits = c(rangeL, rangeH), oob = squish ) + ggtitle(paste0(curPathway, "\nModule score, z-score"))
    #ggsave(paste0("umap_pointAll/", curFilename, "pathwayScoreZscore_umap_pointAll_", curPathway, ".png"), width=dimredPlotWidth, height=dimredPlotHeight, dpi=150)
    #ggsave(paste0("umap_pointAll/", curFilename, "pathwayScoreZscore_umap_pointAll_", curPathway, ".pdf"), width=dimredPlotWidth, height=dimredPlotHeight)
    
    #ggplot(curPathwayScore, aes(x=UMAP_1, y=UMAP_2, fill=Score, color=Score)) + geom_point(alpha=0.7) + theme_cem + scale_color_gradientn(colors=pathwayColors, limits = c(rangeL, rangeH), oob = squish ) + scale_fill_gradientn(colors=pathwayColors, limits = c(rangeL, rangeH), oob = squish ) + ggtitle(paste0(curPathway, "\nModule score, z-score")) + facet_grid(~Type)
    #ggsave(paste0("umap_pointByType/", curFilename, "pathwayScoreZscore_umap_pointByType_", curPathway, ".png"), width=dimredPlotWidth*3-1, height=dimredPlotHeight, dpi=150)
    #ggsave(paste0("umap_pointByType/", curFilename, "pathwayScoreZscore_umap_pointByType_", curPathway, ".pdf"), width=dimredPlotWidth*3-1, height=dimredPlotHeight)
    
    
    sm = smoothScore2d(curPathwayScore$Score, curPathwayScore$UMAP_1, curPathwayScore$UMAP_2, numGrid=100, knn=50, m=2)
    #ggplot(sm) + geom_tile(aes(x = x, y = y, fill = score)) + theme_cem + scale_fill_gradientn(name = "Score", colors = pathwayColors) + ggtitle(paste0(curPathway, "\nModule score, z-score")) + geom_polygon(data = hulls2, aes(x = x, y = y, group = CellType), alpha = 0.3, fill=NA, color="#666666") + geom_text(data = clusterMedian, aes(x=UMAP_1, y=UMAP_2, label = CellType, group = as.factor(CellType)), size = 5, color = "Black")
    #ggsave(paste0("umap_densityAll/", curFilename, "pathwayScoreZscore_umap_densityAll_", curPathway, ".png"), width = dimredPlotWidth, height = dimredPlotHeight, dpi=150)
    #ggsave(paste0("umap_densityAll/", curFilename, "pathwayScoreZscore_umap_densityAll_", curPathway, ".pdf"), width = dimredPlotWidth, height = dimredPlotHeight)
    
    
    sm2 = ddply(curPathwayScore, "Type", function(x) { smoothScore2d( x$Score, x$UMAP_1, x$UMAP_2, numGrid=100, knn=50, m=2, xrng=range(curPathwayScore$UMAP_1), yrng=range(curPathwayScore$UMAP_2)) } )	
    #ggplot(sm2) + geom_tile(aes(x = x, y = y, fill = score)) + theme_cem + scale_fill_gradientn(name = "Score", colors = pathwayColors) + ggtitle(paste0(curPathway, "\nModule score, z-score")) + geom_polygon(data = hulls2, aes(x = x, y = y, group = CellType), alpha = 0.3, fill=NA, color="#666666") + geom_text(data = clusterMedian, aes(x=UMAP_1, y=UMAP_2, label = CellType, group = as.factor(CellType)), size = 5, color = "Black") + facet_grid(~Type)
    #ggsave(paste0("umap_densityByType/", curFilename, "pathwayScoreZscore_umap_densityByType_", curPathway, ".png"), width = dimredPlotWidth*3-1, height = dimredPlotHeight, dpi=150)
    #ggsave(paste0("umap_densityByType/", curFilename, "pathwayScoreZscore_umap_densityByType_", curPathway, ".pdf"), width = dimredPlotWidth*3-1, height = dimredPlotHeight)	
    
    
    #rangeL = quantile(curPathwayScore$Score, 0.001, na.rm = T)
    #rangeH = quantile(curPathwayScore$Score, 0.999, na.rm = T)
    #ggplot(curPathwayScore, aes_string(x = "Score", fill = "Type", color = "Type")) + geom_density(alpha=0.5) + facet_grid(Cluster~.) + theme_cem + ggtitle(paste0(curPathway, "\nModule score, z-score by cluster"))+ scale_fill_manual(values = colorsType) + scale_color_manual(values = colorsType) + coord_trans(limx=c(rangeL, rangeH))
    #ggsave(paste0("byCluster_density/", curFilename, "pathwayScoreZscore_byCluster_density_", curPathway, ".png"), width = dimredPlotHeight*1.5, height = dimredPlotWidth, dpi=150)
    #ggsave(paste0("byCluster_density/", curFilename, "pathwayScoreZscore_byCluster_density_", curPathway, ".pdf"), width = dimredPlotHeight*1.5, height = dimredPlotWidth)	
    
    ##
    
    sm3 = smoothScore2d(curPathwayScore$Score, curPathwayScore$UMAP_1, curPathwayScore$UMAP_2, type=curPathwayScore$Type, numGrid=100, knn=50, m=2)
    sm3W = dcast(sm3, x+y~type, value.var="score")
    sm3W$Diff = sm3W$CPT1aKO - sm3W$Control
    rangeDiff = max(abs(sm3W$Diff))
    
    ggplot(sm3W) + geom_tile(aes(x = x, y = y, fill = Diff)) + theme_cem + 
      scale_fill_gradientn(name = "CPT1aKO vs Control difference", colors = pathwayColorsDiff, limits = c(-rangeDiff, rangeDiff)) + 
      geom_polygon(data = hulls2, aes(x = x, y = y, group = CellType), alpha = 0.3, fill=NA, color="#666666", size = .1) + 
      geom_text_repel(data = clusterMedian, aes(x = UMAP_1, y = UMAP_2, label = CellType, group = as.factor(CellType)), size = 2, box.padding = .1, force = 5, color = "Black") +
      theme(legend.position = "none", text = element_text(size=5),  axis.text.x = element_text(size =  6), axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7))
    ggsave(file= paste0(curFilename, "_umap_densityDifference_", curPathway, ".pdf"), width=2.5, height=2.5, units="in")
    #ggsave(paste0("umap_densityDifference/", curFilename, "pathwayScoreZscore_umap_densityDifference_", curPathway, "_CPT1aKOVsControl.png"), width = dimredPlotWidth, height = dimredPlotHeight, dpi=150)
    #ggsave(paste0("umap_densityDifference/", curFilename, "pathwayScoreZscore_umap_densityDifference_", curPathway, "_CPT1aKOVsControl.pdf"), width = dimredPlotWidth, height = dimredPlotHeight)		
  }
  setwd(curPath)
  
}

## Misc integration: pathway score & Cell type & pseudotime 
{	
  pathwayScores = pathwayScoresZscore
  
  pathwayScoresW = dcast(pathwayScores, Cell+Cluster+Sample+Type+CellType+Pseudotime+UMAP_1+UMAP_2~Pathway, value.var="Score")
  colnames(pathwayScoresW)[9:ncol(pathwayScoresW)] = paste0("Score_", colnames(pathwayScoresW)[9:ncol(pathwayScoresW)])
  
  dimredP = pathwayScoresW
  
  #cellColors = c("Stem 1"="#377eb8", "Stem 2"="#7570b3", "Transit Amplifying"="#e41a1c", "Enterocyte Progenitor"="#ffd92f","Enterocyte (Proximal)"="#4daf4a", "Enterocyte (Distal)"="#a65628", "Enteroendocrine"="#e78ac3")
  cellColors = iwanthue(length(unique(pathwayScoresW$CellType)))
  names(cellColors) = unique(pathwayScoresW$CellType)
  
  ggplot(pathwayScoresW, aes(x=UMAP_1, y=UMAP_2, color=CellType)) + geom_point(size=1, alpha=0.5) + theme_cem + scale_color_manual(values=cellColors)
  ggsave(paste0(curFilename, "_CellType_all.png"), width=12, height=9)
  ggplot(pathwayScoresW, aes(x=UMAP_1, y=UMAP_2, color=CellType)) + geom_point(size=1, alpha=0.5) + theme_cem + scale_color_manual(values=cellColors) + facet_wrap(~CellType)
  ggsave(paste0(curFilename, "_CellType_facet.png"), width=16, height=12)
  
  annColors = list(Pseudotime = brewer.pal(11, "Spectral"), Type = colorsType)
  ggplot(pathwayScoresW, aes(x=UMAP_1, y=UMAP_2)) + geom_point(aes(fill = Pseudotime, color = Pseudotime)) + theme_cem + scale_fill_gradientn(colors = annColors$Pseudotime) + scale_color_gradientn(colors = annColors$Pseudotime) + ggtitle(paste0(curTitle, ", UMAP by Pseudotime"))
  ggsave(paste0(curFilename, "_Pseudotime.png"), width=12, height=9)
  
  ggplot(pathwayScoresW, aes(x=UMAP_1, y=UMAP_2)) + geom_point(aes(fill = Pseudotime, color = Pseudotime), size=.1) + theme_cem + 
    scale_fill_gradientn(colors = annColors$Pseudotime) + scale_color_gradientn(colors = annColors$Pseudotime) + 
    geom_polygon(data = hulls2, aes(x = x, y = y, group=CellType), fill=NA, color="black", size=0.1, alpha = 0.3) +
    geom_text_repel(data = clusterMedian, aes(label = CellType, group = as.factor(CellType)), size = 2, box.padding = .1, force = 5, color = "Black") + 
    labs(fill = "Pseudotime",  direction = 'y') + 
    theme(text = element_text(size=5), legend.key.size = unit(0.0, "cm"), legend.text=element_text(size=0), legend.title =element_text(size=0), axis.text.x = element_text(size =  6), axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7))#+ geom_shadowtext(data = clusterMedian, aes(x = UMAP_1,y= UMAP_2, label = CellType, group = as.factor(CellType)), size = .87, bg.colour="black")
  ggsave(file= 'CPT1aKO_Pseudotime_umap_rainbow.pdf', width=2.5, height=2.5, units="in")
  ggsave(paste0(curFilename, "_Pseudotime_hull.png"), width=12, height=9)
  
  ggplot(dimredP, aes(x=Pseudotime, fill=Type)) + geom_density(alpha=0.8) + theme_cem + scale_fill_manual(values=unlist(colorsType))
  ggsave(paste0(curFilename, "_CellType_PseudoTimeGenotype.png"), width=12, height=5)
  
  ggplot(dimredP, aes(x=Pseudotime, fill=CellType)) + geom_density(alpha=0.8) + theme_cem + facet_grid(CellType~., scales="free_y") + scale_fill_manual(values=cellColors)
  ggsave(paste0(curFilename, "_CellType_PseudoTime.png"), width=12, height=12)
  
  ggplot(dimredP, aes(x=Pseudotime, fill=CellType)) + geom_density(alpha=0.4) + theme_cem + scale_fill_manual(values=cellColors) + ylim(0, 2)
  ggsave(paste0(curFilename, "_CellType_PseudoTime2.png"), width=12, height=5)
  color123 <-  as.factor(dimred$Cluster)
  
  ggplot(dimredP, aes(x=Pseudotime, fill=Type)) + geom_density(alpha=0.4, bw=0.9) + theme_cem + facet_grid(CellType~., scales="free_y") + scale_fill_manual(values=unlist(colorsType)) 
  ggsave(paste0(curFilename, "_CellType_PseudoTimeGenotypeCellType.png"), width=12, height=12)
  
  cluster.ids <- c('Stem 1','Stem 2', 'Transit Amplifying', 'Enterocyte Progenitor','Tuft')
  dimredP_half1 = dimredP %>% filter(CellType == cluster.ids)
  ggplot(dimredP_half1, aes(x=Pseudotime, fill=Type)) + geom_density(alpha=0.4, bw=0.5, size = .1) + theme_cem + 
    facet_grid(CellType~., scales="free_y") + scale_fill_manual(values=unlist(colorsType)) +
    scale_y_continuous(breaks = scales::breaks_extended(n = 3)) + 
    theme(strip.background = element_blank(), strip.text.y = element_blank(), legend.position = "none", text = element_text(size=5), axis.text.x = element_text(size =  6), axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7))
  ggsave("CPT1aKO_Pseudodensity_half_1.pdf", width = 1.6, height = 2, units="in")	
  
  
  cluster.ids <- c('Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Paneth', 'Goblet')
  dimredP_half2 = dimredP %>% filter(CellType == cluster.ids)
  ggplot(dimredP_half2, aes(x=Pseudotime, fill=Type)) + geom_density(alpha=0.4, bw=0.5, size = .1) + theme_cem + 
    facet_grid(CellType~., scales="free_y") + scale_fill_manual(values=unlist(colorsType)) +
    scale_y_continuous(breaks = scales::breaks_extended(n = 3)) + 
    theme(strip.background = element_blank(), strip.text.y = element_blank(), legend.position = "none", text = element_text(size=5), axis.text.x = element_text(size =  6), axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7))
  ggsave("CPT1aKO_Pseudodensity_half_2.pdf", width = 1.6, height = 2, units="in")	
  
  
  
  dimredP2 = dimredP[, c(1:8, grep("Score_", colnames(dimredP)))]
  dimredP2 = melt(dimredP2, colnames(dimredP2)[1:8], variable.name="Pathway", value.name="Score")
  
  ggplot(dimredP2, aes(x=Pseudotime, y=Score, color=Type)) + geom_point(size=1, alpha=0.1) + geom_smooth() + theme_cem + facet_grid(Pathway~., scales="free_y") + scale_color_manual(values=unlist(colorsType))
  ggsave(paste0(curFilename, "_CellType_PseudoTimeScore.png"), width=12, height=12)
  
  
  
  ggplot(dimredP, aes(x = UMAP_1, y = UMAP_2, fill = as.factor(CellType), color = as.factor(CellType))) + 
    geom_polygon(data = hulls2, size=1, aes(x = x, y = y, fill = as.factor(CellType), col = as.factor(CellType)), alpha = 0.3) + 
    geom_point(size=1, alpha=0.1) + 
    geom_shadowtext(data = clusterMedian, aes(label = CellType, group = as.factor(CellType), color=CellType), size = 7, family = "mono", fontface = "bold", bg.colour="black") + 
    theme_cem +
    theme(legend.position = "none")  + 
    scale_fill_manual(values = unlist(cellColors)) + 
    scale_color_manual(values = unlist(cellColors)) + 
    ggtitle(paste0(curTitle))
  
  ggsave(paste0(curFilename, "_CellType_hull.png"), width = dimredPlotWidth, height = dimredPlotHeight)
  ggsave(paste0(curFilename, "_CellType_hull.pdf"), width = dimredPlotWidth, height = dimredPlotHeight)
  
  #
  
  ggplot(dimredP, aes(x = UMAP_1, y = UMAP_2, fill = as.factor(CellType), color = as.factor(CellType))) + 
    geom_polygon(data = hulls2, size=1, aes(x = x, y = y, fill = as.factor(CellType), col = as.factor(CellType)), alpha = 0.3) + 
    geom_point(size=1, alpha=0.1) + 
    geom_shadowtext(data = clusterMedian, aes(label = CellType, group = as.factor(CellType), color=CellType), size = 7, family = "mono", fontface = "bold", bg.colour="black") + 
    theme_cem +
    theme(legend.position = "none")  + 
    scale_fill_manual(values = unlist(cellColors)) + 
    scale_color_manual(values = unlist(cellColors)) + 
    ggtitle(paste0(curTitle)) +
    facet_wrap(~CellType)
  
  ggsave(paste0(curFilename, "_CellType_hullFacet.png"), width = dimredPlotWidth, height = dimredPlotHeight)
  ggsave(paste0(curFilename, "_CellType_hullFacet.pdf"), width = dimredPlotWidth, height = dimredPlotHeight)
  
  
  ##
  
  
  curCol = "CellType"
  celltype_cluster = as.data.frame(table( dimredP[, curCol], dimredP$Cluster ))
  celltype_sample = as.data.frame(table( dimredP[, curCol], dimredP$Sample ))
  
  celltype_cluster = dcast(celltype_cluster, Var1~Var2, value.var = "Freq")
  celltype_sample = dcast(celltype_sample, Var1~Var2, value.var = "Freq")
  
  celltype_cluster_norm = celltype_cluster
  celltype_sample_norm = celltype_sample
  
  celltype_cluster_norm[2:ncol(celltype_cluster_norm)] = celltype_cluster_norm[2:ncol(celltype_cluster_norm)] / rep.row(colSums(as.matrix(celltype_cluster_norm[2:ncol(celltype_cluster_norm)])), nrow(celltype_cluster_norm))
  celltype_sample_norm[2:ncol(celltype_sample_norm)] = celltype_sample_norm[2:ncol(celltype_sample_norm)] / rep.row(colSums(as.matrix(celltype_sample_norm[2:ncol(celltype_sample_norm)])), nrow(celltype_sample_norm))
  
  
  celltype_sample_norm_filt = celltype_sample_norm[ rowMaxs(as.matrix(celltype_sample_norm[,2:ncol(celltype_sample_norm)])) > 0.005, ]
  
  write.tsv(celltype_cluster, paste0(curFilename, "_CellType_ByCluster.txt"), row.names = F)
  write.tsv(celltype_sample, paste0(curFilename, "_CellType_BySample.txt"), row.names = F)
  
  write.tsv(celltype_cluster_norm, paste0(curFilename, "_CellType_ByCluster_Norm.txt"), row.names = F)
  write.tsv(celltype_sample_norm, paste0(curFilename, "_CellType_BySample_Norm.txt"), row.names = F)
  
  celltype_sample_norm2 = celltype_sample_norm[ rowMeans(celltype_sample_norm[, 2:ncol(celltype_sample_norm)]) >= 0.001, ]
  celltype_sample_norm2 = melt(celltype_sample_norm2, by = "Var1")
  colnames(celltype_sample_norm2) = c("CellType", "Sample", "Proportion")
  celltype_sample_norm2$Type = splitGet(as.character(celltype_sample_norm2$Sample), "_", 1)
  celltype_sample_norm2$Type = factor(celltype_sample_norm2$Type, levels=levels(dimred$Type))
  
  ggplot(celltype_sample_norm2, aes(x = Type, y = Proportion, fill = Type)) + geom_boxplot() + geom_point(position = position_dodge(width = 0.75)) + theme_cem + scale_fill_manual(values = colorsType) + facet_wrap(~CellType, scales = "free") + expand_limits(y = c(0))
  ggsave(paste0(curFilename, "_CellType_Box_Type_BySample.png"), width = 20, height = 11)
  
  ggplot(celltype_sample_norm2, aes(x = CellType, y = Proportion, fill = Type)) + geom_boxplot() + geom_point(position = position_dodge(width = 0.75)) + theme_cem + scale_fill_manual(values = colorsType) + expand_limits(y = c(0)) +  facet_wrap(~CellType, scales="free", ncol=5) + theme(axis.text.x = element_blank())
  ggsave(paste0(curFilename, "_CellType_Box2_Type_BySample.png"), width = 13, height = 6)
  
  celltype_sample_norm2$CellType = factor(celltype_sample_norm2$CellType, levels = rev(names(cellColors)))
  ggplot(celltype_sample_norm2, aes(x = CellType, y = Proportion, fill = Type)) + geom_boxplot() + geom_point(position = position_dodge(width = 0.75)) + theme_cem + scale_fill_manual(values = colorsType) + expand_limits(y = c(0)) + coord_flip()
  ggsave(paste0(curFilename, "_CellType_Box3_Type_BySample.png"), width = 10, height = 7)
  
  celltype_sample_norm3 = summarySE(celltype_sample_norm2, measurevar="Proportion", groupvars=c("Type","CellType"))
  celltype_sample_norm3$CellType = factor(celltype_sample_norm3$CellType, levels = rev(names(cellColors)))
  my_levels <- c('Stem 1', 'Stem 2' ,'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft')
  celltype_sample_norm3$CellType <- factor(celltype_sample_norm3$CellType, levels = my_levels)
  ggplot(celltype_sample_norm3, aes(x = CellType, y = Proportion, fill = Type)) + geom_bar(stat="identity", position="dodge") + 
    geom_errorbar(aes(ymin=Proportion-se, ymax=Proportion+se),position=position_dodge(width = 0.85),width=0.3) + theme_cem + 
    scale_fill_manual(values = colorsType) + expand_limits(y = c(0)) + 
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, size =  6, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black"), axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7)) + 
    xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Sample") + scale_y_continuous(expand = expansion(mult = c(0, .1)))
  ggsave( "CPT1aKO_frac_of_cells.pdf", width=4, height=2.5, units="in")
  
  ggplot(celltype_sample_norm3, aes(x = Type, y = Proportion, fill = CellType)) + geom_bar(stat = "identity") + theme_cem + scale_fill_manual(values = iwanthue(42)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + theme(text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title =element_text(size=0), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7))
  ggsave(paste0(curFilename, "_stacked_bar.pdf"), width=1.5, height=2.5, units="in")
  
  ggplot(celltype_sample_norm2, aes(x = Type, y = Proportion, fill = Type)) + geom_bar(stat = "identity") + geom_point(position = position_dodge(width = 0.75)) + theme_cem + scale_fill_manual(values = colorsType) + facet_wrap(~CellType, scales = "free") + expand_limits(y = c(0))
  ggsave(paste0(curFilename, "_CellType_Bar_Type_BySample2.png"), width = 20, height = 11)
  
  celltype_sample_norm3 = melt(celltype_sample_norm_filt, by = "Var1")
  colnames(celltype_sample_norm3) = c("CellType", "Sample", "Proportion")
  ggplot(celltype_sample_norm3, aes(x = Sample, y = Proportion, fill = CellType)) + geom_bar(stat = "identity") + theme_cem + scale_fill_manual(values = iwanthue(42))
  ggsave(paste0(curFilename, "_CellType_Bar_CellType_BySample2.png"), width = 15, height = 7)
  
}


## Plot gene expression on dimred embedding
{
  assay = Matrix::t(seuratObj@assays$MAGIC_RNA@data)
  #assay = t(seuratObj@assays$RNA@counts)
  
  geneList = unique(c('Ifi208','Cxcr3','Cxcl10','Cxcl9','Gzmc','Gzmb','Xcl1', 'Nkg7', 'Trbc1','Klrd1', 'H2-T23','Cd3g','Ctla2a','Ccl5', 'Gzma','Ctsw', 'Cd3d','Trbc2','Cd3e', 'Cd8b1','Cd8a', 'Ifng','Gapdh','Il2rb', 'Ikzf2', 'Trac', 'Klre1', 'Klrk1', 'Bcl2', 'Klrc1', 'Icos', 'Id2','Ltb','Bcl2a1b','Lag3', 'Cxcr6','Ctla4', 'Prf1', 'Cd2', 'H2-Q7','Lck','Apobec3','Klra7', 'H2-T22','Mrc1', 'Tgfbi', 'Tgm2','Cxcl2','Ccl2','Ccl3','Ccl9','Ccl6')) # geneList3
  #geneList <-  unique(c('Ptger1','Ptger2','Ptger3','Ptger4','Ptges', 'Ptges2', 'Ptgs2', 'Ptgs1', 'Alox12b', 'Ptgis','Pla2g2e','Ggt1', 'Gpx3', 'Pla2g10','Pla2g12b','Ephx2', 'Alox12','Gpx1','Gpx4', 'Gpx2', 'Alox5', 'Alox15', 'Pla2g2d'))
  
  #geneList = c("Ptger1", "Ptger2", "Ptger3", "Ptger4", "Ptgs1", "Ptgs2", "Ptges", "Ptges2", "Lgr5", "Olfm4", "Lyz1", "Dclk1", "Myc", "Tcf", "Lef", "Ppar", "Eph", "Il18", "Il33", "S100a6", "Pdk4", "Fabp1", "Cpt1a", "Defa21", "Defa23", "Defa35", "Defa25", "Defa30", "Defa22", "Defa31", "Defa3", "Defa5", "Defa27", "Defa29", "Defa2", "Defa33", "Defa36", "Defa32", "Defa37", "Defa28", "Defa26", "Defa17", "Defa34", "Defa24", "Muc15", "Muc1", "Muc3", "Muc3a", "Muc6", "Muc2", "Muc5ac", "Muc5b", "Muc16", "Muc19", "Mucl1", "Mucl2", "Muc4", "Muc20", "Muc13")
  
  
  curPath = getwd()
  plotPath = paste0(curPath, "/GenePlots_", curFilename)
  dir.create(plotPath)
  setwd(plotPath)
  
  geneList = intersect(geneList, colnames(assay))
  dimredG = merge(dimredP, as.matrix(assay[, geneList, drop=F]), by.x="Cell", by.y="row.names")
  write.csv(dimredG, 'Organoid_vln_obj.csv')
  geneColors = rev(viridis::magma(10))#brewer.pal(9, "YlOrRd")[-c(1)] ## brewer.pal(9, "YlOrRd") ## rev(viridis::magma(10))
  geneColorsDiff = rev(brewer.pal(9, "RdBu"))
  
  #geneColors = rev(c("#552A35","#5B3A51","#544E6B","#3E647C","#237981","#268B77","#40A47B","#62BC7A","#8AD276","#B7E871","#EAFB6D"))
  for(curGene in geneList)
  {
    rangeL = quantile(dimredG[, curGene], 0.01, na.rm = T)
    rangeH = quantile(dimredG[, curGene], 0.95, na.rm = T)
    
    if(rangeL == rangeH) { rangeL = min(dimredG[, curGene]); rangeH = max(dimredG[, curGene]); }
    
    #ggplot(dimredG, aes_string(x = "UMAP_1", y = "UMAP_2", fill = curGene, color = curGene)) + geom_point(alpha=0.7) + theme_cem + scale_fill_gradientn(name = "Score", colors = geneColors, limits = c(rangeL, rangeH), oob = squish) + scale_color_gradientn(name = "Score", colors = geneColors, limits = c(rangeL, rangeH), oob = squish) + ggtitle(paste0(curGene)) + scale_alpha(range = c(0.3, 1))
    #ggsave(paste0(curFilename, "_dimred_geneExpression_umapAll_pointAll_", curGene, ".png"), width = dimredPlotWidth, height = dimredPlotHeight, dpi=150)
    #ggsave(paste0(curFilename, "_dimred_geneExpression_umapAll_pointAll_", curGene, ".pdf"), width = dimredPlotWidth, height = dimredPlotHeight)
    
    #ggplot(dimredG, aes_string(x = "UMAP_1", y = "UMAP_2", fill = curGene, color = curGene)) + geom_point(alpha=0.7) + theme_cem + scale_fill_gradientn(name = "Score", colors = geneColors, limits = c(rangeL, rangeH), oob = squish) + scale_color_gradientn(name = "Score", colors = geneColors, limits = c(rangeL, rangeH), oob = squish) + ggtitle(paste0(curGene)) + scale_alpha(range = c(0.3, 1)) + facet_grid(Type~.)
    #ggsave(paste0(curFilename, "_dimred_geneExpression_umapAll_pointByType_", curGene, ".png"), width = dimredPlotWidth, height = dimredPlotHeight*2+1, dpi=150)
    #ggsave(paste0(curFilename, "_dimred_geneExpression_umapAll_pointByType_", curGene, ".pdf"), width = dimredPlotWidth, height = dimredPlotHeight*2+1)
    
    sm = smoothScore2d(dimredG[, curGene], dimredG$UMAP_1, dimredG$UMAP_2, numGrid=100, knn=50, m=2)
    #ggplot(sm) + geom_tile(aes(x = x, y = y, fill = score)) + theme_cem + scale_fill_gradientn(name = "Score", colors = geneColors) + ggtitle(paste0(curGene)) + geom_polygon(data = hulls2, aes(x = x, y = y, group = CellType), alpha = 0.3, fill=NA, color="#666666") + geom_text(data = clusterMedian, aes(x=UMAP_1, y=UMAP_2, label = CellType, group = as.factor(CellType)), size = 5, color = "Black")
    #ggsave(paste0(curFilename, "_dimred_geneExpression_umapAll_densityAll_", curGene, ".png"), width = dimredPlotWidth, height = dimredPlotHeight, dpi=150)
    #ggsave(paste0(curFilename, "_dimred_geneExpression_umapAll_densityAll_", curGene, ".pdf"), width = dimredPlotWidth, height = dimredPlotHeight)
    
    
    sm2 = ddply(dimredG, "Type", function(x) { smoothScore2d( x[,  curGene], x$UMAP_1, x$UMAP_2, numGrid=100, knn=50, m=2, xrng=range(dimredG$UMAP_1), yrng=range(dimredG$UMAP_2)) } )	
    #ggplot(sm2) + geom_tile(aes(x = x, y = y, fill = score)) + theme_cem + scale_fill_gradientn(name = "Score", colors = geneColors) + ggtitle(paste0(curGene)) + geom_polygon(data = hulls2, aes(x = x, y = y, group = CellType), alpha = 0.3, fill=NA, color="#666666") + geom_text(data = clusterMedian, aes(x=UMAP_1, y=UMAP_2, label = CellType, group = as.factor(CellType)), size = 5, color = "Black") + facet_grid(~Type)
    #ggsave(paste0(curFilename, "_dimred_geneExpression_umapAll_densityByType_", curGene, ".png"), width = dimredPlotWidth*3-1, height = dimredPlotHeight, dpi=150)
    #ggsave(paste0(curFilename, "_dimred_geneExpression_umapAll_densityByType_", curGene, ".pdf"), width = dimredPlotWidth*3-1, height = dimredPlotHeight)	
    
    
    sm3 = smoothScore2d(dimredG[, curGene], dimredG$UMAP_1, dimredG$UMAP_2, type=dimredG$Type, numGrid=100, knn=50, m=2)
    sm3W = dcast(sm3, x+y~type, value.var="score")
    sm3W$CPT1aKO_vs_Control = sm3W$CPT1aKO - sm3W$Control
    rangeDiff = max(c(abs(sm3W$CPT1aKO_vs_Control)))
    
    ggplot(sm3W) + geom_tile(aes(x = x, y = y, fill = CPT1aKO_vs_Control)) +
      scale_fill_gradientn(name = "CPT1aKO vs Control", colors = geneColorsDiff, limits = c(-rangeDiff, rangeDiff)) +
      geom_polygon(data = hulls2, aes(x = x, y = y, group = CellType), alpha = 0.3, fill=NA, color="#666666", size=0.1) + ylab('UMAP_2') + xlab('UMAP_1') + theme_cem +
      theme(legend.position = "none", text = element_text(size=5),  axis.text.x = element_text(size =  6), axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7), panel.background = element_rect(fill = "white", colour = "Black"))
    
    ggsave(file= paste0(curFilename, "_dimred_geneExpression_umapAll_densityDifference_VP16vsControl_", curGene, ".pdf"), width=2.5, height=2.5, units="in")
  }
  
  setwd(curPath)
}


cell_types <- c('Stem 1', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft')

for(i in cell_types){
  cell = i
  dimredP <- dimredP %>% mutate(quantile = ntile(Pseudotime, 10))
  dimredaub = dimredP %>% filter(CellType == i)
  quart <- data.frame(var1 = tapply(dimredP$Pseudotime, dimredP$quantile, mean), var2 = tapply(dimredP$Pseudotime, dimredP$quantile, max))
  p.vals.prop = ''
  dimredc <- dimredaub %>% filter(Type == 'Control')
  dimreda <- dimredaub %>% filter(Type == 'CPT1aKO')
  
  for(i in 1:10){
    prop_table1 <- c(sum(dimreda$quantile == i), sum(dimredc$quantile == i))
    prop_table2 <- c(length(dimreda$quantile) - sum(dimreda$quantile == i), length(dimredc$quantile) - sum(dimredc$quantile == i)  )
    
    prop_table_fin <- cbind(prop_table1, prop_table2)
    p.vals.prop = c(p.vals.prop, prop.test(prop_table_fin)$p.value)
    p.vals.prop = p.vals.prop[2:11]
  }
  
  quart$p_asterik = ''
  quart$p_asterik[p.vals.prop < 0.05] = "*"
  quart$p_asterik[p.vals.prop < 0.01] = "**"
  quart$p_asterik[p.vals.prop < 0.001] = "***"
  quart$p_asterik[is.na(quart$p_asterik)] <- ""
  quart$p_asterik
  
  plot11 <- dimredaub %>% group_by(Type) %>% 
    # calculate densities for each group over same range; store in list column
    summarise(d = list(density(Pseudotime, from = min(.$Pseudotime), to = max(.$Pseudotime)))) %>% 
    # make a new data.frame from two density objects
    do(data.frame(x = .$d[[1]]$x,    # grab one set of x values (which are the same)
                  y = .$d[[2]]$y - .$d[[1]]$y))
  quart$maxi = max(plot11$y)
  
  ggplot(plot11, aes(x, y, fill=y>=0)) + geom_col(position = 'dodge', width = .0978619) + ggtitle(cell) + theme_cem + scale_fill_manual(values = c('#1b9e77', '#d95f02'))+ theme(legend.position = "none") + ylab('Density') + xlab('Pseudotime') + geom_vline(data = quart, aes(xintercept = var2), size = .075) + geom_text(data = quart, aes(var1 , maxi, label = p_asterik), position=position_dodge2(0.2), size = 1.5, vjust=1, inherit.aes = FALSE) + theme(text = element_text(size=1),axis.text = element_text(size = 2) ,plot.title = element_text(size = 6), axis.text.x = element_text(size = 3), axis.text.y = element_text(size = 3), axis.title.x = element_text(size = 4), axis.title.y = element_text(size = 4))
  ggsave(paste0("CPT1aKO_pseudo_density_", cell, ".pdf"), width = 2.2, height = 1, units="in")	
}


ggplot(dimredP, aes(x=Pseudotime, fill=Type)) + geom_density(alpha=0.4, bw=0.9) + theme_cem + facet_grid(CellType~., scales="free_y") + scale_fill_manual(values=unlist(colorsType)) + geom_vline(data = quart, aes(xintercept = var2), size = .25) + geom_text(data = quart, aes(var1 , maxi, label = p_asterik, fontface="bold"), position=position_dodge2(0.75), vjust=.5, inherit.aes = FALSE)
i = 'Stem 1'

cell = i
dimredP <- dimredP %>% mutate(quantile = ntile(Pseudotime, 10))
dimredaub = dimredP %>% filter(CellType == i)
quart <- data.frame(var1 = tapply(dimredP$Pseudotime, dimredP$quantile, mean), var2 = tapply(dimredP$Pseudotime, dimredP$quantile, max))
p.vals.prop = ''
dimredc <- dimredaub %>% filter(Type == 'Control')
dimreda <- dimredaub %>% filter(Type == 'CPT1aKO')

for(i in 1:10){
  prop_table1 <- c(sum(dimreda$quantile == i), sum(dimredc$quantile == i))
  prop_table2 <- c(length(dimreda$quantile) - sum(dimreda$quantile == i), length(dimredc$quantile) - sum(dimredc$quantile == i)  )
  
  prop_table_fin <- cbind(prop_table1, prop_table2)
  p.vals.prop = c(p.vals.prop, fisher.test(prop_table_fin)$p.value)
  p.vals.prop = p.vals.prop[2:11]
}

quart$p_asterik = ''
quart$p_asterik[p.vals.prop < 0.05] = "*"
quart$p_asterik[p.vals.prop < 0.01] = "**"
quart$p_asterik[p.vals.prop < 0.001] = "***"
quart$p_asterik[is.na(quart$p_asterik)] <- ""
quart$p_asterik

plot11 <- dimredaub  %>% group_by(Type) %>% 
  # calculate densities for each group over same range; store in list column
  summarise(d = list(density(Pseudotime, from = min(.$Pseudotime), to = max(.$Pseudotime)))) %>% 
  # make a new data.frame from two density objects
  do(data.frame(x = .$d[[1]]$x,    # grab one set of x values (which are the same)
                y = .$d[[2]]$y - .$d[[1]]$y))
quart$maxi = max(plot11$y)

ggplot(plot11, aes(x, y, fill=y>=0)) + geom_col(position = 'dodge', width = .0978619) + ggtitle('Stem 1') + theme_cem + scale_fill_manual(values = c('#1b9e77', '#d95f02'))+ theme(legend.position = "none") + ylab('Density') + xlab('Pseudotime') + geom_vline(data = quart, aes(xintercept = var2), size = .075) + geom_text(data = quart, aes(var1 , maxi, label = p_asterik), position=position_dodge2(0.2), size = 1.5, vjust=1, inherit.aes = FALSE) + theme(text = element_text(size=1),axis.text = element_text(size = 2) ,plot.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 0), axis.title.y = element_text(size = 0))
ggsave(paste0("CPT1aKO_pseudo_density_Stem_1.pdf"), width = 2.2, height = 1, units="in")	

ggplot(dimredP, aes(x=Pseudotime, fill=Type)) + geom_density(alpha=0.4, bw=0.5, size = .1) + theme_cem + 
  facet_grid(CellType~., scales="free_y") + scale_fill_manual(values=unlist(colorsType)) + 
  scale_y_continuous(expand = c(0.2,0), breaks = scales::breaks_extended(n = 3)) +
  geom_vline(data = quart, aes(xintercept = var2), size = .075) +
  theme(strip.background = element_blank(), strip.text.y = element_text(angle = 0, size = 6), text = element_text(size=1),axis.text = element_text(size = 2) ,plot.title = element_blank(), axis.text.x = element_text(size = 6), legend.position = "none", axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(.1, .1, 0, 0), "cm"))
ggsave(paste0("CPT1aKO_pseudo_density_all_cell_type.pdf"), width = 3, height = 4, units="in")	


{
  #stem 1
  i = 'Stem 1'
  #i = 'Enterocyte Progenitor'
  cell = i
  dimredP <- dimredP %>% mutate(quantile = ntile(Pseudotime, 10))
  dimredaub = dimredP %>% filter(CellType == i)
  quart <- data.frame(var1 = tapply(dimredP$Pseudotime, dimredP$quantile, mean), var2 = tapply(dimredP$Pseudotime, dimredP$quantile, max))
  p.vals.prop = ''
  dimredc <- dimredaub %>% filter(Type == 'Control')
  dimreda <- dimredaub %>% filter(Type == 'CPT1aKO')
  
  for(i in 1:10){
    prop_table1 <- c(sum(dimreda$quantile == i), sum(dimredc$quantile == i))
    prop_table2 <- c(length(dimreda$quantile) - sum(dimreda$quantile == i), length(dimredc$quantile) - sum(dimredc$quantile == i)  )
    
    prop_table_fin <- cbind(prop_table1, prop_table2)
    p.vals.prop = c(p.vals.prop, fisher.test(prop_table_fin)$p.value)
    p.vals.prop = p.vals.prop[2:11]
  }
  #p.vals.prop <- p.adjust(p.vals.prop, method = 'BH')
  quart$p_asterik = ''
  quart$p_asterik[p.vals.prop < 0.05] = "*"
  quart$p_asterik[p.vals.prop < 0.01] = "**"
  quart$p_asterik[p.vals.prop < 0.001] = "***"
  quart$p_asterik[is.na(quart$p_asterik)] <- ""
  quart$p_asterik
  
  plot11 <- dimredaub %>% group_by(Type) %>% 
    # calculate densities for each group over same range; store in list column
    summarise(d = list(density(Pseudotime, from = min(.$Pseudotime), to = max(.$Pseudotime)))) %>% 
    # make a new data.frame from two density objects
    do(data.frame(x = .$d[[1]]$x,    # grab one set of x values (which are the same)
                  y = .$d[[2]]$y - .$d[[1]]$y))
  quart$maxi = max(plot11$y)
  #quart$maxi2 = .005
  quart$maxi2 = .1
  Density_diff_stem_1 <- ggplot(plot11, aes(x, y, fill=y>=0)) + geom_col(position = 'dodge', width = .0978619) + ggtitle('Stem 1') + theme_cem + scale_fill_manual(values = c('#1b9e77', '#d95f02'))+ theme(legend.position = "none") + ylab('Density') + xlab('Pseudotime') + geom_vline(data = quart, aes(xintercept = var2), size = .075) + geom_text(data = quart, aes(var1 , maxi, label = p_asterik), position=position_dodge2(0.2), size = 1.5, vjust=1, inherit.aes = FALSE) + theme(text = element_text(size=1),axis.text = element_text(size = 2) ,plot.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(.3, 0, 0, 0), "cm"))
  Density_all_stem_1 <-ggplot(dimredaub, aes(x=Pseudotime, fill=Type)) + geom_density(alpha=0.6, size = .1) + 
    theme_cem + scale_fill_manual(values=unlist(colorsType)) + theme(legend.position = "none") + ylab('Density') + 
    xlab('Pseudotime') + geom_vline(data = quart, aes(xintercept = var2), size = .075) + geom_text(data = quart, aes(var1 , maxi2, label = p_asterik), position=position_dodge2(0.2), size = 1.5, vjust=.35, inherit.aes = FALSE) + 
    scale_y_continuous(expand = c(0.2,0), breaks = scales::breaks_extended(n = 6)) +
    theme(text = element_text(size=1),axis.text = element_text(size = 2) ,plot.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(.2, .1, 0, .2), "cm"))
  
  #all cells
  
  dimredaub = dimredP
  quart <- data.frame(var1 = tapply(dimredP$Pseudotime, dimredP$quantile, mean), var2 = tapply(dimredP$Pseudotime, dimredP$quantile, max))
  p.vals.prop = ''
  dimredc <- dimredaub %>% filter(Type == 'Control')
  dimreda <- dimredaub %>% filter(Type == 'CPT1aKO')
  
  for(i in 1:10){
    prop_table1 <- c(sum(dimreda$quantile == i), sum(dimredc$quantile == i))
    prop_table2 <- c(length(dimreda$quantile) - sum(dimreda$quantile == i), length(dimredc$quantile) - sum(dimredc$quantile == i)  )
    
    prop_table_fin <- cbind(prop_table1, prop_table2)
    p.vals.prop = c(p.vals.prop, fisher.test(prop_table_fin)$p.value)
    p.vals.prop = p.vals.prop[2:11]
  }
  p.vals.prop <- p.adjust(p.vals.prop, method = 'BH')
  quart$p_asterik = ''
  quart$p_asterik[p.vals.prop < 0.05] = "*"
  quart$p_asterik[p.vals.prop < 0.01] = "**"
  quart$p_asterik[p.vals.prop < 0.001] = "***"
  quart$p_asterik[is.na(quart$p_asterik)] <- ""
  quart$p_asterik
  
  plot11 <- dimredaub %>% group_by(Type) %>% 
    # calculate densities for each group over same range; store in list column
    summarise(d = list(density(Pseudotime, from = min(.$Pseudotime), to = max(.$Pseudotime)))) %>% 
    # make a new data.frame from two density objects
    do(data.frame(x = .$d[[1]]$x,    # grab one set of x values (which are the same)
                  y = .$d[[2]]$y - .$d[[1]]$y))
  quart$maxi = max(plot11$y)
  quart$maxi2 = .04
  
  Density_diff_all <- ggplot(plot11, aes(x, y, fill=y>=0)) + geom_col(position = 'dodge', width = .0978619) + ggtitle('All Cells') + theme_cem + scale_fill_manual(values = c('#1b9e77', '#d95f02'))+ theme(legend.position = "none") + ylab('Density') + xlab('Pseudotime') + geom_vline(data = quart, aes(xintercept = var2), size = .075) + geom_text(data = quart, aes(var1 , maxi, label = p_asterik), position=position_dodge2(0.2), size = 1.5, vjust=1, inherit.aes = FALSE) + theme(text = element_text(size=1),axis.text = element_text(size = 2) ,plot.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(.1, 0, 0, 0), "cm"))
  Density_all <-ggplot(dimredP, aes(x=Pseudotime, fill=Type)) + geom_density(alpha=0.6, size = .1) + 
    theme_cem + scale_fill_manual(values=unlist(colorsType)) + scale_y_continuous(expand = c(0.2,0), breaks = scales::breaks_extended(n = 5)) +
    theme(legend.position = "none") + ylab('Density') + xlab('Pseudotime') + geom_vline(data = quart, aes(xintercept = var2), size = .075) + geom_text(data = quart, aes(var1 , maxi2, label = p_asterik), position=position_dodge2(0.2), size = 1.5, vjust=.35, inherit.aes = FALSE) + theme(text = element_text(size=1),axis.text = element_text(size = 2) ,plot.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(.2, .1, 0, .2), "cm"))
  
  
  #stem 2
  i = 'Stem 2'
  cell = i
  dimredaub = dimredP %>% filter(CellType == i)
  quart <- data.frame(var1 = tapply(dimredP$Pseudotime, dimredP$quantile, mean), var2 = tapply(dimredP$Pseudotime, dimredP$quantile, max))
  p.vals.prop = ''
  dimredc <- dimredaub %>% filter(Type == 'Control')
  dimreda <- dimredaub %>% filter(Type == 'CPT1aKO')
  
  for(i in 1:10){
    prop_table1 <- c(sum(dimreda$quantile == i), sum(dimredc$quantile == i))
    prop_table2 <- c(length(dimreda$quantile) - sum(dimreda$quantile == i), length(dimredc$quantile) - sum(dimredc$quantile == i)  )
    
    prop_table_fin <- cbind(prop_table1, prop_table2)
    p.vals.prop = c(p.vals.prop, fisher.test(prop_table_fin)$p.value)
    p.vals.prop = p.vals.prop[2:11]
  }
  p.vals.prop <- p.adjust(p.vals.prop, method = 'BH')
  quart$p_asterik = ''
  quart$p_asterik[p.vals.prop < 0.05] = "*"
  quart$p_asterik[p.vals.prop < 0.01] = "**"
  quart$p_asterik[p.vals.prop < 0.001] = "***"
  quart$p_asterik[is.na(quart$p_asterik)] <- ""
  quart$p_asterik
  
  plot11 <- dimredaub %>% 
    # calculate densities for each group over same range; store in list column
    summarise(d = list(density(Pseudotime, from = min(.$Pseudotime), to = max(.$Pseudotime)))) %>% 
    # make a new data.frame from two density objects
    do(data.frame(x = .$d[[1]]$x,    # grab one set of x values (which are the same)
                  y = .$d[[1]]$y))
  quart$maxi = max(plot11$y)
  quart$maxi2 = 1.4
  Density_diff_stem_2 <- ggplot(plot11, aes(x, y, fill=y>=0)) + geom_col(position = 'dodge', width = .0978619) + theme_cem + scale_fill_manual(values = c('#d95f02'))+ 
    theme(legend.position = "none") + ylab('Density') + xlab('Pseudotime') + geom_vline(data = quart, aes(xintercept = var2), size = .075) + geom_text(data = quart, aes(var1 , maxi, label = p_asterik), position=position_dodge2(0.2), size = 1.5, vjust=1, inherit.aes = FALSE) + theme(text = element_text(size=1),axis.text = element_text(size = 2) ,plot.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(.1, 0, 0, 0), "cm"))
  Density_all_stem_2 <-ggplot(dimredaub, aes(x=Pseudotime, fill=Type)) + geom_density(alpha=0.6, size = .1) + 
    theme_cem + scale_fill_manual(values=unlist(colorsType)) + scale_y_continuous(expand = c(0.2,0), breaks = scales::breaks_extended(n = 5)) +
    theme(legend.position = "none") + ylab('Density') + xlab('Pseudotime') + geom_vline(data = quart, aes(xintercept = var2), size = .075) + geom_text(data = quart, aes(var1 , maxi, label = p_asterik), position=position_dodge2(0.2), size = 1.5, vjust=1, inherit.aes = FALSE) + 
    theme(text = element_text(size=1),axis.text = element_text(size = 2) ,plot.title = element_text(size = 6), axis.text.x = element_text(size = 6, hjust = 1), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(.2, 0, 0, .2), "cm"))
  
  
  plot_grid(Density_all, Density_diff_all, Density_diff_stem_1, Density_diff_stem_2, labels = c('','All', 'Stem 1', 'Stem 2'), label_size = 6, ncol = 1, align = "v", hjust = -.005) #'All Cells', 'Stem 1', 'Stem 2'
  
  ggsave(paste0("CPT1aKO_pseudo_density.pdf"), width = 2.2, height = 3.1, units="in")	
  
  plot_grid(Density_all, Density_all_stem_1, Density_all_stem_2, labels = c('All', 'Stem 1', 'Stem 2'), label_size = 6, ncol = 1, align = "v", hjust = -.02) #'All Cells', 'Stem 1', 'Stem 2'
  ggsave(paste0("CPT1aKO_pseudo_density_all.pdf"), width = 2.2, height = 3.1, units="in")	
  
  
}
# Plot the density difference
{
  library(MASS)
  library(reshape2)
  library(scales)
  
  # Calculate the common x and y range for geyser1 and geyser2
  xrng = range(dimred$UMAP_1)
  yrng = range(dimred$UMAP_2)
  
  extendRange = 0.06 # extend range by %6
  xDiff = (xrng[2] - xrng[1])
  yDiff = (yrng[2] - yrng[1])
  
  xrng[1] = xrng[1] - xDiff* extendRange
  xrng[2] = xrng[2] + xDiff * extendRange
  yrng[1] = yrng[1] - yDiff * extendRange
  yrng[2] = yrng[2] + yDiff * extendRange
  
  
  clusterMedian = dimredP %>%
    group_by(CellType) %>%
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
  
  hulls2 = concaveHull(dimred, "UMAP_1", "UMAP_2", "CellType", alpha = 0.5, extend = T, minPoint = 20)
  
  
  
  for(caseType in unique(dimred$Type))
  {
    for(ctrlType in unique(dimred$Type))
    {
      if(caseType == ctrlType) next
      
      caseVsCtrlName = paste0(caseType, "_vs_", ctrlType)
      
      colorLow = colorsType[ctrlType]
      colorMid = "white"
      colorHigh = colorsType[caseType]
      
      d_Case = kde2d( dimred$UMAP_1[dimred$Type == caseType], dimred$UMAP_2[dimred$Type == caseType ], lims = c(xrng, yrng), n = 500)
      d_Ctrl = kde2d( dimred$UMAP_1[dimred$Type == ctrlType  ], dimred$UMAP_2[dimred$Type == ctrlType ], lims = c(xrng, yrng), n = 500)
      
      # Confirm that the grid points for each density estimate are identical
      identical(d_Case$x, d_Ctrl$x) # TRUE
      identical(d_Case$y, d_Ctrl$y) # TRUE
      
      # Calculate the difference between the 2d density estimates
      diff_CaseVsCtrl = d_Ctrl 
      diff_CaseVsCtrl$z = d_Case$z - d_Ctrl$z
      diff_CaseVsCtrl$z = diff_CaseVsCtrl$z / max(diff_CaseVsCtrl$z)
      
      rownames(diff_CaseVsCtrl$z) = diff_CaseVsCtrl$x
      colnames(diff_CaseVsCtrl$z) = diff_CaseVsCtrl$y
      
      diff_CaseVsCtrlM = melt(diff_CaseVsCtrl$z, id.var = rownames(diff_CaseVsCtrl))
      names(diff_CaseVsCtrlM) = c("UMAP_1", "UMAP_2", caseVsCtrlName)
      
      ggplot(diff_CaseVsCtrlM, aes(x = UMAP_1, y = UMAP_2)) +
        ggrastr::rasterise(geom_tile(aes_string(fill = caseVsCtrlName), alpha = 1), dpi = 200) + #use hannah's raster code here GGrastr
        scale_fill_gradient2(low = colorLow, mid = colorMid, high = colorHigh, midpoint = 0) +
        coord_cartesian(xlim = xrng, ylim = yrng) +
        scale_color_manual(values = colorsType) +
        guides(colour = FALSE) + theme_cem +
        geom_polygon(data = hulls2, aes(x = x, y = y, group = CellType), alpha = 0.1, size=0.1, fill=NA, color="Grey50") + 
        geom_text_repel(data = clusterMedian, aes(label = CellType, group = as.factor(CellType)), size = 2, box.padding = .1, force = 5, color = "Black") + 
        ggtitle(paste0("Density comparison of ", caseType, " vs ", ctrlType)) + theme(text = element_text(size=5), legend.key.size = unit(0.0, "cm"), legend.text=element_text(size=0), legend.title =element_text(size=0), axis.text.x = element_text(size = 6), plot.title = element_text(size = 0), axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7))
      
      ggsave(paste0(curFilename, "_dimred_densityDiff_cellType_", caseVsCtrlName, ".pdf"), width = 2.5, height = 2.5, units="in")
    }
  }
}



## Trajectory analysis
{
  
  library(Seurat)
  library(scales)
  library(dplyr)
  library(dyno)
  library(tidyverse)
  library(Matrix)
  library(pheatmap)
  
  
  seuratExp = t(seuratObj@assays$SCT@scale.data)
  #seuratExp = t(seuratObj@assays$integrated@scale.data)
  
  seuratCt = t(seuratObj@assays$SCT@counts)
  commonCells = intersect(rownames(seuratExp), rownames(seuratCt))
  commonGenes = intersect(colnames(seuratExp), colnames(seuratCt))
  seuratExp2 = seuratExp[ commonCells, commonGenes ]
  seuratCt2 = seuratCt[ commonCells, commonGenes ]
  
  dataset = wrap_expression(
    counts = seuratCt2,
    expression = seuratExp2
  )
  dynwrap::test_docker_installation(detailed = TRUE)
  ti_methods <- get_ti_methods()
  model = infer_trajectory(dataset, ti_slingshot(), verbose = TRUE)
  model_CPT1aKO = model
  model <- readRDS('CPT1aKO_pseudotime.rds', refhook = NULL)
  overall_feature_importances = dynfeature::calculate_overall_feature_importance(model, expression_source = dataset$expression)
  
  topN = 100
  features_overall = overall_feature_importances %>% 
    top_n(topN, importance) %>% 
    pull(feature_id)
  
  groupingCluster = factor(dimred$Cluster)
  names(groupingCluster) = dimred$Cell
  
  groupingType = factor(dimred$Type)
  names(groupingType) = dimred$Cell
  
  identities <- levels(CPT1aKO.obj$Cell_type)
  my_color_palette <- hue_pal()(length(identities))
  
  plot_heatmap(
    model, 
    expression_source = dataset$expression, 
    features_oi = features_overall,
    grouping = groupingCluster,
    
  ) + scale_fill_manual(values = my_color_palette)
  
  ggsave(paste0(curFilename, "_heatmapFeatureOverall.png"), width=12, height=14, dpi=150)
  ggsave(paste0(curFilename, "_heatmapFeatureOverall.pdf"), width=12, height=14)
  
  
  
  pseudotime = dynwrap::calculate_pseudotime(model)
  features_overall = overall_feature_importances[ !grepl("^Rp", overall_feature_importances$feature_id), ] %>% top_n(topN, importance)  %>% pull(feature_id)
  features_pseudo_genes <- c( "S100a6", "Ascl2", "Lgr5", "Ly6a")
  overall_feature_impartace <- overall_feature_importances[ !grepl("^Rp", overall_feature_importances$feature_id),]
  features_pseudo_genes = overall_feature_impartace[overall_feature_impartace$feature_id %in% features_pseudo_genes,]  %>% pull(feature_id)
  
  features_overall = c(features_overall, features_pseudo_genes)
  
  seuratExpSel = t(seuratExp[ names(pseudotime)[ order(pseudotime) ], features_overall ] )
  rownames(seuratExpSel)
  geneOrder = getTSPOrder(seuratExpSel)
  seuratExpSel = seuratExpSel[geneOrder, ]
  #seuratExpSelControl = seuratExpSel[, colnames(seuratExpSel) %in% dimred$Cell[ dimred$Type == "Control"] ]
  #seuratExpSelCPT1aKO = seuratExpSel[, colnames(seuratExpSel) %in% dimred$Cell[ dimred$Type == "CPT1aKO"] ]
  
  
  annColumn = data.frame(row.names=names(pseudotime), Pseudotime=pseudotime, Type=splitGet(names(pseudotime), "_", 2))
  annColors = list(Pseudotime = brewer.pal(11, "Spectral"), Type = colorsType)
  
  colScale = viridis::inferno(50)
  maxVal = max(seuratExpSel)
  minVal = min(seuratExpSel)
  breaks = c(minVal, seq(quantile(seuratExpSel, 0.2), quantile(seuratExpSel, 0.99), length.out = 49), maxVal)
  
  heatmapWidth = 18
  heatmapHeight = 13
  
  
  ######
  
  pheatmap(seuratExpSel, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn, filename = paste0(curFilename, "_heatmapByPseudotime_byCells_all.png"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nAll Cells")
  pheatmap(seuratExpSel, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn, filename = paste0(curFilename, "_heatmapByPseudotime_byCells_all.pdf"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nAll Cells")
  
  pheatmap(seuratExpSelControl, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn,  filename = paste0(curFilename, "_heatmapByPseudotime_byCells_Control.png"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nControl Cells")
  pheatmap(seuratExpSelControl, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn,  filename = paste0(curFilename, "_heatmapByPseudotime_byCells_Control.pdf"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nControl Cells")
  
  pheatmap(seuratExpSelCPT1aKO, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn, filename = paste0(curFilename, "_heatmapByPseudotime_byCells_CPT1aKO.png"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nCPT1aKO Cells")
  pheatmap(seuratExpSelCPT1aKO, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn,  filename = paste0(curFilename, "_heatmapByPseudotime_byCells_CPT1aKO.pdf"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nCPT1aKO Cells")
  
  
  ######
  
  
  source("~/analysis/scRNA_Intestine/Single_Cell_Functions/cem_polyOutline.R")
  
  dimredP = merge(dimred, annColumn[, c("Pseudotime"), drop=F], by.x="Cell", by.y="row.names")
  
  clusterMedian = dimred %>%
    group_by(Cluster) %>%
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
  clusterMedian = as.data.frame(clusterMedian)
  hulls = concaveHull(dimred, "UMAP_1", "UMAP_2", "Cluster", alpha = 0.5, extend = T, minPoint = 3)
  
  
  
  ddimred = data.frame(row.names=dimred$Cell, comp_1=dimred$UMAP_1, comp_2=dimred$UMAP_2)
  cell_positions = ddimred %>% as.data.frame() %>% rownames_to_column("cell_id")
  cell_positions = as.data.frame(cell_positions)
  cell_positions = left_join(cell_positions, model$milestone_percentages %>% 
                               group_by(cell_id) %>% arrange(desc(percentage)) %>% 
                               filter(dplyr::row_number() == 1) %>% select(cell_id, milestone_id), "cell_id")
  color_trajectory = "none"	
  waypoints = dynwrap::select_waypoints(model)
  trajectory_projection_sd = sum(model$milestone_network$length) * 0.05
  waypoint_projection = dynplot:::project_waypoints(trajectory = model, 
                                                    cell_positions = cell_positions, waypoints = waypoints, 
                                                    trajectory_projection_sd = trajectory_projection_sd, 
                                                    color_trajectory = color_trajectory)
  
  arrow = if (any(model$milestone_network$directed))	arrow(type = "closed", length = (unit(0.1, "inches")))
  
  edges = as.data.frame(waypoint_projection$edges)
  edgesArrow = edges[ which(edges$arrow), ]
  ggplot(dimredP, aes(x = UMAP_1, y = UMAP_2)) + geom_point(aes(fill = Pseudotime, color = Pseudotime)) + theme_cem + scale_fill_gradientn(colors = annColors$Pseudotime) + scale_color_gradientn(colors = annColors$Pseudotime) + ggtitle(paste0(curTitle, ", UMAP by Pseudotime")) + geom_text(data = clusterMedian, aes(label = Cluster, group = as.factor(Cluster)), size = 5, color = "Black") + geom_polygon(data = hulls, aes(x = x, y = y, group = Cluster), alpha = 0.2, fill=NA, color="grey50") + geom_segment(aes(comp_1_from, comp_2_from, xend = comp_1_to, yend = comp_2_to), data = edgesArrow , arrow = arrow, color = "#333333", size = 1, linejoin = "mitre", lineend = "butt") + geom_segment(aes(comp_1_from, comp_2_from, xend = comp_1_to, yend = comp_2_to), data = edges, size = 1, color = "#333333")
  ggsave(paste0(curFilename, "_dimred_pseudotime_hull.png"), width = dimredPlotWidth, height = dimredPlotHeight)
  ggsave(paste0(curFilename, "_dimred_pseudotime_hull.pdf"), width = dimredPlotWidth, height = dimredPlotHeight)
  
  ggplot(dimredP, aes(x = UMAP_1, y = UMAP_2)) + geom_point(aes(fill = Pseudotime, color = Pseudotime)) + theme_cem + scale_fill_gradientn(colors = annColors$Pseudotime) + scale_color_gradientn(colors = annColors$Pseudotime) + ggtitle(paste0(curTitle, ", UMAP by Pseudotime")) + geom_segment(aes(comp_1_from, comp_2_from, xend = comp_1_to, yend = comp_2_to), data = edgesArrow, arrow = arrow, color = "#333333", size = 1, linejoin = "mitre", lineend = "butt") + geom_segment(aes(comp_1_from, comp_2_from, xend = comp_1_to, yend = comp_2_to), edges, size = 1, color = "#333333")
  ggsave(paste0(curFilename, "_dimred_pseudotime.png"), width = dimredPlotWidth, height = dimredPlotHeight)
  ggsave(paste0(curFilename, "_dimred_pseudotime.pdf"), width = dimredPlotWidth, height = dimredPlotHeight)
  
  write.tsv(dimredP, "_dimred_pseudotime.txt", row.names=F)
  
  
  ############
  
  
  
  #seuratExpSel2 = seuratExpSel - rep.col(rowMins(seuratExpSel), ncol(seuratExpSel))
  seuratExpSelMed = as.data.frame(seuratExpSel)
  seuratExpSelMed$Gene = rownames(seuratExpSelMed)
  seuratExpSelMed = melt(seuratExpSelMed, "Gene", variable.name = "Cell", value.name = "Expression")
  seuratExpSelMed = merge(seuratExpSelMed, data.frame(Cell=names(pseudotime), Pseudotime=pseudotime), by="Cell")
  seuratExpSelMed = merge(seuratExpSelMed, dimred[, c("Cell", "Type")], by="Cell")
  #seuratExpSelMed$PseudotimeFactor = cut(seuratExpSelMed$Pseudotime, 100)
  seuratExpSelMed$PseudotimeFactor = Hmisc::cut2(seuratExpSelMed$Pseudotime, g=100)
  seuratExpSelMedSum1 = seuratExpSelMed %>% group_by(Gene, PseudotimeFactor) %>% summarise(Expression=median(Expression, na.rm=T))
  seuratExpSelMedSum2 = seuratExpSelMed %>% group_by(Gene, PseudotimeFactor, Type) %>% summarise(Expression=median(Expression, na.rm=T))
  seuratExpSelMedSumW = as.matrix(setRemoveRownames(dcast(seuratExpSelMedSum1, Gene~PseudotimeFactor, value.var="Expression")))
  seuratExpSelMedSumW_Ctrl = as.matrix(setRemoveRownames(dcast(seuratExpSelMedSum2[ seuratExpSelMedSum2$Type == "Control", ], Gene~PseudotimeFactor, value.var="Expression")))
  seuratExpSelMedSumW_Case = as.matrix(setRemoveRownames(dcast(seuratExpSelMedSum2[ seuratExpSelMedSum2$Type == "CPT1aKO", ], Gene~PseudotimeFactor, value.var="Expression")))
  
  cellCounts = seuratExpSelMed %>% dplyr::select(Cell, Type, PseudotimeFactor) %>% distinct()  #%>% group_by(Type, PseudotimeFactor) %>% summarise(CellCount = n())
  cellCounts$Pseudotime = as.numeric(gsub("[\\[\\(\\ ]", "", splitGet(as.character(cellCounts$PseudotimeFactor), ",", 1)))
  
  #clust = hclust(as.dist( 1 - (cor(seuratExpSelMedSumW)+1)/2 ), method = "ward.D2")
  #clust = hclust(dist(seuratExpSelMedSumW), method = "ward.D2")
  #geneOrder = clust$order
  
  #geneOrder = getTSPOrder(seuratExpSelMedSumW)
  seuratExpSelMedSumW = seuratExpSelMedSumW[geneOrder, ]
  seuratExpSelMedSumW_Ctrl = seuratExpSelMedSumW_Ctrl[geneOrder, ]
  seuratExpSelMedSumW_Case = seuratExpSelMedSumW_Case[geneOrder, ]
  seuratExpSelMedSumW_Diff = seuratExpSelMedSumW_Case - seuratExpSelMedSumW_Ctrl
  
  annColumn2 = setRemoveRownames(unique(cellCounts[, c("PseudotimeFactor", "Pseudotime")]), 1)
  
  colScale = viridis::inferno(50)
  maxVal = max(seuratExpSelMedSumW)
  minVal = min(seuratExpSelMedSumW)
  breaks = c(minVal, seq(quantile(seuratExpSelMedSumW, 0.2), quantile(seuratExpSelMedSumW, 0.99), length.out = 49), maxVal)
  
  pheatmap(seuratExpSelMedSumW, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn2, annotation_colors = annColors, filename = paste0(curFilename, "_heatmapByPseudotime_medianExpression_all.png"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nMedian expression across pseudotime, All Cells")
  pheatmap(seuratExpSelMedSumW, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn2, annotation_colors = annColors, filename = paste0(curFilename, "_heatmapByPseudotime_medianExpression_all_prominent.pdf"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of genes by pseudotime importance\nMedian expression across pseudotime, All Cells")
  
  pheatmap(seuratExpSelMedSumW_Case, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn2, annotation_colors = annColors, filename = paste0(curFilename, "_heatmapByPseudotime_medianExpression_CPT1aKO.png"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nMedian expression across pseudotime, CPT1aKO Cells")
  pheatmap(seuratExpSelMedSumW_Case, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn2, annotation_colors = annColors, filename = paste0(curFilename, "_heatmapByPseudotime_medianExpression_CPT1aKO.pdf"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nMedian expression across pseudotime, CPT1aKO")
  
  
  pheatmap(seuratExpSelMedSumW_Ctrl, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn2, annotation_colors = annColors, filename = paste0(curFilename, "_heatmapByPseudotime_medianExpression_Control.png"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nMedian expression across pseudotime, Control Cells")
  pheatmap(seuratExpSelMedSumW_Ctrl, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn2, annotation_colors = annColors, filename = paste0(curFilename, "_heatmapByPseudotime_medianExpression_Control.pdf"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nMedian expression across pseudotime, Control Cells")
  
  
  
  colScale = colorRampPalette(rev(c("#e66101", "#fdb863", "#f7f7f7", "#b2abd2", "#5e3c99")))(50)
  maxVal = max(abs(seuratExpSelMedSumW_Diff))
  breaks = c(-maxVal, seq(-2, 2, length.out = 49), maxVal)
  
  annColumn3 = annColumn2[, "Pseudotime", drop = F]
  pheatmap(seuratExpSelMedSumW_Diff, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn3, annotation_colors = annColors, filename = paste0(curFilename, "_heatmapByPseudotime_medianExpressionDifference_CPT1aKOvsCtrl.png"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nCPT1aKO - Control median expression difference across pseudotime")
  pheatmap(seuratExpSelMedSumW_Diff, scale="none", cluster_rows = F, cluster_cols = F, show_colnames = F, breaks = breaks, treeheight_row = 0, treeheight_col = 0, color = colScale, annotation_col = annColumn3, annotation_colors = annColors, filename = paste0(curFilename, "_heatmapByPseudotime_medianExpressionDifference_CPT1aKOvsCtrl.pdf"), width=heatmapWidth, height=heatmapHeight, border_color = NA, main = "Heatmap of top 100 genes by pseudotime importance\nCPT1aKO - Control median expression difference across pseudotime")
  
  ggplot(cellCounts, aes(x=Pseudotime, fill=Type, color=Type)) + geom_density(alpha=0.7, adjust=0.5) + theme_cem + scale_fill_manual(values=colorsType) + scale_color_manual(values=colorsType)
  ggsave(paste0(curFilename, "_densityByPseudotimeFactor.pdf"),  width=heatmapWidth, height=heatmapHeight)
  ggsave(paste0(curFilename, "_densityByPseudotimeFactor.png"),  width=heatmapWidth, height=heatmapHeight, dpi=150)
  
  ggplot(annColumn, aes(x=Pseudotime, fill=Type, color=Type)) + geom_density(alpha=0.7, adjust=0.5) + theme_cem + scale_fill_manual(values=colorsType) + scale_color_manual(values=colorsType)
  ggsave(paste0(curFilename, "_densityByPseudotime.pdf"),  width=heatmapWidth, height=heatmapHeight)
  ggsave(paste0(curFilename, "_densityByPseudotime.png"),  width=heatmapWidth, height=heatmapHeight, dpi=150)
  
}


