# Vyom Shah, November 2020
# Expression analysis of single cell RNA Sequencing data for IL17 genes in various cell types
# Project title: "IL-17RA signaling via Lgr5 stem cells and progenitors regulates intestinal secretory cell lineage commitment"

Idents(Gene_de) <- Gene_de$orig.ident
SBO4 <- subset(Gene_de, idents =c('Control4'))
Idents(SBO4) <- SBO4$Cell_type
immune <- subset(SBO4, idents =c("Stem 1", "Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying", "Enterocyte Progenitor" , "Enterocyte", "Enteroendocrine", "Goblet", "Paneth" ))

current.cluster.ids <- c("Stem 1","Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying", "Enterocyte Progenitor" , "Enterocyte", "Enteroendocrine", "Goblet", "Paneth" )
new.cluster.ids <- c("Stem","Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying", "Enterocyte Progenitor" , "Enterocyte", "Enteroendocrine", "Goblet", "Paneth")
immune$Cell_type <- plyr::mapvalues(x = immune@active.ident, from = current.cluster.ids, to = new.cluster.ids)

DotPlot_Sig <- c('Il17ra','Il17rc','Il17rb',' Il17re','Il17rd')
DotPlot(immune, features = DotPlot_Sig, assay = 'RNA', group.by = 'Cell_type') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(y= "Cell Type", x="") + scale_colour_distiller( palette ="RdYlBu")

Il17ra <- VlnPlot(immune, group.by = 'Cell_type', assay = 'RNA', features = c('Il17ra'), ncol = 1, log = TRUE, pt.size = 0) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
Il17rc <- VlnPlot(immune, group.by = 'Cell_type', assay = 'RNA', features = c('Il17rc'), ncol = 1, log = TRUE, pt.size = 0) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
Il17rb <- VlnPlot(immune, group.by = 'Cell_type', assay = 'RNA', features = c('Il17rb'), ncol = 1, log = TRUE, pt.size = 0) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')
Il17rd <- VlnPlot(immune, group.by = 'Cell_type', assay = 'RNA', features = c('Il17rd'), ncol = 1, log = TRUE, pt.size = 0) + geom_boxplot(width=0.1, fill='#A4A4A4',outlier.shape = NA) + theme(legend.position = 'none')

plot_grid(Il17ra, Il17rb, Il17rc, Il17rd)
