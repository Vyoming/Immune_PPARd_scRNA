library(nichenetr)
library(Seurat)
library(tidyverse)
library(circlize)

Immune_VP16@meta.data %>% head()

Immune_VP16@meta.data$Type %>% table()

# load in Nichenet's data
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
head(lr_network)

weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network

head(weighted_networks$gr) # interactions and their weights in the gene regulatory network


lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()

ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()


## receiver
levels(Immune_VP16$Cell_Type)
levels(Immune_VP16$Type)
Idents(Immune_VP16) <- Immune_VP16$Cell_Type
DefaultAssay(Immune_VP16) <-  "RNA"
receiver = c("CD8+ T Cell",'Ex CD8+ T Cell')
expressed_genes_receiver = get_expressed_genes(receiver, Immune_VP16, pct = 0.10)

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## sender
sender_celltypes = c("Neutrophil","TReg", "M1 Macrophage", "Proliferating TReg", "M2 Macrophage", "B Cell",'DC','NK Cell','Monocyte')

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, Immune_VP16, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

seurat_obj_receiver= subset(Immune_VP16, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["Type"]])

condition_oi = "VP16PPARd"
condition_reference = "Control" 

DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_logFC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

#see which cell population expresses which of these top-ranked ligands
DotPlot(Immune_VP16, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()

#identify active target gene inference
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network

#Receptors of top ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

#Receptors of top-ranked ligands, but after considering only bona fide ligand-receptor interactions documented in literature and publicly available databases

lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))

lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()

p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
p_ligand_receptor_network_strict


#Add log fold change information of ligands from sender cells
# DE analysis for each sender cell type
# this uses a new nichenetr function - reinstall nichenetr if necessary!
DE_table_all = Idents(Immune_VP16) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seurat_obj = Immune_VP16, condition_colname = "Type", condition_oi = condition_oi, condition_reference = condition_reference, celltype_col = NULL, expression_pct = 0.10) %>% reduce(full_join)
DE_table_all[is.na(DE_table_all)] = 0

# Combine ligand activities with DE information
ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0

# make LFC heatmap
lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,]

colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()
p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
p_ligand_lfc

order_receptors
order_ligands_receptor
best_upstream_receptors
#run the same for receiver cells
DE_table_all = Idents(Immune_VP16) %>% levels() %>% intersect(receiver) %>% lapply(get_lfc_celltype, seurat_obj = Immune_VP16, condition_colname = "Type", condition_oi = condition_oi, condition_reference = condition_reference, celltype_col = NULL, expression_pct = 0.10) %>% reduce(full_join)
DE_table_all[is.na(DE_table_all)] = 0

# Combine ligand activities with DE information
receptor_probability <- data.frame(vis_ligand_receptor_network_strict)
receptor_probability <-  receptor_probability[c('Cxcr3', 'Ifngr1', 'Ctla4','Cd28','Pdcd1'),]
receptor_probability$receptors <- rownames(receptor_probability)
receptor_activities_de = receptor_probability %>% select(receptors) %>% rename(receptor = receptors) %>% left_join(DE_table_all %>% rename(receptor = gene))
receptor_activities_de[is.na(receptor_activities_de)] = 0

# make LFC heatmap
lfc_matrix = receptor_activities_de  %>% select(-receptor) %>% as.matrix() %>% magrittr::set_rownames(receptor_activities_de$receptor)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix

colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()
p_receptor_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized Receptors","LFC in Receiver", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"), legend.text = element_text(size = 9)) + coord_flip()
p_receptor_lfc

#summary of aforementioned visualizations in the Nichenet analysis

# ligand activity heatmap
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)


rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
order_ligands <- c('Cxcl10', 'Ifng', 'Cd86','Icam1','Ccl12')
vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")

p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9)) 

# ligand expression Seurat dotplot
order_ligands_adapted = order_ligands
order_ligands_adapted[order_ligands_adapted == "H2.M3"] = "H2-M3" # cf required use of make.names for heatmap visualization | this is not necessary if these ligands are not in the list of prioritized ligands!
order_ligands_adapted[order_ligands_adapted == "H2.T23"] = "H2-T23" # cf required use of make.names for heatmap visualization | this is not necessary if these ligands are not in the list of prioritized ligands!
rotated_dotplot = DotPlot(Immune_VP16 %>% subset(Cell_Type %in% sender_celltypes), features = order_ligands_adapted, cols = "RdYlBu") + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) # flip of coordinates necessary because we want to show ligands in the rows when combining all plots

#combine everything together
figures_without_legend = cowplot::plot_grid(
  p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
  p_ligand_lfc + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
  p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
  align = "hv",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_pearson)+6, ncol(vis_ligand_lfc) + 7, ncol(vis_ligand_lfc) + 8, ncol(vis_ligand_target)))

legends = cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),
  ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_lfc)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
  nrow = 1,
  align = "h", rel_widths = c(1.5, 1, 1, 1))

combined_plot = cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
combined_plot

#Circos Plot to show how everything reacts together

sender_celltypes = c("Neutrophil","TReg", "M1 Macrophage", "Proliferating TReg", "M2 Macrophage", "B Cell",'DC','NK Cell','Monocyte')

Idents(Immune_VP16) <- Immune_VP16$Cell_Type
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = Immune_VP16, 
  receiver = c("CD8+ T Cell",'Ex CD8+ T Cell'), 
  condition_colname = "Type", condition_oi = "VP16PPARd", condition_reference = "Control", 
  sender = sender_celltypes,
  assay_oi = 'RNA',
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, 
  organism = "mouse")


avg_expression_ligands = AverageExpression(Immune_VP16, features = best_upstream_ligands)
sender_ligand_assignment = avg_expression_ligands$RNA %>% apply(1, function(ligand_expression){
  ligand_expression > (ligand_expression %>% mean() + ligand_expression %>% sd())
}) %>% t()
typeof(sender_ligand_assignment)
sender_ligand_assignment['Cxcl10','TReg'] = TRUE
sender_ligand_assignment['Cxcl10','Proliferating TReg'] = TRUE
sender_ligand_assignment['Cxcl10','DC'] = TRUE
sender_ligand_assignment = sender_ligand_assignment %>% apply(2, function(x){x[x == TRUE]}) %>% purrr::keep(function(x){length(x) > 0})
names(sender_ligand_assignment)

all_assigned_ligands = sender_ligand_assignment %>% lapply(function(x){names(x)}) %>% unlist()
unique_ligands = all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
general_ligands = best_upstream_ligands %>% setdiff(unique_ligands)
sender_celltypes = c("Neutrophil","TReg", "M1 Macrophage", "Proliferating TReg", "M2 Macrophage", "B Cell",'DC','NK Cell','Monocyte')
B_specific_ligands = sender_ligand_assignment[['B Cell']] %>% names() %>% setdiff(general_ligands)
Ptreg_specific_ligands = sender_ligand_assignment[['Proliferating TReg']] %>% names() %>% setdiff(general_ligands)
Treg_specific_ligands = sender_ligand_assignment[['TReg']] %>% names() %>% setdiff(general_ligands)
M1_specific_ligands = sender_ligand_assignment[['M1 Macrophage']] %>% names() %>% setdiff(general_ligands)
M2_specific_ligands = sender_ligand_assignment[['M2 Macrophage']] %>% names() %>% setdiff(general_ligands)
DC_specific_ligands = sender_ligand_assignment[['DC']] %>% names() %>% setdiff(general_ligands)
NK_specific_ligands = sender_ligand_assignment[['NK Cell']] %>% names() %>% setdiff(general_ligands)
Monocyte_specific_ligands = sender_ligand_assignment[['Monocyte']] %>% names() %>% setdiff(general_ligands)
Neutr_specific_ligands = sender_ligand_assignment[['Neutrophil']] %>% names() %>% setdiff(general_ligands)


ligand_type_indication_df = tibble(
  ligand_type = c(rep("B-specific", times = B_specific_ligands %>% length()),
                  rep("NK-specific", times = NK_specific_ligands %>% length()),
                  rep("Mono-specific", times = Monocyte_specific_ligands %>% length()),
                  rep("DC-specific", times = DC_specific_ligands %>% length()),
                  rep("Ptreg_specific", times = Ptreg_specific_ligands %>% length()),
                  rep("Treg_specific", times = Treg_specific_ligands %>% length()),
                  rep("M1_specific", times = M1_specific_ligands %>% length()),
                  rep("M2_specific", times = M2_specific_ligands %>% length()),
                  rep("Neutr_specific", times = Neutr_specific_ligands %>% length()),
                  rep("General", times = general_ligands %>% length())),
  ligand = c(B_specific_ligands, NK_specific_ligands, Monocyte_specific_ligands, DC_specific_ligands, Ptreg_specific_ligands, Treg_specific_ligands, M1_specific_ligands, M2_specific_ligands, Neutr_specific_ligands, general_ligands))
 
