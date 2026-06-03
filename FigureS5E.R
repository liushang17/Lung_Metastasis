library(tidyverse)
library(Seurat)
library(ggplot2)
library(nichenetr)

args<-commandArgs(T)

ligand_target_matrix = readRDS("ligand_target_matrix.rds") # target genes in rows, ligands in columns
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols() 
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols() 
ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]


#1 Define expressed genes in sender and receiver cell populations:至少 10% 的细胞中表达
## receiver
objtcell <- readRDS("Tumor cell.rds")
Idents(objtcell) <- objtcell$celltype
receiver = ct
expressed_genes_receiver = get_expressed_genes(receiver, objtcell, pct = 0.1)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
length(background_expressed_genes)
## sender
tme <- readRDS("TME.rds")
Idents(tme) <- tme$orig.ident
sender_celltypes = unique(as.character(tme$orig.ident))
expressed_genes_sender = get_expressed_genes(sender_celltypes, tme, pct = 0.02)

#2 Define the gene set of interest and a background of genes
geneset <- read.table("Spatial.gene.xls",sep = "\t",header = T)
pos <- which(geneset$Status %in% "Phase2")
geneset1 <- geneset[pos,]

#3 Define a set of potential ligands
lr_network = readRDS("lr_network.rds")
lr_network$from <- lr_network$from %>% convert_human_to_mouse_symbols() 
lr_network$to <- lr_network$to %>% convert_human_to_mouse_symbols() 
pos <- which(lr_network$from %in% NA | lr_network$to %in% NA)
lr_network <- lr_network[-pos,]

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)
receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)
lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)
potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
head(potential_ligands)   #从数据库中筛出所有包括在>10%表达量gene中的geneligand-receptor对

#4 ligand activity analysis:predict whether a gene belongs to the interest program or not,person相关系数
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities %>% arrange(-pearson) 
# show histogram of ligand activity scores
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
ggsave(paste0(cl,'_',ct,'_lig_activity_hist.pdf'),plot = p_hist_lig_activity,width = 5,height = 5)
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)

#5 ligand推测：Infer target genes:regulatory potential scores between ligands and target genes of interest
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()  #250 most strongly predicted targets
active_ligand_target_links_df <- na.omit(active_ligand_target_links_df)
nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)
nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)
#配体活性预测heatmap:推断排名靠前的靶基因的配体,配体和靶基因之间的潜在调节评分
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
order_targets <- intersect(order_targets,rownames(active_ligand_target_links))
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized Tumor-ligands","migrition genes in T cells", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))
ggsave(paste0(cl,'_',ct,'_ligand_target_network_heat.pdf'),p_ligand_target_network,width = 7,height = 5)

#6 receptor推测：Ligand-receptor network inference for top-ranked ligands
# get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
# get the weights of the ligand-receptor interactions as used in the NicheNet model
weighted_networks = readRDS("weighted_networks.rds")

weighted_networks_lr <- weighted_networks$lr_sig
weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()



lr_network_top_df = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
# convert to a matrix
lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
# perform hierarchical clustering to order the ligands and receptors
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]
dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]  #ligand重新排序
#显示配体-受体相互作用调节评分的热图
vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Prioritized Tumor-ligands","Receptors expressed by T cells", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
ggsave(paste0(cl,'_',ct,'_ligand_receptor_network_heat.pdf'),p_ligand_receptor_network,width = 7,height = 5)

#可视化top预测ligand及其靶基因的表达:显示ligand活性、ligand表达、靶基因表达和ligand-target gene调节潜力的组合图
library(RColorBrewer)
library(cowplot)
library(ggpubr)
#准备ligand活性矩阵
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized Tumor-ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)")
ggsave(paste0(cl,'_ligand_activity_pearson.pdf'),p_ligand_pearson,width = 7,height = 5)

#配体的表达矩阵:显示每个样本中配体的平均表达
expression_tumor <- as.matrix(objtumor@assays$RNA@data)
metatumor <- objtumor@meta.data
colnames(metatumor)[2] <- 'cell'
expression_df_Tumor = t(expression_tumor)[,order_ligands] %>% data.frame() %>% rownames_to_column("cell") %>% as_tibble() %>% inner_join(metatumor %>% select(cell,Sample), by =  "cell")
aggregated_expression_Tumor = expression_df_Tumor %>% group_by(Sample) %>% select(-cell) %>% summarise_all(mean)
aggregated_expression_df_Tumor = aggregated_expression_Tumor %>% select(-Sample) %>% t() %>% magrittr::set_colnames(aggregated_expression_Tumor$Sample) %>% data.frame() %>% rownames_to_column("ligand") %>% as_tibble() 
aggregated_expression_matrix_Tumor = aggregated_expression_df_Tumor %>% select(-ligand) %>% as.matrix() %>% magrittr::set_rownames(aggregated_expression_df_Tumor$ligand)
#vis_ligand_expression = aggregated_expression_matrix_Tumor[order_ligands,]
vis_ligand_expression = aggregated_expression_matrix_Tumor[intersect(order_ligands,rownames(aggregated_expression_matrix_Tumor)),]
library(RColorBrewer)
color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
p_ligand_expression = vis_ligand_expression %>% make_heatmap_ggplot("Prioritized Tumor-ligands","Tumor", color = color[100],legend_position = "top", x_axis_position = "top", legend_title = "Expression\n(averaged over\nsingle cells)") + theme(axis.text.y = element_text(face = "italic"))
ggsave(paste0(cl,'_ligand_expression_heat.pdf'),p_ligand_expression,width = 7,height = 5)

#靶基因表达矩阵
expression_tcell <- as.matrix(objtcell@assays$RNA@data)
metatcell <- objtcell@meta.data
colnames(metatcell)[2] <- 'cell'
expression_df_target = t(expression_tcell)[,geneset_oi] %>% data.frame() %>% rownames_to_column("cell") %>% as_tibble() %>% inner_join(metatcell %>% select(cell,Sample), by =  "cell") 
aggregated_expression_target = expression_df_target %>% group_by(Sample) %>% select(-cell) %>% summarise_all(mean)
aggregated_expression_df_target = aggregated_expression_target %>% select(-Sample) %>% t() %>% magrittr::set_colnames(aggregated_expression_target$Sample) %>% data.frame() %>% rownames_to_column("target") %>% as_tibble() 
aggregated_expression_matrix_target = aggregated_expression_df_target %>% select(-target) %>% as.matrix() %>% magrittr::set_rownames(aggregated_expression_df_target$target)
vis_target_expression_scaled = aggregated_expression_matrix_target %>% t() %>% scale_quantile()
vis_target_expression_scaled <- vis_target_expression_scaled[,intersect(order_targets,colnames(vis_target_expression_scaled))]
p_target_scaled_expression = vis_target_expression_scaled  %>% make_threecolor_heatmap_ggplot("Sample","Target", low_color = color[1],mid_color = color[50], mid = 0.5, high_color = color[100], legend_position = "top", x_axis_position = "top" , legend_title = "Scaled expression\n(averaged over\nsingle cells)") + theme(axis.text.x = element_text(face = "italic"))
ggsave(paste0(cl,'_targetgene_tcell_scaled_expression_heat.pdf'),p_target_scaled_expression,width = 7,height = 5)


##将不同的热图组合在一个概览图中
figures_without_legend = plot_grid(
  p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""), 
  NULL,
  NULL,
  p_target_scaled_expression + theme(legend.position = "none", axis.ticks = element_blank()) + xlab(""), 
  align = "hv",
  nrow = 2,
  rel_widths = c(ncol(vis_ligand_pearson)+ 8, ncol(vis_ligand_expression)+8, ncol(vis_ligand_target)) -2,
  rel_heights = c(nrow(vis_ligand_pearson)+20, nrow(vis_target_expression_scaled) + 30))
legends = plot_grid(
  as_ggplot(get_legend(p_ligand_pearson)),
  as_ggplot(get_legend(p_ligand_target_network)),
  nrow = 2,
  align = "h")

p <- plot_grid(figures_without_legend, 
               legends, 
               rel_heights = c(20,6), nrow = 2, align = "hv")
ggsave(paste0(cl,'_grid.pdf'),plot = p,width = 10,height = 8)
ggsave(paste0(cl,'_grid.pdf'),plot = p,width = 10,height = 8)
save(objtumor,objtcell,
     expressed_genes_sender,expressed_genes_receiver,geneset_oi,
     expressed_ligands,expressed_receptors,lr_network_expressed,potential_ligands,
     ligand_activities,best_upstream_ligands,
     active_ligand_target_links_df,active_ligand_target_links,
     order_ligands,order_receptors,
     lr_network_top,best_upstream_receptors,lr_network_top_df,dist_receptors,hclust_receptors,order_receptors,dist_ligands,hclust_ligands,order_ligands_receptor,
     ligand_pearson_matrix,expression_df_Tumor,vis_ligand_expression,
     expression_df_target,vis_target_expression_scaled,
     p_ligand_pearson,p_ligand_expression,p_ligand_target_network,p_target_scaled_expression,p_ligand_receptor_network,
     file = paste0(cl,"_nichenet.RData"))
write.table(data.frame(rev(order_ligands)),paste0(cl,"_ligands.txt"),sep = "\t",quote = F,row.names = F,col.names = F)
