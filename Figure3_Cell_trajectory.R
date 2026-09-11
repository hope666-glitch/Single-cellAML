library(monocle)
library(RColorBrewer)
n_types <- length(unique(pData(cds)$celltype))
color_palette <- brewer.pal(n_types, "Set1")  
shape_palette <- 15:(15 + n_types - 1)  

p <- plot_cell_trajectory(
  cds,
  color_by = "celltype",
  label_branch_points = TRUE,
  cell_size = 1.5,
  palette = brewer.pal(length(unique(pData(cds)$celltype)), "Set1")
) +
  facet_wrap(~group, ncol = 2) +
  ggtitle("") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 12),  
    legend.position = "none", 
    axis.line = element_line(color = "black", linewidth = 1.5),
    axis.ticks = element_line(color = "black", linewidth = 1.2),
    axis.text = element_text(color = "black", size = 10,face ="bold" ),  
    axis.title = element_text(face = "bold", size = 11,color = "black"),  
    panel.border = element_rect(color = "black", linewidth = 1.2)
  )

print(p)
library(ggplot2)
library(RColorBrewer)
library(monocle3)
cell_types <- cell_type_order
fixed_palette <- celltype_strip_cols
names(fixed_palette) <- cell_types

p <- plot_cell_trajectory(
  cds,
  color_by = "celltype",
  label_branch_points = TRUE,
  cell_size = 2.2,       
  show_branch_points = F
) +
  scale_color_manual(values = fixed_palette) +  
  facet_wrap(~group, ncol = 2) +
  
  ggtitle("") +
  
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    strip.text = element_blank(),    
    strip.background = element_blank(), 
    legend.position = "none",
    axis.text = element_text(color = "black", size = 13, face = "bold"),
    axis.title = element_text(face = "bold", size = 16, color = "black"),
    axis.line = element_line(color = "black", linewidth = 1.2),
    axis.ticks = element_line(color = "black", linewidth = 1.5),
    panel.border = element_rect(color = "black", linewidth = 2, fill = NA),
    panel.grid = element_blank(),  
    plot.margin = margin(10,10,10,10)
  )

print(p)
p <- plot_cell_trajectory(
  cds,
  color_by = "celltype",
  label_branch_points = TRUE,  
  cell_size = 1.5,
  palette = brewer.pal(length(unique(pData(cds)$celltype)), "Set1")  
) +
  facet_wrap(~celltype, ncol = 2)+
  ggtitle("")+
  theme(
    strip.text = element_text(face = "bold", size = 12),  
    legend.position = "none", 
    panel.background = element_blank(),  
    axis.line = element_line(color = "black", linewidth = 1.5),
    axis.ticks = element_line(color = "black", linewidth = 1.2),
    axis.text = element_text(color = "black", size = 10,face ="bold" ),  
    axis.title = element_text(face = "bold", size = 11,color = "black")  
  )

print(p)
p <- plot_cell_trajectory(
  cds,
  color_by = "celltype",          
  label_branch_points = TRUE,     
  cell_size = 1.5,
  palette = brewer.pal(length(unique(pData(cds)$EC)), "Set1")  
) +
  facet_wrap(~celltype, ncol = 2)+ 
  ggtitle("")+
  theme(
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "none",
    panel.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1.5),
    axis.ticks = element_line(color = "black", linewidth = 1.2),
    axis.text = element_text(color = "black", size = 10, face = "bold"),  
    axis.title = element_text(face = "bold", size = 11, color = "black")  
  )

p_CR <- p %+% dplyr::filter(p$data, group == "CR") + 
  ggtitle("CR")

p_NR <- p %+% dplyr::filter(p$data, group == "NR") + 
  ggtitle("NR")
library(patchwork)
p_final <- p_CR / p_NR

print(p_final)
library(dplyr)
library(tibble)
pseudotime_deg <- differentialGeneTest(
  cds,
  fullModelFormulaStr = "~ Pseudotime",  
  reducedModelFormulaStr = "~ 1",        
  cores = 8, 
  verbose = FALSE  
)

top1900_degs <- pseudotime_deg %>%
  rownames_to_column("gene") %>%
  filter(qval < 0.05) %>%  
  arrange(qval) %>%
  head(1900) %>%
  pull(gene)  

cds <- setOrderingFilter(cds, ordering_genes = top1900_degs)

cds <- reduceDimension(cds, method = "DDRTree")
library(monocle)
library(ggplot2)
library(pheatmap)
expression_matrix <- seurat_obj[["RNA"]]$counts
cell_metadata <- seurat_obj@meta.data
gene_annotation <- data.frame(
  gene_short_name = rownames(seurat_obj),
  row.names = rownames(seurat_obj)
)

cds <- newCellDataSet(
  expression_matrix,
  phenoData = new("AnnotatedDataFrame", data = cell_metadata),
  featureData = new("AnnotatedDataFrame", data = gene_annotation),
  expressionFamily = negbinomial.size()
)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr=1) 
cds <- cds[fData(cds)$num_cells_expressed > 10, ]

disp_table <- dispersionTable(cds)
high_var_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1)
rownames(high_var_genes) <- high_var_genes$gene_id
top_var_genes <- row.names(high_var_genes)[order(high_var_genes$dispersion_empirical, decreasing = TRUE)]

cds <- setOrderingFilter(cds, top_var_genes)

cds <- reduceDimension(cds, method = "DDRTree")
cds <- orderCells(cds)
saveRDS(cds,"cds.rds")
cds <- orderCells(cds,root_state = )
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds)
plot_pc_variance_explained(cds, return_all = F)
cds <- reduceDimension(cds, max_components = 2, num_dim = 6,
                       reduction_method = 'tSNE', verbose = T)
cds <- clusterCells(cds, num_clusters = 6)
plot_cell_clusters(cds, 1, 2 )
table(pData(cds)$Cluster) 
colnames(pData(cds))
save(cds,file = 'input_cds.Rdata')
pData(cds)$Cluster <- pData(cds)$celltype  
table(pData(cds)$Cluster)
diff_test_res <- differentialGeneTest(
  cds,
  fullModelFormulaStr = "~Cluster",
  cores = 1,
  verbose = FALSE
)
sig_genes <- subset(diff_test_res, qval < 0.1)
sig_genes <- sig_genes[order(sig_genes$pval), ]
cg <- as.character(head(sig_genes$gene_short_name, 6))

plot_genes_jitter(cds[cg,], grouping = "Cluster", color_by = "Cluster", nrow = 3)

high_var_genes <- high_var_genes[order(high_var_genes$dispersion_empirical, decreasing = TRUE), ]
top_var_genes <- head(high_var_genes$gene_id, 2000)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
cds <- setOrderingFilter(cds, ordering_genes)