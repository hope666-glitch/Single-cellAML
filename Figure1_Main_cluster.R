library(ClusterGVis)

package_path <- "ClusterGVis_0.1.4.tar.gz"

Idents(seurat_obj) <- "seurat_clusters"
st.data <- prepareDataFromscRNA(object = seurat_obj,
                                diffData = top5seurat.markers,
                                showAverage= TRUE)

library(org.Hs.eg.db)
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 5,
                        seed = 5201314)
head(enrich)
library(dplyr)
top5_markers_per_cluster <- seurat_obj.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%  
  ungroup()

markGenes <- unique(top5_markers_per_cluster$gene)  
markGenes
stopifnot(all(unique(st.data$cluster) %in% c(0:4)))
stopifnot(all(unique(enrich$cluster) %in% c(0:4)))
visCluster(object = st.data,
           plot.type ="line")
visCluster(object = st.data,
           plot.type = "heatmap",
           column_names_rot = 45,
           markGenes = markGenes,
           cluster.order = c(0:4))

visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           cluster.order = c(1:5),
           textbar.pos = c(0.7,0.1),
           add.bar = T)
theme_set(theme(
  text = element_text(face = "bold"),
  axis.text = element_text(face = "bold"),
  legend.text = element_text(face = "bold")
))
cell_colors <- c("#6495ED","#DD3F4E","#FFBF00","#40B5AD","#F3B0C3")
pdf('sc_optimized_layout4.pdf', height = 6, width = 12, onefile = FALSE)

visCluster(
  object = st.data,
  plot.type = "both",           
  column_names_rot = 45,        
  show_row_dend = FALSE,        
  markGenes = markGenes,        
  markGenes.side = "left",      
  annoTerm.data = enrich,ctAnno.col =cell_colors,       
  line.side = "left",           
  cluster.order = c(1:5),       
  go.col = rep(ggsci::pal_d3()(5), each = 5),  
  textbar.pos = c(0.9, 0.1),  
  add.bar = TRUE
)

dev.off()