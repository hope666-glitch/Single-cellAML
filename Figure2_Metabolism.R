cytotoxic_genes <- c("KLRF1", "GNLY", "CTSW", "NKG7", "KLRD1", "GZMA", "ADGRG1", 
                     "CST7", "KLRK1", "FASLG", "HCST", "KLRB1", "ITGB1", "GZMB", "PRF1")

exhaustion_genes <- c("PDCD1","LAG3","HAVCR2","CD244","CD160",
                      "EOMES","NR4A2","PTGER4","TOX",
                      "TOX2","TIGIT","CTLA4","ENTPD1")

seurat_obj <- AddModuleScore(
  seurat_obj,
  features = list(cytotoxic_genes),
  name = "Cytotoxicity"
)

seurat_obj <- AddModuleScore(
  seurat_obj,
  features = list(exhaustion_genes),
  name = "Exhaustion"
)

meta_data <- seurat_obj@meta.data %>%
  dplyr::select(
    Cytotoxicity = Cytotoxicity1,  
    Exhaustion = Exhaustion1,
    group,
    celltype,EC
  )

plot_list <- list()
cytotoxic_genes <- c("KLRF1", "GNLY", "CTSW", "NKG7", "KLRD1", "GZMA", "ADGRG1", 
                     "CST7", "KLRK1", "FASLG", "HCST", "KLRB1", "ITGB1", "GZMB", "PRF1")

exhaustion_genes <- c("HAVCR2", "LAG3", "TIGIT", "CTLA4", "PDCD1", "LAYN", "BTLA", 
                      "TOX", "CXCL13")
cytotoxic_genes <-c("GZMA", "GZMB", "GZMH", "GZMK", "GZMH", "NKG7", "PRF1", "KLRF1", "KLRK1", "FCGR3A")
exhaustion_genes <- c("PDCD1", "LAYN", "HAVCR2", "LAG3", "CTLA4", "TOX", "VSIR", "BTLA", "ENTPD1", "CD160")
exhaustion_genes <- c("PDCD1","LAG3","HAVCR2","CD244","CD160",
                      "EOMES","NR4A2","PTGER4","TOX",
                      "TOX2","TIGIT","CTLA4","ENTPD1")
lactate_genes<- c("LDHA","LDHB","MCT1","MCT4","PDHA1","PDHB","PDK1")
glutamine_genes <- c("ALDH18A1","GLS","GLS2","GLUD1","GLUD2","GLUL","GOT2","KYAT1","OAT","PYCR1","PYCR2","PYCR3","RIMKLA","RIMKLB")
a <- read.csv("mit.csv")
b <- read.csv("magene.csv")
genes <- a$Mitochondrial.Genes..n.2030.
genes <- b$Gene.symbol
seurat_obj <- AddModuleScore(
  seurat_obj,
  features = list(cytotoxic_genes),
  name = "Cytotoxicity"
)

seurat_obj <- AddModuleScore(
  seurat_obj,
  features = list(exhaustion_genes),
  name = "Exhaustion"
)
Mitochondria
seurat_obj <- AddModuleScore(
  seurat_obj,
  features = list(genes),
  name = "Mitochondria"
)
seurat_obj <- AddModuleScore(
  seurat_obj,
  features = list(genes),
  name = "Mitochondria"
)
seurat_obj <- AddModuleScore(
  seurat_obj,
  features = list(metabolic_genes),
  name = "Mitometa"
)
seurat_obj <- AddModuleScore(
  seurat_obj,
  features = list(lactate_genes),
  name = "Lactate"
)
meta_data <- seurat_obj@meta.data %>%
  dplyr::select(
    Cytotoxicity = Cytotoxicity1,  
    Exhaustion = Exhaustion1,
    group,
    celltype,EC
  )
meta_data <- seurat_obj@meta.data %>%
  dplyr::select(
    Mitochondria=Mitochondria1,
    group,
    celltype,EC
  )

plot_list <- list()
for (module in c("Cytotoxicity", "Exhaustion")) {
  p <- ggplot(meta_data, aes(x = group, y = .data[[module]], fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +  
    geom_boxplot(width = 0.2, fill = "white") +  
    scale_fill_manual(values = c("CR" = "steelblue", "NR" = "red")) +
    labs(y = module, x = "Group") +
    theme_classic() +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.signif",  
      hide.ns = T       
    )
  plot_list[[module]] <- p
}
library(patchwork)
plotA <- plot_list[["Cytotoxicity"]] + plot_list[["Exhaustion"]] +
  plot_layout(guides = "collect")  
print(plotA)
library(ggplot2)
library(ggpubr)
library(patchwork)

plot_list_c <- list()
for (module in c("Cytotoxicity", "Exhaustion")) {
  means <- tapply(meta_data[[module]], meta_data$celltype, mean, na.rm=TRUE)
  celltype_order <- names(sort(means, decreasing = TRUE))
  meta_data$celltype <- factor(meta_data$celltype, levels = celltype_order)
  
  p <- ggplot(meta_data, aes(x = celltype, y = .data[[module]], fill = celltype)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.2, fill = "white") +
    labs(y = module, x = "Cell Type") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    stat_compare_means(
      method = "kruskal.test",
      label = "p.format",  
      hide.ns = TRUE
    )
  
  plot_list_c[[module]] <- p
}

plotC <- plot_list_c[["Cytotoxicity"]] + plot_list_c[["Exhaustion"]] +
  plot_layout(guides = "collect")

print(plotC)
deg <- FindMarkers(
  seurat_obj,
  ident.1 = "NFKBIZ+γδT",        
  ident.2 = "IFI44L+γδT",
  group.by = "celltype",
  min.pct = 0.1,          
  logfc.threshold = 0.25,
  test.use = "wilcox",        
  verbose = FALSE
)
ggplot(deg_combined, aes(x = avg_log2FC, y = log10_padj, color = significance)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c(
    "NFKBIZ+γδT up (sig)" = "#E53935",    
    "IFI44L+γδT up (sig)" = "#1E88E5",    
    "Not significant" = "#9E9E9E"         
  )) +
  geom_text_repel(
    aes(label = highlight),
    na.rm = TRUE,
    color = "black",
    size = 3.2,
    box.padding = 0.6,
    max.overlaps = 30,  
    segment.color = "gray50",
    segment.size = 0.3
  ) +
  geom_vline(xintercept = fc1.5_log2, linetype = "dashed", color = "gray30", linewidth = 0.8) +
  geom_vline(xintercept = -fc1.5_log2, linetype = "dashed", color = "gray30", linewidth = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray30", linewidth = 0.8) +
  xlim(-7, 7) +
  ylim(0, min(200, max(deg_combined$log10_padj, na.rm = TRUE) + 10)) +
  labs(
    x = "Log2FoldChange",
    y = "-Log10(P.adjust)",
    title = "DEG",
    color = "Group"  
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  annotate(
    "text",
    x = c(-3.5, 3.5),  
    y = max(deg_combined$log10_padj, na.rm = TRUE) + 5,
    label = c(
      paste0("IFI44L+γδT up (", sum(deg_combined$significance == "IFI44L+γδT up (sig)"), ")"),
      paste0("NFKBIZ+γδT up (", sum(deg_combined$significance == "NFKBIZ+γδT up (sig)"), ")")
    ),
    color = c("#1E88E5", "#E53935"),
    size = 4.5,
    fontface = "bold"
  ) 

ggplot(patient_lr_sig, aes(x = level, y = pathway)) +
  geom_point(
    data = subset(patient_lr_results, p.value < 0.05),  
    aes(size = neg_log10_p, color = direction)
  ) +
  scale_color_manual(values = c("Positive" = "#CC3333", "Negative" = "#0099CC")) +
  scale_size_continuous(range = c(2, 8)) +
  theme_classic() +
  labs(x = "", y = "", size = "-log10(P)", color = "NR/CR") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )
df_plot <- patient_lr_sig %>%
  mutate(
    sig = case_when(
      p.value < 0.05 ~ "p<0.05",
      TRUE ~ "NS"
    ),
    fill_group = paste0(direction, ".", sig),
    log10_OR = log10(OR),
    pathway = reorder(pathway, log10_OR)
  )

p_heat <- ggplot(df_plot, aes(x = 1, y = pathway)) +
  geom_tile(aes(fill = fill_group), color = "white", linewidth = 0.4) +  
  scale_fill_manual(
    values = c(
      "Positive.p<0.01" = "#800026",
      "Positive.p<0.05" = "#FF6B6B",
      "Positive.NS" = "white",
      "Negative.p<0.01" = "#2166ac",
      "Negative.p<0.05" = "#4ECDC4",
      "Negative.NS" = "white"
    ),
    breaks = c("Positive.p<0.01", "Positive.p<0.05", "Negative.p<0.01", "Negative.p<0.05"),
    labels = c("OR>1 P<0.01", "OR>1", "OR<1 P<0.01", "OR<1")
  ) +
  
  theme_classic() +
  labs(x = "", y = "", fill = "") +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 14, face = "bold", color = "black"),
    legend.key.size = unit(0.8, "cm"),  
    legend.position = "none",
    panel.border = element_blank(),
    axis.line.y = element_line(linewidth = 1.2, color = "black"),
    axis.line.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black")
  )

p_heat
ggsave("metaforheat.pdf",
       p_heat, width = 10, height = 5)