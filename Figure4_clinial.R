library(dplyr)
library(ggplot2)
library(pROC)
g_vars <- c("G1","G2","G3","G4")                  
meta_vars <- c(
  "Lysine degradation",
  "D-Glutamine and D-glutamate metabolism"
)                                                  
cell_vars <- c(
  "NFKBIZ+γδT",
  "IFI44L+γδT",
  "GZMK+γδT",
  "FGFBP2+γδT",
  "MYOM2+γδT"
)                                                  
outcome <- "Y"                                    
add_backtick <- function(x) paste0("`", x, "`")

g_vars_bt <- add_backtick(g_vars)
meta_vars_bt <- add_backtick(meta_vars)
cell_vars_bt <- add_backtick(cell_vars)
all_single <- c(g_vars_bt, meta_vars_bt, cell_vars_bt)

comb_g_meta    <- apply(expand.grid(g_vars_bt, meta_vars_bt), 1, paste, collapse = " + ")
comb_g_cell    <- apply(expand.grid(g_vars_bt, cell_vars_bt), 1, paste, collapse = " + ")
comb_meta_cell <- apply(expand.grid(meta_vars_bt, cell_vars_bt), 1, paste, collapse = " + ")
all_double <- c(comb_g_meta, comb_g_cell, comb_meta_cell)

all_triple <- apply(expand.grid(g_vars_bt, meta_vars_bt, cell_vars_bt), 1, paste, collapse = " + ")

all_formulas <- c(all_single, all_double, all_triple)
length(all_formulas) 

result <- data.frame()

for (f in all_formulas) {
  form <- paste0(outcome, " ~ ", f)
  fit <- glm(as.formula(form), data = df, family = binomial)
  
  aic  <- AIC(fit)
  pred <- predict(fit, type = "response")
  auc  <- as.numeric(auc(df[[outcome]], pred))
  
  result <- rbind(result, data.frame(
    formula = f,
    AIC = aic,
    AUC = auc
  ))
}

result <- result %>% arrange(AIC)
head(result, 10)

result$model_id <- 1:nrow(result)
best_aic <- result$AIC[5]

ggplot(result, aes(x = model_id, y = AIC)) +
  geom_point(size = 3, color = "black", shape = 16) +
  geom_hline(yintercept = best_aic, color = "red", linewidth = 1.2, linetype = 1) +
  labs(
    title = "",
    x = "Best models",
    y = "AIC"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.line = element_line(linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2.5)
  )
p <- ggplot(result, aes(x = model_id, y = AIC)) +
  geom_point(size = 2.5, color = "#2E86AB", shape = 16) +  
  geom_hline(yintercept = best_aic, color = "red", linewidth = 1.2, linetype = 1) +
  labs(
    title = "",
    x = "Best models",
    y = "AIC"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.line = element_line(linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2.5),
    panel.background = element_rect(fill = "#FBF9E4", color = NA)
  )
p
ggsave(
  filename = "lgAIC2.pdf",
  plot = p,
  width = 4,
  height = 3
)
ggplot(result, aes(x = model_id, y = AUC)) +
  geom_point(size = 3, color = "black", shape = 16) +
  labs(
    title = "",
    x = "Best models",
    y = "AUC"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.line = element_line(linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2.5)
  )
p1 <- ggplot(result, aes(x = model_id, y = AUC)) +
  geom_point(aes(color = AUC), size = 3, shape = 16) +
  scale_color_gradientn(
    colors = c('#4ECDC4','#83E377','#FFD166', '#FF9E6B', '#FF6B6B'),
    name = "AUC"
  ) +
  
  labs(
    title = "",
    x = "Best models",
    y = "AUC"
  ) +
  
  theme_classic() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 13, face = "bold"),
    axis.line = element_line(linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    legend.title = element_text(size = 14, face = "bold"), 
    legend.text = element_text(size = 12, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2.5),
    panel.background = element_rect(fill = "#FBF9E4", colour = NA)
  )
p1
ggsave(
  filename = "lgAIC1.pdf",
  plot = p,
  width = 4,
  height = 2.5
)
library(ComplexHeatmap)
library(grid)
clin_matrix <- as.matrix(clin_data)
p <- Heatmap(
  t(clin_matrix),
  col = cols_vec,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 14),
  rect_gp = gpar(col = "white", lwd = 1),
  show_column_names = FALSE,
  show_heatmap_legend = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE
)
draw(p)

p <- Heatmap(
  t(clin_matrix),
  col = cols_vec,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10, fontface = "bold"), 
  rect_gp = gpar(col = "white", lwd = 2),
  show_column_names = FALSE,
  show_heatmap_legend = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 14, fontface = "bold"),
    labels_gp = gpar(fontsize = 10, fontface = "bold")
  )
)

draw(p, padding = unit(c(2,2,2,2), "mm"))
k = 1
lgd = list()
for(i in 1:ncol(clin_data)){
  un = sort(unique(clin_data[,i]))
  ti = colnames(clin_data)[i]
  lgd[[i]] = Legend(
    labels = un,
    title = ti,
    legend_gp = gpar(fill = cols_vec[un]),
    title_gp = gpar(fontsize = 12, fontface = "bold"),  
    labels_gp = gpar(fontsize = 10, fontface = "bold")  
  )
}

col1 <- packLegend(
  list = lgd[c(1,2)], 
  direction = "vertical", 
  row_gap = unit(0.3, "cm")
)

col2 <- packLegend(
  list = lgd[c(3,4)], 
  direction = "vertical", 
  row_gap = unit(0.3, "cm")
)
col3 <- lgd[[5]]

final_legend <- packLegend(
  list = list(col1, col2, col3),  
  direction = "horizontal",    
  column_gap = unit(10, "mm")  
)

draw(pd)  

pdf("clinical_heatmap3.pdf", height = 1, width = 12)
draw(p, padding = unit(c(2,2,2,2), "mm"))
dev.off()