rm(list = ls())
library(stringr)
library(Seurat)
library(dplyr)
seurat_obj <- PercentageFeatureSet(seurat_obj, 
                                   pattern = "^MT-", 
                                   col.name = "PMT")

seurat_obj <- PercentageFeatureSet(seurat_obj, 
                                   pattern = "^HBA|^HBB", 
                                   col.name = "PhB")
VlnPlot(seurat_obj,
        features = c("nCount_RNA","nFeature_RNA","PMT","PhB"),
        log = T,
        pt.size = 0,group.by = "orig.ident")

seurat_obj <- subset(seurat_obj, 
                     subset = nFeature_RNA > 200 &
                       nFeature_RNA < 5000 &
                       nCount_RNA > 200 & 
                       nCount_RNA < 30000 &
                       PMT < 20 &
                       PhB<5&
                       !is.na(PMT) &
                       !is.na(PhB))

seurat_obj <- NormalizeData(seurat_obj)
library(tidyverse)
library(patchwork)
seurat_obj <- FindVariableFeatures(seurat_obj)
scale.genes <- VariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
library(harmony)
seurat_obj <- RunHarmony(seurat_obj, group.by.vars ="orig.ident")
seurat_obj <- RunTSNE(seurat_obj, dims = 1:15,reduction = "harmony")
seurat_obj <- RunUMAP(seurat_obj, dims = 1:15,reduction = "harmony")
plot5 <- DimPlot(seurat_obj,reduction = "umap",label = T)
for (res in c(0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,1,2)) {
  seurat_obj <- FindClusters(seurat_obj, graph.name ="RNA_snn", resolution = res, algorithm =1)}
apply(seurat_obj@meta.data[,grep("RNA_snn_res",colnames(seurat_obj@meta.data))],2,table)

library(clustree)
p2_tree <- clustree(seurat_obj@meta.data, prefix ="RNA_snn_res.")
p2_tree

seurat_obj <- FindClusters(seurat_obj, resolution = 0.01)

seurat_obj.markers <- FindAllMarkers(seurat_obj, only.pos =TRUE, min.pct =0.25, logfc.threshold =0.25, assay ="RNA")

seurat_obj.markers <- seurat_obj.markers[seurat_obj.markers$p_val_adj<0.05,]
seurat_obj.markers %>%
  group_by(cluster)%>%
  top_n(n =100, wt = avg_log2FC)-> top150

top5seurat.markers <- seurat_obj.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)

