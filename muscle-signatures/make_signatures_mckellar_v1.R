# explore RDS from McKellar et al (muscle dataset)
library(tidyverse)
library(Seurat)
library(hdf5r)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(cowplot)
library(ape)
library(scales)
library(ggthemes)
library(ggrepel)

temp <- load( "Dryad download/myo_slim_seurat_v1-1.RData" )

date='20230621'

# Phate1Bins
myo.slim.seurat <- SetIdent(myo.slim.seurat, value = myo.slim.seurat@meta.data[["phate1.bins"]])
pdf( paste0("figures/umap_phate1bins_factorIDs_",date,".pdf" ), height=8, width=10 ) #height=10, width=12 )
DimPlot(myo.slim.seurat, reduction="phate_harmony", label=T, repel=T, pt.size=1 ) #+ theme(legend.position="none")
dev.off()
markers <- FindAllMarkers(myo.slim.seurat, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 1)
for (cluster in unique(markers$cluster)){
  print(cluster)
  glist = markers[markers$cluster==cluster,'gene']
  write_delim( as.data.frame(glist), paste0("text_outputs/phate1bins",cluster,"_logfc1_pct10_",date,".txt"), col_names=F, delim="\t" )
}

# Harmony
myo.slim.seurat <- SetIdent(myo.slim.seurat, value = myo.slim.seurat@meta.data[["harmony_factorIDs"]])
pdf( paste0("figures/umap_harmony_factorIDs_",date,".pdf" ), height=8, width=10 ) #height=10, width=12 )
DimPlot(myo.slim.seurat, reduction="phate_harmony", label=T, repel=T, pt.size=1 ) #+ theme(legend.position="none")
dev.off()
markers <- FindAllMarkers(myo.slim.seurat, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 1)
for (cluster in unique(markers$cluster)){
  print(cluster)
  glist = markers[markers$cluster==cluster,'gene']
  write_delim( as.data.frame(glist), paste0("text_outputs/harmony",gsub('/','_',cluster),"_logfc1_pct10_",date,".txt"), col_names=F, delim="\t" )
}

#prefilter_factorIDS
myo.slim.seurat <- SetIdent(myo.slim.seurat, value = myo.slim.seurat@meta.data[["prefilter_factorIDs"]])
pdf( paste0("figures/umap_prefilter_factorIDs_",date,".pdf" ), height=8, width=10 ) #height=10, width=12 )
DimPlot(myo.slim.seurat, reduction="phate_harmony", label=T, repel=T, pt.size=1 ) #+ theme(legend.position="none")
dev.off()

#scanorama
myo.slim.seurat <- SetIdent(myo.slim.seurat, value = myo.slim.seurat@meta.data[["scanorama_factorIDs"]])
pdf( paste0("figures/umap_scanorama_factorIDs_",date,".pdf" ), height=8, width=10 ) #height=10, width=12 )
DimPlot(myo.slim.seurat, reduction="phate_harmony", label=T, repel=T, pt.size=1 ) #+ theme(legend.position="none")
dev.off()
markers <- FindAllMarkers(myo.slim.seurat, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 1)
for (cluster in unique(markers$cluster)){
  print(cluster)
  glist = markers[markers$cluster==cluster,'gene']
  write_delim( as.data.frame(glist), paste0("text_outputs/scanorama",cluster,"_logfc1_pct10_",date,".txt"), col_names=F, delim="\t" )
}

#bbknn
myo.slim.seurat <- SetIdent(myo.slim.seurat, value = myo.slim.seurat@meta.data[["bbknn_factorIDs"]])
pdf( paste0("figures/umap_bbknn_factorIDs_",date,".pdf" ), height=8, width=10 ) #height=10, width=12 )
DimPlot(myo.slim.seurat, reduction="phate_harmony", label=T, repel=T, pt.size=1 ) #+ theme(legend.position="none")
dev.off()
markers <- FindAllMarkers(myo.slim.seurat, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 1)
for (cluster in unique(markers$cluster)){
  print(cluster)
  glist = markers[markers$cluster==cluster,'gene']
  write_delim( as.data.frame(glist), paste0("text_outputs/bbknn",cluster,"_logfc1_pct10_",date,".txt"), col_names=F, delim="\t" )
}
