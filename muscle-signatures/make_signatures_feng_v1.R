# explore RDS from Feng et al (muscle dataset)
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

seurat_feng <- readRDS( "geo download/GSE193346_C57BL6_seurat_object.rds" )

date='20230630'

seurat_feng <- SetIdent(seurat_feng, value = seurat_feng@meta.data[["cell_type_spec"]])

pdf( paste0("figures/umap_feng_",date,".pdf" ), height=8, width=10 ) #height=10, width=12 )
DimPlot(seurat_feng, reduction="umap", label=T, repel=T, pt.size=1 ) #+ theme(legend.position="none")
dev.off()

markers <- FindAllMarkers(seurat_feng, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 1)
for (cluster in unique(markers$cluster)){
  print(cluster)
  glist = markers[markers$cluster==cluster,'gene']
  write_delim( as.data.frame(glist), paste0("text_outputs/signature_",gsub('/','_',cluster),"_logfc1_pct10_",date,".txt"), col_names=F, delim="\t" )
}

