#####
# Visualize predictions from scGPT

#####
library(tidyverse)
library(Seurat)
library(Matrix)
library(hdf5r)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(cowplot)
library(ape)
library(scales)
library(ggthemes)
library(future)
library(tidyverse)

plan("sequential")

date="20240314"

cmap_hu_for_ms <- c("#59A14F","#B6992D","#d3d3d3","#E15759","#FFBE7D","#F28E2B","#F1CE63","#D4A6C8","#FF9D9A",
                    "#000000","#D37295","#8CD17D","#499894","#86BCB6","#BAB0AC","#FABFD2","#B07AA1","#9D7660")

cmap_hu_for_zf <- c("#59A14F","#FABFD2","#BAB0AC","#d3d3d3","#E15759","#FFBE7D","#F28E2B","#F1CE63","#D4A6C8",
                    "#FF9D9A","#000000","#499894","#D37295","#8CD17D","#B6992D","#86BCB6","#B07AA1","#9D7660")

#####
#read in scGPT predictions
pred_df_zf <- read.csv("/zebrafish_obs_info_mimetics-only_version.csv",header = T, row.names = 1, quote = "")
pred_df_ms <- read.csv("/mouse_obs_info_mimetics-only_version.csv",header = T, row.names = 1, quote = "")

human_compartments = c("Transitional","Keratinocyte","HMX2+ ionocyte","ASCL3+ ionocyte",
                       "Tuft","Early muscle","Late muscle","lncRNA-enriched muscle","CUX2+ neuro",     
                       "FEZF2+ neuro","NKX6-2+ neuro","SHOX2+ neuro","DRGX+ sensory neuro","Cochlear hair",'Novel')

########################

combine_predictions_for_all_cells <- function(pred_df, orig){
  all_cells = orig@assays[["RNA"]]@counts@Dimnames[[2]] #list of all cells in original seurat object
  combined_pred_df = data.frame( matrix(ncol=1, nrow=length(all_cells)) )
  row.names(combined_pred_df) = all_cells
  colnames(combined_pred_df)='scgpt_prediction'
  combined_pred_df['scgpt_prediction'] = 'excluded' #set all to 'excluded' then overwrite the cells with their prediction
  
  y = pred_df['predictions'] 
  y['cells'] = row.names(pred_df)
  combined_pred_df['cells'] = row.names(combined_pred_df)
  colnames(y)=colnames(combined_pred_df)
  combined_pred_df2 = rows_update(combined_pred_df, y, by='cells')
  return(combined_pred_df2)
}

process_fn <- function(seurat_table,nfeatures = 2000, res = 0.8){
  seurat_table <- NormalizeData(seurat_table, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_table <- FindVariableFeatures(seurat_table, selection.method = "vst", nfeatures = nfeatures)
  top10 <- head(VariableFeatures(seurat_table), 10)
  
  #####
  #run dimensionality reduction
  all.genes <- rownames(seurat_table)
  seurat_table <- ScaleData(seurat_table, features = VariableFeatures(object = seurat_table))
  seurat_table <- RunPCA(seurat_table, features = VariableFeatures(object = seurat_table))
  
  #####
  # cluster and umap
  seurat_table <- FindNeighbors(seurat_table, dims = 1:30)
  seurat_table <- FindClusters(seurat_table, resolution = res)
  seurat_table <- RunUMAP(seurat_table, dims = 1:30, seed.use = 123)
  return(seurat_table)
}

###########
zebra_seurat <- readRDS("zebrafish_seurat.rds")
mouse_seurat <- readRDS("mouse_seurat.rds")

################
# Repeat normalization and plot predictions on new UMAP of just mimetic cells
### Zebrafish
zebra_seurat@meta.data[["clusterlabel"]] = zebra_seurat@active.ident
keep <- zebra_seurat@active.ident!='aire-stage' & zebra_seurat@active.ident!='Immature' & zebra_seurat@active.ident!='cTEC' 
zebra_seurat <- zebra_seurat[,keep]
DimPlot(zebra_seurat, reduction="umap", label=T, repel=T, pt.size=1) + theme(legend.position="none")

zebra_seurat = process_fn(zebra_seurat)
zebra_seurat@active.ident = zebra_seurat@meta.data[["clusterlabel"]] 

DimPlot(zebra_seurat, reduction="umap", label=T, repel=T, pt.size=1 ) + theme(legend.position="none")

combined_pred_df2 = combine_predictions_for_all_cells(pred_df_zf, zebra_seurat)
zebra_seurat@meta.data$scgpt_prediction <- combined_pred_df2$scgpt_prediction

DimPlot(zebra_seurat, group.by="scgpt_prediction", cols=cmap_hu_for_zf, pt.size=1)#, label=T, repel=T, pt.size=1) + theme(legend.position="none")

for (i in human_compartments) {
  print(DimPlot(zebra_seurat, cells.highlight= list(zebra_seurat@assays[["RNA"]]@counts@Dimnames[[2]] [ zebra_seurat@meta.data[["scgpt_prediction"]]==i ]), label=F, repel=T, pt.size=1, cols.highlight = "#3F007D") +ggtitle(i)+ theme(legend.position="none")) #+ plot_annotation(title = i) 
}


### Mouse
mouse_seurat@meta.data[["clusterlabel"]] = mouse_seurat@active.ident
keep <- mouse_seurat@active.ident!='perinatal cTEC' & mouse_seurat@active.ident!='adult cTEC' &
  mouse_seurat@active.ident!='TA MEC' & mouse_seurat@active.ident!='Immature MEC' & mouse_seurat@active.ident!='Aire-stage'
mouse_seurat <- mouse_seurat[,keep]
DimPlot(mouse_seurat, reduction="umap", label=T, repel=T, pt.size=1) + theme(legend.position="none")

mouse_seurat = process_fn(mouse_seurat)
mouse_seurat@active.ident = mouse_seurat@meta.data[["clusterlabel"]]

DimPlot(mouse_seurat, reduction="umap", label=T, repel=T, pt.size=1 ) + theme(legend.position="none")

combined_pred_df2 = combine_predictions_for_all_cells(pred_df_ms, mouse_seurat)
mouse_seurat@meta.data$scgpt_prediction <- combined_pred_df2$scgpt_prediction

DimPlot(mouse_seurat, group.by="scgpt_prediction", cols=cmap_hu_for_ms, pt.size=0.5)#, label=T, repel=T, pt.size=1) + theme(legend.position="none")

for (i in human_compartments) {
  print(DimPlot(mouse_seurat, cells.highlight= list(mouse_seurat@assays[["RNA"]]@counts@Dimnames[[2]] [ mouse_seurat@meta.data[["scgpt_prediction"]]==i ]), label=F, repel=T, cols.highlight = "#3F007D") +ggtitle(i)+ theme(legend.position="none")) #+ plot_annotation(title = i) 
}

