#####
# Do initial preprocessing and save RDS
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
date="20230511"

cmap <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]$`Tableau 20`$value
cmap <- c(cmap, "firebrick3", "navyblue")

##### Import and merge data
dat <- Read10X(data.dir='BD1/merged_wBroadRun230120/filtered_feature_bc_matrix/') #data from this and prior sequencing of this sample
seurat_table_1 <- CreateSeuratObject(dat$`Gene Expression`, min.cells = 3, project="HT2", min.features=0) 
dat <- Read10X(data.dir='BD2/merged_wBroadRun230120/filtered_feature_bc_matrix/') #data from this and prior sequencing of this sample
seurat_table_2 <- CreateSeuratObject(dat$`Gene Expression`, min.cells = 3, project="HT3", min.features=0)
dat <- Read10X(data.dir='BD3/filtered_feature_bc_matrix/')
seurat_table_3 <- CreateSeuratObject(dat$`Gene Expression`, min.cells = 3, project="HT4", min.features=0)
dat <- Read10X(data.dir='BD4/filtered_feature_bc_matrix/')
seurat_table_4 <- CreateSeuratObject(dat$`Gene Expression`, min.cells = 3, project="HT5", min.features=0)
dat <- Read10X(data.dir='BD5/filtered_feature_bc_matrix/')
seurat_table_5 <- CreateSeuratObject(dat$`Gene Expression`, min.cells = 3, project="HT6", min.features=0)

seurat_table <- merge( seurat_table_1, y = c( seurat_table_2,seurat_table_3,seurat_table_4,seurat_table_5 ) )

#print some QC general info
for (obj in c(seurat_table_1, seurat_table_2,seurat_table_3,seurat_table_4,seurat_table_5)){
  print(c(ncol(x=obj),nrow(x=obj)))
  print(c(median(obj$nCount_RNA),mean(obj$nCount_RNA),median(obj$nFeature_RNA),mean(obj$nFeature_RNA)))
}

###### filter data
seurat_table[["percent.mt"]] <- PercentageFeatureSet(seurat_table, pattern = "^MT-")

# QC figure
x <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=nFeature_RNA ) ) +
  geom_hline( yintercept=1000, color="red" ) + #line where cutoff is
  geom_hline( yintercept=7000, color="red" ) + #line where cutoff is
  geom_violin( fill="steelblue" ) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2, outlier.shape=NA) +
  theme_bw() +
  xlab("")

y <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=nCount_RNA ) ) +
  geom_violin( fill="pink" ) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2, outlier.shape=NA) +
  theme_bw() +
  xlab("")

z <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=percent.mt ) ) +
  geom_hline( yintercept=12.5, color="red" ) +
  geom_violin( fill="orchid4" ) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2, outlier.shape=NA) +
  theme_bw() +
  xlab("")

pdf( paste0("figures/qc_violin_prefilter_",date,".pdf" ), height=3, width=8 )
plot_grid( x,y,z, ncol=3 )
dev.off()
##### 
seurat_table <- subset(seurat_table, subset = percent.mt < 12.5)
seurat_table <- subset(seurat_table, subset = nFeature_RNA > 1000) 
seurat_table <- subset(seurat_table, subset = nFeature_RNA < 7000) 
##### 
x <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=nFeature_RNA ) ) +
  geom_hline( yintercept=1000, color="red" ) +
  geom_hline( yintercept=7000, color="red" ) +
  geom_violin( fill="steelblue" ) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2, outlier.shape=NA) +
  theme_bw() +
  xlab("")

y <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=nCount_RNA ) ) +
  geom_violin( fill="pink" ) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2, outlier.shape=NA) +
  theme_bw() +
  xlab("")

z <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=percent.mt ) ) +
  geom_hline( yintercept=12.5, color="red" ) +
  geom_violin( fill="orchid4" ) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2, outlier.shape=NA) +
  theme_bw() +
  xlab("")

pdf( paste0("figures/qc_violin_postfilter_",date,".pdf" ), height=3, width=8 )
plot_grid( x,y,z, ncol=3 )
dev.off()

#####
#normalize data
seurat_table <- NormalizeData(seurat_table, normalization.method = "LogNormalize", scale.factor = 10000)

nfeatures = 2000
seurat_table <- FindVariableFeatures(seurat_table, selection.method = "vst", nfeatures = nfeatures)

top10 <- head(VariableFeatures(seurat_table), 10)

plot1 <- VariableFeaturePlot(seurat_table)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf( paste0("figures/variable_feature_",date,".pdf" ), height=4, width=6 )
plot1
plot2
dev.off()


#####
#run dimensionality reduction
all.genes <- rownames(seurat_table)
seurat_table <- ScaleData(seurat_table, features = VariableFeatures(object = seurat_table))
seurat_table <- RunPCA(seurat_table, features = VariableFeatures(object = seurat_table))

sink( paste0("text_outputs/top_PC_genes_1-30_",date,".txt"), append=F )
print(seurat_table[["pca"]], dims = 1:30, nfeatures = 5)
sink()

pdf( paste0("figures/PC_dim_loading_",date,".pdf" ), height=40, width=12 )
VizDimLoadings(seurat_table, dims = 1:30, reduction = "pca")
dev.off()

pdf( paste0("figures/pca_",date,".pdf" ), height=8, width=12 )
DimPlot(seurat_table, reduction = "pca", group.by = "orig.ident")
dev.off()

pdf( paste0("figures/PC_dim_heatmap_",date,".pdf" ), height=40, width=12 )
DimHeatmap(seurat_table, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()

pdf( paste0("figures/PC_elbow_postfilter_",date,".pdf" ), height=5, width=6 )
ElbowPlot(seurat_table, ndims=50)
dev.off()

#####
# cluster and umap
seurat_table <- FindNeighbors(seurat_table, dims = 1:30)
seurat_table <- FindClusters(seurat_table, resolution = 2)
seurat_table <- RunUMAP(seurat_table, dims = 1:30, seed.use = 123)

pdf( paste0("figures/umap_prefiltering_",date,".pdf" ), height=8, width=12 )
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1)
dev.off()

pdf( paste0("figures/umap_features_prefiltering_", date,".pdf"), height=15, width=15 )
FeaturePlot(seurat_table, features=c("EPCAM","PDPN","ITGB4","CAV1","CD4",
                                     "CCL19","AIRE","HLA-DRA","KRT5",
                                     "MYOG","CKM","DLK1",
                                     "FOXI1","CFTR","CGA",
                                     "NEUROD1","NKX6-2","ATOH1","OTOF",
                                     "DRGX","NTRK3","IVL","POU2F3","IL18","NMU","IL2","IL21"),order=TRUE)
dev.off()


## Predict doublets: 
library(DoubletFinder)

doubletfinder_fn <- function(test,out_prefix){
  nExp <- round(ncol(test) * 0.04)  # expect X% doublets
  test <- doubletFinder_v3(test, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
  
  DF.name = colnames(test@meta.data)[grepl("DF.classification", colnames(test@meta.data))]
  pdf( paste0("figures/umap_doubletFinder_",out_prefix,'_',date,".pdf" ), height=15, width=25 )
  print(cowplot::plot_grid(ncol = 2, DimPlot(test, group.by = "orig.ident") + NoAxes(),
                           DimPlot(test, group.by = DF.name) + NoAxes()))
  dev.off()
  return(test)
}

#separate seurat object to use doublet finder
for (donor in c('HT2','HT3','HT4','HT5','HT6')){
  test = doubletfinder_fn(subset(seurat_table,subset=orig.ident==donor),donor)
  
  DF.name = colnames(test@meta.data)[grepl("DF.classification", colnames(test@meta.data))] #name of column in metadata where doubletfinder (DF) results stored
  test = as.data.frame(test@meta.data[DF.name]) #extract metadata, save as dataframe
  colnames(test) <- 'df'
  
  if(donor=='HT2'){ merged_metadata = test } #first one, save as variable
  if(donor!='HT2'){ merged_metadata = rbind(merged_metadata, test) } #after the first, concatenate variable
}

seurat_table <- AddMetaData(seurat_table, merged_metadata, col.name='Doublet.Finder')

pdf( paste0("figures/umap_doubletFinder_all_",date,".pdf" ), height=15, width=25 )
print(cowplot::plot_grid(ncol = 2, DimPlot(seurat_table, group.by = "orig.ident") + NoAxes(),
                         DimPlot(seurat_table, group.by = 'Doublet.Finder') + NoAxes()))
dev.off()

pdf( paste0("figures/umap_prefiltering_postDF_",date,".pdf" ), height=8, width=12 )
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1)
dev.off()

##### 
# Analyze cluster composition - look for 'contaminating' fibroblasts/T cells
pdf( paste0("figures/vln_contam-features_prefiltering_", date,".pdf"), height=10, width=20 )
VlnPlot(seurat_table, features = c("CAV1","CD4","LCK","CD1B","CD247","CD3E"), slot = "counts", log = TRUE)
dev.off()
pdf( paste0("figures/umap_contam-features_prefiltering_", date,".pdf"), height=15, width=15 )
FeaturePlot(seurat_table, features=c("CAV1","CD4","LCK","CD1B","CD247","CD3E"), pt.size=0.5, order=T)
dev.off()

markers_24 <- FindMarkers(seurat_table, ident.1=c(24), min.pct=0.1, logfc.threshold=1, only.pos=T, test.use="wilcox")
markers_35 <- FindMarkers(seurat_table, ident.1=c(35), min.pct=0.1, logfc.threshold=1, only.pos=T, test.use="wilcox")
markers_24_35 <- FindMarkers(seurat_table, ident.1=c(24,35), min.pct=0.1, logfc.threshold=1, only.pos=T, test.use="wilcox")
markers_26 <- FindMarkers(seurat_table, ident.1=c(26), min.pct=0.1, logfc.threshold=1, only.pos=T, test.use="wilcox")
markers_36 <- FindMarkers(seurat_table, ident.1=c(36), min.pct=0.1, logfc.threshold=1, only.pos=T, test.use="wilcox")
markers_26_36 <- FindMarkers(seurat_table, ident.1=c(26,36), min.pct=0.1, logfc.threshold=1, only.pos=T, test.use="wilcox")

#####
#remove contaminating cells
keep <- seurat_table@active.ident!=26 & seurat_table@active.ident!=36 #fib
seurat_table <- seurat_table[,keep]

keep <- seurat_table@active.ident!=24 & seurat_table@active.ident!=35 #t
seurat_table <- seurat_table[,keep]

DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1)

############### Repeat normalization, find var features, etc
#normalize data
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
seurat_table <- FindClusters(seurat_table, resolution = 2)
seurat_table <- RunUMAP(seurat_table, dims = 1:30, seed.use = 123)


pdf( paste0("figures/umap_",date,".pdf" ), height=8, width=12 )
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1)
dev.off()


pdf( paste0("figures/umap_doubletFinder_all_postfilter_",date,".pdf" ), height=15, width=25 )
print(cowplot::plot_grid(ncol = 2, DimPlot(seurat_table, group.by = "orig.ident") + NoAxes(),
                         DimPlot(seurat_table, group.by = 'Doublet.Finder') + NoAxes()))
dev.off()


pdf( paste0("figures/umap_split_", date,".pdf"), height=5, width=15 )
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1,split.by='orig.ident' ) + theme(legend.position="none")
dev.off()

pdf( paste0("figures/umap_group_", date,".pdf"), height=5, width=8 )
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=0.1,group.by="orig.ident" ) 
dev.off()

pdf( paste0("figures/umap_features_", date,".pdf"), height=15, width=15 )
FeaturePlot(seurat_table, features=c("EPCAM","PDPN","ITGB4","CAV1","CD4",
                                     "CCL19","AIRE","HLA-DRA","KRT5",
                                     "MYOG","CKM","DLK1",
                                     "FOXI1","CFTR","CGA",
                                     "NEUROD1","NKX6-2","ATOH1","OTOF",
                                     "DRGX","NTRK3","IVL","POU2F3","IL18","NMU","IL2","IL21"),order=TRUE)
dev.off()



#####
#remove cluster 34  - almost exclusively predicted doublets
keep <- seurat_table@active.ident!=34
seurat_table <- seurat_table[,keep]

############### Repeat normalization, find var features, etc
#normalize data
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
seurat_table <- FindClusters(seurat_table, resolution = 2)
seurat_table <- RunUMAP(seurat_table, dims = 1:30, seed.use = 123)


pdf( paste0("figures/umap_post34removed_",date,".pdf" ), height=8, width=12 )
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1)
dev.off()


pdf( paste0("figures/umap_doubletFinder_all_post34removed_",date,".pdf" ), height=15, width=25 )
print(cowplot::plot_grid(ncol = 2, DimPlot(seurat_table, group.by = "orig.ident") + NoAxes(),
                         DimPlot(seurat_table, group.by = 'Doublet.Finder') + NoAxes()))
dev.off()


pdf( paste0("figures/umap_split_post34removed_", date,".pdf"), height=5, width=15 )
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1,split.by='orig.ident' ) + theme(legend.position="none")
dev.off()

pdf( paste0("figures/umap_group_post34removed_", date,".pdf"), height=5, width=8 )
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=0.1,group.by="orig.ident" ) 
dev.off()

pdf( paste0("figures/umap_features_post34removed_", date,".pdf"), height=15, width=15 )
FeaturePlot(seurat_table, features=c("EPCAM","PDPN","ITGB4","CAV1","CD4",
                                     "CCL19","AIRE","HLA-DRA","KRT5",
                                     "MYOG","CKM","DLK1",
                                     "FOXI1","CFTR","CGA",
                                     "NEUROD1","NKX6-2","ATOH1","OTOF",
                                     "DRGX","NTRK3","IVL","POU2F3","IL18","NMU","IL2","IL21"),order=TRUE)
dev.off()


#####
#save work
saveRDS(seurat_table, file = "Rdata/meclo_ht2-ht6_seurat_table_20230511-clust34removed.rds")
