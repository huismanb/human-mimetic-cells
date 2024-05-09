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

date="20231017"

cmap <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]$`Tableau 20`$value
cmap <- c(cmap, "firebrick3", "navyblue")

#####
#read gene expression data
dat <- Read10X(data.dir='4wpfThymus_1/filtered_feature_bc_matrix')
seurat_table_1 <- CreateSeuratObject(dat, min.cells = 3, project="4wpfThymus_1", min.features=0) 
dat <- Read10X(data.dir='4wpfThymus_3/filtered_feature_bc_matrix')
seurat_table_2 <- CreateSeuratObject(dat, min.cells = 3, project="4wpfThymus_3", min.features=0)
dat <- Read10X(data.dir='4wpfThymus_5/filtered_feature_bc_matrix')
seurat_table_3 <- CreateSeuratObject(dat, min.cells = 3, project="4wpfThymus_5", min.features=0)
dat <- Read10X(data.dir='4wpfThymus_6/filtered_feature_bc_matrix')
seurat_table_4 <- CreateSeuratObject(dat, min.cells = 3, project="4wpfThymus_6", min.features=0)
dat <- Read10X(data.dir='adult_thy_liberase_1/filtered_feature_bc_matrix')
seurat_table_5 <- CreateSeuratObject(dat, min.cells = 3, project="adult_thy_liberase_1", min.features=0)
dat <- Read10X(data.dir='adult_thy_liberase_2/filtered_feature_bc_matrix')
seurat_table_6 <- CreateSeuratObject(dat, min.cells = 3, project="adult_thy_liberase_2", min.features=0)
dat <- Read10X(data.dir='adult_thy_untreated/filtered_feature_bc_matrix')
seurat_table_7 <- CreateSeuratObject(dat, min.cells = 3, project="adult_thy_untreated", min.features=0)
dat <- Read10X(data.dir='Fish_22_Thymus_/filtered_feature_bc_matrix')
seurat_table_8 <- CreateSeuratObject(dat, min.cells = 3, project="Fish_22_Thymus_", min.features=0)
dat <- Read10X(data.dir='Fish_24_Thymus_/filtered_feature_bc_matrix')
seurat_table_9 <- CreateSeuratObject(dat, min.cells = 3, project="Fish_24_Thymus_", min.features=0)
dat <- Read10X(data.dir='lck-eGFP_Thymus_1_/filtered_feature_bc_matrix')
seurat_table_10 <- CreateSeuratObject(dat, min.cells = 3, project="lck-eGFP_Thymus_1_", min.features=0)
dat <- Read10X(data.dir='lck-eGFP_Thymus_2_/filtered_feature_bc_matrix')
seurat_table_11 <- CreateSeuratObject(dat, min.cells = 3, project="lck-eGFP_Thymus_2_", min.features=0)

seurat_table <- merge( seurat_table_1, y = c( seurat_table_2,seurat_table_3,seurat_table_4,seurat_table_5,seurat_table_6,
                                              seurat_table_7,seurat_table_8,seurat_table_9,seurat_table_10,seurat_table_11) )

#####
#filter data
seurat_table[["percent.mt"]] <- PercentageFeatureSet(seurat_table, pattern = "^mt-")

x <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=nFeature_RNA ) ) +
  geom_hline( yintercept=300, color="red" ) + 
  geom_violin( fill="steelblue" ) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2, outlier.shape=NA) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("")

y <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=nCount_RNA ) ) +
  geom_violin( fill="pink" ) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2, outlier.shape=NA) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
  xlab("")

z <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=percent.mt ) ) +
  geom_hline( yintercept=10, color="red" ) +
  geom_violin( fill="orchid4" ) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2, outlier.shape=NA) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
  xlab("")

pdf( paste0("figures/qc_violin_prefilter_",date,".pdf" ), height=6, width=10 )
plot_grid( x,y,z, ncol=3 )
dev.off()

seurat_table <- subset(seurat_table, subset = percent.mt < 10)
seurat_table <- subset(seurat_table, subset = nFeature_RNA > 300)

#####
#normalize data
seurat_table <- NormalizeData(seurat_table, normalization.method = "LogNormalize", scale.factor = 10000)
nfeatures = 2000
seurat_table <- FindVariableFeatures(seurat_table, selection.method = "vst", nfeatures = nfeatures)
top10 <- head(VariableFeatures(seurat_table), 10)

#####
#run dimensionality reduction
all.genes <- rownames(seurat_table)
seurat_table <- ScaleData(seurat_table, features = VariableFeatures(object = seurat_table))
seurat_table <- RunPCA(seurat_table, features = VariableFeatures(object = seurat_table))

DimPlot(seurat_table, reduction = "pca", group.by = "orig.ident")

# Do elbow plot or JackStraw here
pdf( paste0("figures/PC_elbow_postfilter_",date,".pdf" ), height=5, width=6 )
ElbowPlot(seurat_table, ndims=50)
dev.off()

#####
# cluster and umap
seurat_table <- FindNeighbors(seurat_table, dims = 1:30)
seurat_table <- FindClusters(seurat_table, resolution = 2)
seurat_table <- RunUMAP(seurat_table, dims = 1:30, seed.use = 123)

pdf( paste0("figures/umap_",date,".pdf" ), height=8, width=12 )
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1)
dev.off()

DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1,split.by='orig.ident')

pdf( paste0("figures/umap_mimetic-features-test_",date,".pdf" ), height=8, width=12 )
FeaturePlot(seurat_table, features = c("krt5","foxi1","neurod1","pou2f3",
                                       "myog", "dmd","drgx","epcam"), pt.size=0.3, order=T )
dev.off()


#### Predict doublets: 
id_counts = as.data.frame(table(seurat_table@meta.data[["orig.ident"]])) 

library(DoubletFinder)

doubletfinder_fn <- function(test,out_prefix){
  nExp <- round(ncol(test) * 0.01)  
  test <- doubletFinder_v3(test, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
  
  DF.name = colnames(test@meta.data)[grepl("DF.classification", colnames(test@meta.data))]
  pdf( paste0("figures/umap_doubletFinder_",out_prefix,'_',date,".pdf" ), height=15, width=25 )
  print(cowplot::plot_grid(ncol = 2, DimPlot(test, group.by = "orig.ident") + NoAxes(),
                           DimPlot(test, group.by = DF.name) + NoAxes()))
  dev.off()
  return(test)
}

#separate seurat object by individual run to use doublet finder
for (donor in c("4wpfThymus_1", "4wpfThymus_3", "4wpfThymus_5", "4wpfThymus_6", "adult_thy_liberase_1",
                "adult_thy_liberase_2", "adult_thy_untreated", "Fish_22_Thymus_","Fish_24_Thymus_",
                "lck-eGFP_Thymus_1_", "lck-eGFP_Thymus_2_")){
  test = doubletfinder_fn(subset(seurat_table,subset=orig.ident==donor),donor)
  
  DF.name = colnames(test@meta.data)[grepl("DF.classification", colnames(test@meta.data))] #name of column in metadata where doubletfinder (DF) results stored
  test = as.data.frame(test@meta.data[DF.name]) #extract metadata, save as dataframe
  colnames(test) <- 'df'
  
  if(donor=='4wpfThymus_1'){ merged_metadata = test } #first one, save as variable
  if(donor!='4wpfThymus_1'){ merged_metadata = rbind(merged_metadata, test) } #after the first, concatenate variable
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
pdf( paste0("figures/umap_celltypes_allcells_",date,".pdf" ), height=25, width=25 )
FeaturePlot(seurat_table, features = c("rag1","rag2","ipcef1","il2rb","cd4-1","cd8a", #T cells
                                       "tcf7", #T
                                       "spock3","spi1a","ctsbb", #DCs
                                       "pdgfra","pdgfrb","twist1a","prrx1a","prrx1b","col2a1a", #fibroblasts
                                       "havcr2","mfap4", #macrophages
                                       "mpx","cpa5", #granulocytes
                                       "eomesa","fcer1gl", #NK cells
                                       "hbaa1","hbaa2", #erythrocytes
                                       "pax5","cd79a", #B cells
                                       "epcam","cdh1",'krt5'), pt.size=0.3, order=T ) #TEC
dev.off()

#remove contaminating cells
keep <- seurat_table@active.ident!=0 & seurat_table@active.ident!=1 & seurat_table@active.ident!=2 & 
  seurat_table@active.ident!=3 & seurat_table@active.ident!=4 & seurat_table@active.ident!=5 & 
  seurat_table@active.ident!=6 & seurat_table@active.ident!=8 & seurat_table@active.ident!=9 & 
  seurat_table@active.ident!=10 & seurat_table@active.ident!=11 #T
seurat_table <- seurat_table[,keep]
keep <- seurat_table@active.ident!=12 & seurat_table@active.ident!=26 & seurat_table@active.ident!=46 #erythrocytes 
seurat_table <- seurat_table[,keep]
keep <- seurat_table@active.ident!=16 & seurat_table@active.ident!=32 & seurat_table@active.ident!=39 #DC
seurat_table <- seurat_table[,keep]

keep <- seurat_table@active.ident!=24 #DC
seurat_table <- seurat_table[,keep]

keep <- seurat_table@active.ident!=15 & seurat_table@active.ident!=17 & seurat_table@active.ident!=18 & 
  seurat_table@active.ident!=43 & seurat_table@active.ident!=48 & seurat_table@active.ident!=47 & 
  seurat_table@active.ident!=14 & seurat_table@active.ident!=19 & seurat_table@active.ident!=20 #T
seurat_table <- seurat_table[,keep]

keep <- seurat_table@active.ident!=57 & seurat_table@active.ident!=41 #NK
seurat_table <- seurat_table[,keep]

keep <- seurat_table@active.ident!=29 & seurat_table@active.ident!=31 #granulocytes
seurat_table <- seurat_table[,keep]

keep <- seurat_table@active.ident!=53 #T
seurat_table <- seurat_table[,keep]
#
keep <- seurat_table@active.ident!=33 #macrophage
seurat_table <- seurat_table[,keep]

keep <- seurat_table@active.ident!=56 #NK
seurat_table <- seurat_table[,keep]

keep <- seurat_table@active.ident!=21 & seurat_table@active.ident!=23 &
  seurat_table@active.ident!=13 & seurat_table@active.ident!=50 & 
  seurat_table@active.ident!=38 & seurat_table@active.ident!=7 #T
seurat_table <- seurat_table[,keep]

keep <- seurat_table@active.ident!=36 & 
  seurat_table@active.ident!=28 & seurat_table@active.ident!=52 & 
  seurat_table@active.ident!=45 
seurat_table <- seurat_table[,keep]


############
pdf( paste0("figures/umap_aftercut1_",date,".pdf" ), height=8, width=12 )
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1)
dev.off()

## repeat normalization and plot to look for remaining contaminating cells
#####
#normalize data
norm_red <- function(seurat_table,nfeatures = 2000){
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
  
  return(seurat_table)
}

#repeat normalization and reduction after removing contaminating cells
seurat_table <- norm_red(seurat_table)

pdf( paste0("figures/umap_afterfilter1_",date,".pdf" ), height=8, width=12 )
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1)
dev.off()


pdf( paste0("figures/umap_celltypes_afterfilter1_",date,".pdf" ), height=25, width=25 )
FeaturePlot(seurat_table, features = c("rag1","rag2","ipcef1","il2rb","cd4-1","cd8a", #T cells
                                       "tcf7", #T
                                       "spock3","spi1a","ctsbb", #DCs
                                       "pdgfra","pdgfrb","twist1a","prrx1a","prrx1b","col2a1a", #fibroblasts
                                       "havcr2","mfap4", #macrophages
                                       "mpx","cpa5", #granulocytes
                                       "eomesa","fcer1gl", #NK cells
                                       "hbaa1","hbaa2", #erythrocytes
                                       "pax5","cd79a", #B cells
                                       "epcam","cdh1",'krt5'), pt.size=0.3, order=T ) #TEC
dev.off()


## Look for doublets ###########
pdf( paste0("figures/vln_counts_afterfilter1_", date,".pdf"), height=8, width=20 )
VlnPlot(seurat_table, features = c('nFeature_RNA','nCount_RNA'))
dev.off()

pdf( paste0("figures/umap_doubletFinder_afterfilter1_",date,".pdf" ), height=15, width=25 )
print(cowplot::plot_grid(ncol = 2, DimPlot(seurat_table, group.by = "orig.ident") + NoAxes(),
                         DimPlot(seurat_table, group.by = 'Doublet.Finder') + NoAxes()))
dev.off()

#########
saveRDS(seurat_table, file = "Rdata/zebrafish_mec_seurat_table_20231017_afterfilter1-wDoubletFinder.rds")

