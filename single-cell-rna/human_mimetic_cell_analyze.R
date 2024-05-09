#####
# Load saved RDS and continue analyses
#####
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

seurat_table <- readRDS( "Rdata/meclo_ht2-ht6_seurat_table_20230511-clust34removed.rds" )

date='20240213'

cmap <- c("#D7B5A6",
          "#4E79A7","#D37295","#A0CBE8",
          "#F28E2B","#FFBE7D","#59A14F","#8CD17D","#B6992D","#F1CE63",
          "#D4A6C8","#86BCB6","#E15759","#FF9D9A","#499894","#BAB0AC",
          "#FABFD2","#B07AA1","#9D7660"
)

pdf( paste0("figures/umap_orig-labels_",date,".pdf" ), height=8, width=10 ) 
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1 ) 
dev.off()

## Merge and rename clusters based on characteristics
keep <- seurat_table@active.ident==21 | seurat_table@active.ident==6 | seurat_table@active.ident==4 | seurat_table@active.ident==25 
seurat_table@active.ident[keep] <- 4 #Ionocyte
#
keep <- seurat_table@active.ident==16 | seurat_table@active.ident==17 
seurat_table@active.ident[keep] <- 16 #ATOH1+ neuro
# 
keep <- seurat_table@active.ident==23 | seurat_table@active.ident==8 
seurat_table@active.ident[keep] <- 23 #NKX6-2+ neuro
#
keep <- seurat_table@active.ident==14    | seurat_table@active.ident==12 | seurat_table@active.ident==20
seurat_table@active.ident[keep] <- 14 #SHOX2+ neuro 
# 
keep <- seurat_table@active.ident==26 | seurat_table@active.ident==31 
seurat_table@active.ident[keep] <- 26 #Tuft 
# 
keep <- seurat_table@active.ident==7 | seurat_table@active.ident==19
seurat_table@active.ident[keep] <- 7 #Late muscle
# 
keep <- seurat_table@active.ident==0 | seurat_table@active.ident==2 | seurat_table@active.ident==3 | seurat_table@active.ident==24 | seurat_table@active.ident==33 | seurat_table@active.ident==5
seurat_table@active.ident[keep] <- 0 #Early muscle
# 
keep <- seurat_table@active.ident==1 | seurat_table@active.ident==11 | seurat_table@active.ident==30 
seurat_table@active.ident[keep] <- 1 #CCL19
# 
keep <- seurat_table@active.ident==9 | seurat_table@active.ident==18 | seurat_table@active.ident==15 | 
  seurat_table@active.ident==10 
seurat_table@active.ident[keep] <- 9 #Intermediate
# 
seurat_table <- RenameIdents( object=seurat_table,
                              "0"="Early muscle",
                              "7"="Late muscle",
                              "22"="lncRNA-enriched muscle",
                              "13"="HMX2+ ionocyte",
                              "4" ="Ionocyte", #ASCL3+ ionocyte
                              "26"="Tuft",
                              "16"="Cochlear hair cell",#"ATOH1+ neuro",
                              "23"="NKX6-2+ neuro",
                              "14"="SHOX2+ neuro",
                              "29"="FEZF2+ neuro",
                              "32"="DRGX+ sensory neuro", 
                              "1"="Immature", 
                              "34"="Keratinocyte",
                              "35"='Cycling',
                              "28"='CUX2+ neuro',
                              "9"='Intermediate', #Transitional
                              "27"="Aire-stage")

print(levels(seurat_table@active.ident))

seurat_table@active.ident <- factor(seurat_table@active.ident,
                                    levels=c( "Aire-stage","Immature","Intermediate","Cycling",
                                              "Keratinocyte",
                                              "HMX2+ ionocyte","Ionocyte",
                                              "Tuft",
                                              "Early muscle", "Late muscle", "lncRNA-enriched muscle",
                                              "CUX2+ neuro","FEZF2+ neuro","NKX6-2+ neuro", "SHOX2+ neuro",#'PHOX2B+ neuro',
                                              "DRGX+ sensory neuro", "Cochlear hair cell"
                                    ))

print(levels(seurat_table@active.ident))


####################################
pdf( paste0("figures/umap_relabeled_orig_",date,".pdf" ), height=8, width=10 ) 
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1, cols=cmap ) + theme(legend.position="none")
dev.off()

############################
# Refine boundary of intermediate and ionocyte clusters 
# look more closely at what was originally cluster 9
############### Repeat normalization, find var features, etc
#normalize data
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

keep <- seurat_table@meta.data[["seurat_clusters"]]=='9'

pdf( paste0("figures/umap_cluster9_beforesubclustered_",date,".pdf" ), height=4, width=5 ) #height=10, width=12 )
DimPlot(seurat_table[,keep], reduction="umap", label=T, repel=T, pt.size=1)
dev.off()

seurat_9 <- process_fn(seurat_table[,keep],2000,1)#0.5)

pdf( paste0("figures/umap_cluster9_",date,".pdf" ), height=4, width=5 ) #height=10, width=12 )
DimPlot(seurat_9, reduction="umap", label=T, repel=T, pt.size=1)
dev.off()

#add new cluster names to original object
seurat_table$clust9 <- as.character(Idents(seurat_table))
seurat_table$clust9[Cells(seurat_9)] <- paste("c1",Idents(seurat_9))
pdf( paste0("figures/umap_clust9_subclusters_",date,".pdf" ), height=10, width=14 ) #height=10, width=12 )
DimPlot(seurat_table, group.by = "clust9", label=T, repel=T)
dev.off()

FeaturePlot(seurat_9, features=c('FOXI1','POU2F3','IL18','CFTR','ASCL3','HMX2','AQP6','CGA'), order=TRUE, pt.size = 0.1) &NoLegend()# + theme(legend.position="none")

# assign subcluster6 from cluster9 to the ionocyte cluster
keep <- seurat_9@meta.data[["seurat_clusters"]]==6 
seurat_table@active.ident[Cells(seurat_9[,keep])] <- 'Ionocyte'

print(levels(seurat_table@active.ident))

################
####################################
pdf( paste0("figures/umap_relabeled_",date,".pdf" ), height=8, width=10 ) 
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1, cols=cmap ) + theme(legend.position="none")
dev.off()

pdf( paste0("figures/umap_relabeled_split_",date,".pdf" ), height=5, width=20 )
DimPlot(seurat_table, reduction="umap", label=F, repel=T, pt.size=1, cols=cmap, split.by = 'orig.ident' ) + theme(legend.position="none")
dev.off()

pdf( paste0("figures/umap_relabeled_notext_",date,".pdf" ), height=8, width=10 )
DimPlot(seurat_table, reduction="umap", label=F, repel=T, pt.size=1, cols=cmap) + theme(legend.position="none")
dev.off()

################
#Do this to reorder donors
seurat_table@meta.data$orig.identABC = ''
seurat_table@meta.data$orig.identABC[seurat_table@meta.data$orig.ident=='HT2'] = 'D'
seurat_table@meta.data$orig.identABC[seurat_table@meta.data$orig.ident=='HT3'] = 'C'
seurat_table@meta.data$orig.identABC[seurat_table@meta.data$orig.ident=='HT4'] = 'A'
seurat_table@meta.data$orig.identABC[seurat_table@meta.data$orig.ident=='HT5'] = 'B'
seurat_table@meta.data$orig.identABC[seurat_table@meta.data$orig.ident=='HT6'] = 'E'

#Plot distributions of cells among clusters
#group by cluster
pdf( paste0("figures/cluster_distribution_relabeled-clusters_",date,".pdf" ), height=5, width=12 )
cluster_dist = t(table(seurat_table@active.ident, seurat_table@meta.data$orig.identABC))
par(mar=c(12,4,4,4)) #margins
barplot(cluster_dist,beside=TRUE,ylab='Number of cells',las=2, col=rep(cmap, each=5))
dev.off()

#split by ABC
pdf( paste0("figures/umap_relabeled_split_ABCorder_",date,".pdf" ), height=5, width=20 )
DimPlot(seurat_table, reduction="umap", label=F, repel=T, pt.size=1,split.by = 'orig.identABC' , cols=cmap) + theme(legend.position="none")
dev.off()

DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1, cols=cmap ) + theme(legend.position="none")

################
## dotplot of TFs for clusters  ##########
pdf( paste0("figures/dotplot_TFs_",date,".pdf"), height=6, width=10 )
DotPlot(seurat_table, features=(c('AIRE','IRF8','TP63',
                                  "GRHL1","BARX2","RUNX2",
                                  "TFCP2L1",
                                  "FOXI1","HMX2","ASCL3",
                                  "POU2F3",
                                  "MYOG","MYOD1","MEF2C",
                                  'CUX2',
                                  "FEZF2",
                                  "NKX6-2",
                                  'SHOX2',
                                  'DRGX',
                                  "ATOH1")),
        cols=c("lightgray","#3F007D"), dot.min = 0.05, dot.scale=8 ) + 
  theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust=1)) +
  theme(axis.text.x = element_text(face = "italic"))+
  scale_y_discrete(limits=rev) #reverse y-direction
dev.off()

#### 
#downsample cells to plot in heatmap:
print(min(table(Idents(seurat_table)))) # to check count of cells for the ident with the fewest

#save heatmap of top cluster markers -- curate features to plot in heatmap
heatmap_features = c("HLA-DQA1",'HLA-DRB1',
                     "CCL19","KRT14",'KRT17','FN1',
                     'SGCZ',
                     "MKI67","CENPF","ASPM","DEPDC1", 
                     "KRT5","KRT6A","IVL","MUC4","NMU",
                     'IL18','SLC26A4','SLC4A9','INSRR', 
                     'CFTR','ATP6V0D2', 
                     "PRH1","TRPM5",'CHAT',  
                     "RGS21", "PLCB2","AZGP1", 
                     'MYL5',"NEB",'CKM',"DLK1","TTN","TNNT3","TNNC2", 
                     "LINC00326","AC004949.1",'GPC5', 
                     "VWC2","JPH4", 
                     "CADPS","EYS", 
                     "SYT1","FXYD3",'NRG1',
                     "MEGF11", 
                     "SLC6A5","NTRK3", 'FSTL5','GRIK1','GHRH',
                     "USH2A","OTOF") 

#slot=default scaled data
breakslist = seq(1,3,by=0.01) 
pdf( paste0("figures/markers_heatmap_selected_",date,".pdf" ), height=8.5, width=8 ) #height=10, width=7 )
DoHeatmap(subset(seurat_table, downsample=100), features = heatmap_features,angle=90, disp.min=0.5, size=4, group.colors=cmap) + theme(axis.text.y.left = element_text(face = "italic"))+ #size changes cluster label size
  scale_fill_gradientn(colors = colorRampPalette((brewer.pal(n=9, name="Purples")))(length(breakslist)))
dev.off()

##########################################################################################
# Subcluster neuroendocrine clusters
keep <- seurat_table@active.ident=='Cochlear hair cell'|seurat_table@active.ident=='NKX6-2+ neuro'|
  seurat_table@active.ident=='FEZF2+ neuro' |seurat_table@active.ident=='CUX2+ neuro'|
  seurat_table@active.ident=='DRGX+ sensory neuro'|seurat_table@active.ident=='SHOX2+ neuro'

seurat_neuro <- process_fn(seurat_table[,keep],2000,1)

pdf( paste0("figures/umap_neuro_subclustered_",date,".pdf" ), height=4, width=5 ) #height=10, width=12 )
DimPlot(seurat_neuro, reduction="umap", label=T, repel=T, pt.size=1)
dev.off()

#add orig cluster names to subclustered object
seurat_neuro$sub_cluster <- as.character(Idents(seurat_neuro))
seurat_neuro$sub_cluster[Cells(seurat_table)] <- paste("c1",Idents(seurat_table))
pdf( paste0("figures/umap_neuro_subclustered_orignames1_",date,".pdf" ), height=4, width=7 ) #height=10, width=12 )
DimPlot(seurat_neuro, group.by = "sub_cluster", label=T, repel=T, pt.size=1, cols = c("#FABFD2","#86BCB6","#BAB0AC","#E15759","#FF9D9A","#499894","#B07AA1"))#&NoLegend()
dev.off()

#remove cluster 14, 15
markers_14 <- FindMarkers(seurat_neuro, ident.1=c(14), min.pct=0.1, logfc.threshold=1, only.pos=T, test.use="wilcox")
markers_15 <- FindMarkers(seurat_neuro, ident.1=c(15), min.pct=0.1, logfc.threshold=1, only.pos=T, test.use="wilcox")
markers_14_15 <- FindMarkers(seurat_neuro, ident.1=c(14,15), min.pct=0.1, logfc.threshold=1, only.pos=T, test.use="wilcox")

pdf( paste0("figures/umap_neuro_features_",date,".pdf" ), height=6, width=5 ) #height=10, width=12 )
FeaturePlot(seurat_neuro, features=c('CCL19','TP63','CCL21','KRT15','CFTR'), order=TRUE, pt.size = 0.1) &NoLegend()# + theme(legend.position="none")
dev.off()

keep <- seurat_neuro@active.ident=='14'|seurat_neuro@active.ident=='15'

seurat_neuro2 <- process_fn(seurat_neuro[,!keep],2000,1)

pdf( paste0("figures/umap_neuro_subclustered2_",date,".pdf" ), height=4, width=5 ) #height=10, width=12 )
DimPlot(seurat_neuro2, reduction="umap", label=T, repel=T, pt.size=1)
dev.off()

#add orig cluster names to subclustered object
seurat_neuro2$sub_cluster <- as.character(Idents(seurat_neuro2))
seurat_neuro2$sub_cluster[Cells(seurat_table)] <- paste("c1",Idents(seurat_table))
pdf( paste0("figures/umap_neuro_subclustered_orignames2_",date,".pdf" ), height=4, width=7 ) #height=10, width=12 )
DimPlot(seurat_neuro2, group.by = "sub_cluster", label=T, repel=T, pt.size=1, cols = c("#FABFD2","#86BCB6","#BAB0AC","#E15759","#FF9D9A","#499894","#B07AA1"))#&NoLegend()
dev.off()

pdf( paste0("figures/umap_neuro_subclustered_orignames2_notext_",date,".pdf" ), height=4, width=7 ) #height=10, width=12 )
DimPlot(seurat_neuro2, group.by = "sub_cluster", label=F, repel=T, pt.size=1, cols = c("#FABFD2","#86BCB6","#BAB0AC","#E15759","#FF9D9A","#499894","#B07AA1"))#&NoLegend()
dev.off()

#TFs
pdf( paste0("figures/features-neuroTF_",date,".pdf" ), height=12, width=17 ) 
FeaturePlot(seurat_neuro2, features=c('DRGX',
                                      'CUX2',
                                      "FEZF2",'EN1','HAND1','SIM1',
                                      "PHOX2B","CREB5","UNCX",
                                      'MYT1L','NHLH2','NKX2-2',
                                      'LEF1','HLX',
                                      "NKX6-2",
                                      'NEUROD4',
                                      'SHOX2',
                                      "NEUROD1","EBF3",
                                      "ZNF536","SOX11","POU4F1",
                                      "IRX2","SOX2",
                                      "ATOH1","BARHL1","POU4F3","PAX2","LHX3","GFI1"
), order=TRUE, pt.size = 0.5,ncol=6, cols=c("lightgray","#3F007D")) &
  NoLegend() & NoAxes() + theme(plot.title = element_text( face = "bold.italic"))
dev.off()

####################################################################
#split by age, sex of donor
#separate donor with DS, as this has skewed distribution

#individual donor calculations
cluster_dist = t(table(seurat_table@active.ident, seurat_table@meta.data$orig.ident)) #number of cells for each donor
cluster_df = as.data.frame(cluster_dist)
names(cluster_df)[1] = 'donor'
names(cluster_df)[2] = 'celltype'
cluster_df = cluster_df[cluster_df['donor']!='HT4',] #exclude DS donor which has skewed cluster distribution
cluster_df

cluster_df$sex = '_'
cluster_df[cluster_df$donor=='HT2','sex'] = 'Male'
cluster_df[cluster_df$donor=='HT3','sex'] = 'Female'
cluster_df[cluster_df$donor=='HT4','sex'] = 'Female (DS)'
cluster_df[cluster_df$donor=='HT5','sex'] = 'Male'
cluster_df[cluster_df$donor=='HT6','sex'] = 'Female'

cluster_df$age = '_'
cluster_df[cluster_df$donor=='HT2','age'] = '3-4Mo'
cluster_df[cluster_df$donor=='HT3','age'] = '2-4D'
cluster_df[cluster_df$donor=='HT4','age'] = '3-4Mo (DS)'
cluster_df[cluster_df$donor=='HT5','age'] = '2-4D'
cluster_df[cluster_df$donor=='HT6','age'] = '3-4Mo'
cluster_df

avg_sex = aggregate(cluster_df$Freq, list(cluster_df$celltype, cluster_df$sex), FUN=mean)
avg_age = aggregate(cluster_df$Freq, list(cluster_df$celltype, cluster_df$age), FUN=mean)

#sex
pdf( paste0("figures/cluster_distribution_sex_",date,".pdf" ), height=2.5, width=12 )
ggplot(cluster_df, aes(x = factor(celltype), y = Freq, fill = sex)) +
  geom_col(data=avg_sex, aes(x = factor(Group.1), y = x, fill = factor(Group.2)), color = "black", size = 1, width = .8, position = "dodge") +
  geom_point(position = position_jitterdodge(jitter.width = 0), size=1.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_fill_manual(values = c("grey", "white")) + ylab('Number of cells') + 
  theme(axis.text.y = element_text(color = "black", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"))
dev.off()  

#age
pdf( paste0("figures/cluster_distribution_age_",date,".pdf" ), height=2.5, width=12 )
ggplot(cluster_df, aes(x = factor(celltype), y = Freq, fill = age)) +
  geom_col(data=avg_age, aes(x = factor(Group.1), y = x, fill = factor(Group.2)), color = "black", size = 1, width = .8, position = "dodge") +
  geom_point(position = position_jitterdodge(jitter.width = 0), size=1.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_fill_manual(values = c("grey", "white")) + ylab('Number of cells') + 
  theme(axis.text.y = element_text(color = "black", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"))
dev.off()  


#Split UMAP
#separate donor with DS, as this has skewed distribution
seurat_table@meta.data$orig.identSex = ''
seurat_table@meta.data$orig.identSex[seurat_table@meta.data$orig.ident=='HT2'] = 'Male'
seurat_table@meta.data$orig.identSex[seurat_table@meta.data$orig.ident=='HT3'] = 'Female'
seurat_table@meta.data$orig.identSex[seurat_table@meta.data$orig.ident=='HT4'] = 'X (DS) Female'
seurat_table@meta.data$orig.identSex[seurat_table@meta.data$orig.ident=='HT5'] = 'Male'
seurat_table@meta.data$orig.identSex[seurat_table@meta.data$orig.ident=='HT6'] = 'Female'

seurat_table@meta.data$orig.identAge = ''
seurat_table@meta.data$orig.identAge[seurat_table@meta.data$orig.ident=='HT2'] = '3-4Mo'
seurat_table@meta.data$orig.identAge[seurat_table@meta.data$orig.ident=='HT3'] = '2-4D'
seurat_table@meta.data$orig.identAge[seurat_table@meta.data$orig.ident=='HT4'] = 'X (DS) 3-4Mo'
seurat_table@meta.data$orig.identAge[seurat_table@meta.data$orig.ident=='HT5'] = '2-4D'
seurat_table@meta.data$orig.identAge[seurat_table@meta.data$orig.ident=='HT6'] = '3-4Mo'

#UMAP split by sex
pdf( paste0("figures/umap_relabeled_split_sex_",date,".pdf" ), height=5, width=12 )
DimPlot(seurat_table, reduction="umap", label=F, repel=T, pt.size=1,split.by = 'orig.identSex' , cols=cmap) + theme(legend.position="none")
dev.off()

#UMAP split by age
pdf( paste0("figures/umap_relabeled_split_age_",date,".pdf" ), height=5, width=12 )
DimPlot(seurat_table, reduction="umap", label=F, repel=T, pt.size=1,split.by = 'orig.identAge' , cols=cmap) + theme(legend.position="none")
dev.off()

##########################################################################################
#tuft markers
pdf( paste0("figures/features-tuft-related_",date,".pdf" ), height=3, width=12.5 ) 
FeaturePlot(seurat_table, features=c('DCLK1', 'IL19','IL37','NMU','IL25'), order=TRUE, pt.size = 0.3,ncol=5, cols=c("lightgray","#3F007D")) &
  NoLegend() & NoAxes() + theme(plot.title = element_text( face = "bold.italic")) # + theme(legend.position="none")
dev.off()

pdf( paste0("figures/features-tuft-related2_",date,".pdf" ), height=6, width=12.5 ) 
FeaturePlot(seurat_table, features=c('IL37',#'IL17RB',
                                     'IL19','IL25',
                                     'NMU','DCLK1',
                                     'PLCB2','CHAT',
                                     'ALOX5AP','ALOX5',
                                     'GNAT3'
                                     ), order=TRUE, pt.size = 0.3,ncol=5, cols=c("lightgray","#3F007D")) &
  NoLegend() & NoAxes() + theme(plot.title = element_text( face = "bold.italic")) # + theme(legend.position="none")
dev.off()

##########################################################################################
#ionocyte markers
pdf( paste0("figures/features-ionocyte-related_",date,".pdf" ), height=3, width=12.5 ) 
FeaturePlot(seurat_table, features=c('HMX2','ASCL3','SLC26A4','SLC4A9','INSRR'), order=TRUE, pt.size = 0.3,ncol=5, cols=c("lightgray","#3F007D")) &
  NoLegend() & NoAxes() + theme(plot.title = element_text( face = "bold.italic")) 
dev.off()

##########################################################################################
#subcluster muscle clusters
keep <- seurat_table@active.ident=='Early muscle'|seurat_table@active.ident=='Late muscle'|
  seurat_table@active.ident=='lncRNA-enriched muscle' 
pdf( paste0("figures/umap_muscle_beforesubclustered_",date,".pdf" ), height=4, width=5 ) 
DimPlot(seurat_table[,keep], reduction="umap", label=T, repel=T, pt.size=1)
dev.off()

seurat_muscle <- process_fn(seurat_table[,keep],2000,0.5)
pdf( paste0("figures/umap_muscle_subclustered_",date,".pdf" ), height=4, width=5 ) 
DimPlot(seurat_muscle, reduction="umap", label=T, repel=T, pt.size=1)
dev.off()
#####
## clusters 8+7 are enriched for ionocyte/intermediate genes -- exclude and recluster
pdf( paste0("figures/markers_muscle_check_",date,".pdf" ), height=8, width=8 ) 
FeaturePlot(seurat_muscle, features=c('FOXI1','CCL19','CCL21', 'MYOG', 'MYOD1'), order=TRUE, pt.size = 0.1) &NoLegend()# + theme(legend.position="none")
dev.off()

#subcluster muscle clusters
keep <- seurat_muscle@active.ident==8|seurat_muscle@active.ident==7
DimPlot(seurat_muscle[,!keep], reduction="umap", label=T, repel=T, pt.size=1)

#Use 0.3, for 5 clusters
seurat_muscle2 <- process_fn(seurat_muscle[,!keep],2000,0.3) 
DimPlot(seurat_muscle2, reduction="umap", label=T, repel=T, pt.size=1)+ theme(legend.position="none") + NoAxes()

#Flip over axes for aesthetics (late muscle on the right)
seurat_muscle2@reductions[["umap"]]@cell.embeddings[,'UMAP_1'] = -1*seurat_muscle2@reductions[["umap"]]@cell.embeddings[,'UMAP_1']
seurat_muscle2@reductions[["umap"]]@cell.embeddings[,'UMAP_2'] = -1*seurat_muscle2@reductions[["umap"]]@cell.embeddings[,'UMAP_2']

## relabel muscle clusters so that in this desired order for heatmap:
seurat_muscle2 <- RenameIdents( object=seurat_muscle2,
                                "2"="Early", #early
                                "0"="Inter1", #inter 1
                                "1"="Inter2", #inter 2
                                "3"="Late", #late
                                "4"="lncRNA-enriched") #lnc


cmap2 = c("#4E79A7","#F28E2B",  "#59A14F", "#B07AA1", "#FABFD2")

pdf( paste0("figures/umap_muscle_",date,".pdf" ), height=4, width=3.33)#5 )
DimPlot(seurat_muscle2, reduction="umap", label=F, repel=T, pt.size=1, cols=cmap2)+ theme(legend.position="none") + NoAxes()
dev.off()

pdf( paste0("figures/umap_muscle_label_",date,".pdf" ), height=4, width=5 )
DimPlot(seurat_muscle2, reduction="umap", label=T, repel=T, pt.size=1, cols=cmap)+ theme(legend.position="none")
dev.off()

###### 
# HEATMAP of muscle cluster genes to include in supplement: 
markers <- FindAllMarkers(seurat_muscle2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
heatmap_features = top10$gene

breakslist = seq(1,3,by=0.01)
pdf( paste0("figures/markers_muscle_heatmap_selected_purple_",date,".pdf" ), height=10, width=7 )
DoHeatmap(subset(seurat_muscle2, downsample=50), features = heatmap_features,angle=90, disp.min=0.5, size=4, group.colors=cmap2) + theme(axis.text.y.left = element_text(face = "italic"))+ #NoLegend() + #size changes cluster label size
  scale_fill_gradientn(colors = colorRampPalette((brewer.pal(n=9, name="Purples")))(length(breakslist))) 
dev.off()

##############
# Muscle feature plots 
#muscle differentiation genes, selected
pdf( paste0("figures/features-muscle_mtg_",date,".pdf" ), height=3, width=20 ) 
FeaturePlot(seurat_muscle2, features=c('CCL19','SLIT3','MYMX','MYMK','DYSF','DES','DMD','DLK1'), order=TRUE, pt.size = 0.5,ncol=8, cols=c("lightgray","#3F007D")) &
  NoLegend() & NoAxes() + theme(plot.title = element_text( face = "bold.italic")) # + theme(legend.position="none")
dev.off()

#canonical muscle differentiation genes
pdf( paste0("figures/features_classic-muscle_",date,".pdf" ), height=3, width=20 )
FeaturePlot(seurat_muscle2, features=c('PAX3', 'PAX7', 'MYF5','MYOD1', 'MYOG', 'MEF2A', 'MEF2C','MYF6'),ncol=8, order=TRUE, pt.size = 0.5, cols=c("lightgray","#3F007D")) &
  NoLegend() & NoAxes() + theme(plot.title = element_text( face = "bold.italic"))
dev.off()

#developmental myosins: https://skeletalmusclejournal.biomedcentral.com/articles/10.1186/s13395-015-0046-6/tables/1
pdf( paste0("figures/features_myosins_",date,".pdf" ), height=10, width=14 )
FeaturePlot(seurat_muscle2, features=c('MYH1', 'MYH2', 'MYH3', 'MYH4', 'MYH6', 
                                       'MYH7', 'MYH8', 'MYH7B', 'MYH13', 'MYH15', 'MYH16',     
                                       'MYL1',  'MYL2', 'MYL3', 'MYL4',  'MYL6', 
                                       'MYL6B', 'MYLPF'), order=TRUE, pt.size = 0.5,ncol=5)& 
  NoLegend() & NoAxes() + theme(plot.title = element_text( face = "bold.italic"))
dev.off()

#troponins
pdf( paste0("figures/features_troponins_",date,".pdf" ), height=5, width=10 )
FeaturePlot(seurat_muscle2, features=c('TNNI1', 'TNNI2', 'TNNI3', 'TNNT1', 'TNNT2', 'TNNT3', 'TNNC1', 'TNNC2'), order=TRUE, pt.size = 1,ncol=4)& theme(legend.position="none")+ theme(plot.title = element_text( face = "bold.italic"))
dev.off()


################################################################################
### project muscle reference data signatures onto umap
# HARMONY cluster signatures ################
setwd("/McKellar/text_outputs/signatures_human_conversions/")
idx <- list.files()[grepl("harmony", list.files())]
sig_list <- list()
for (i in idx ) {
  sig <- read.delim(i,header=T, sep="\t")
  sig <- sig[,1]
  sig_name <- substr(i, 0, nchar(i)-40)
  sig_list[[sig_name]] <- sig
  seurat_muscle2 <- AddModuleScore(seurat_muscle2, features=as.data.frame(sig), name=sig_name)
}
#combine signatures from 'Myonuclei Type IIb' + 'Myonuclei Type IIx'
sig_list[['harmonyMyonuclei (merged)']] <- unique(c(sig_list[['harmonyMyonuclei (Type IIx)']], sig_list[['harmonyMyonuclei (Type IIb)']] ))
seurat_muscle2 <- AddModuleScore(seurat_muscle2, features=as.data.frame(sig_list[['harmonyMyonuclei (merged)']]), name='harmonyMyonuclei (merged)')

setwd( "/Analysis1" )
pdf( paste0("figures/mckellar_harmony_muscle-signatures_seurat-muscle_",date,".pdf" ), height=3, width=20 )
FeaturePlot(seurat_muscle2, features=c('harmonyMuSCs1',"harmonyMyoblasts_Progenitors1",
                                       'harmonyMyonuclei..merged.1'
                                       ), order=T, min.cutoff="q75",pt.size=0.5, cols=c("lightgray","#3F007D"), ncol=8) &
  NoLegend() & NoAxes() + theme(plot.title = element_text( face = "bold.italic"))
dev.off()
#


# Cardiac, Feng et al cluster signatures ################
setwd("/Feng cardiac/text_outputs/signatures_human_conversions/")
idx <- list.files()[grepl("signature", list.files())]
sig_list <- list()
for (i in idx ) {
  sig <- read.delim(i,header=T, sep="\t")
  sig <- sig[,1]
  sig_name <- substr(i, 11, nchar(i)-40)
  sig_list[[sig_name]] <- sig
  seurat_muscle2 <- AddModuleScore(seurat_muscle2, features=as.data.frame(sig), name=sig_name)
}
setwd( "/Analysis1" )
pdf( paste0("figures/feng_cardiac_muscle-signatures_seurat-muscle_",date,".pdf" ), height=3, width=20 )
FeaturePlot(seurat_muscle2, features=c('atrial_cm1','ventricular_cm1'
                                       ), order=T, min.cutoff="q75",pt.size=0.5, cols=c("lightgray","#3F007D"), ncol=8) &
  NoLegend() & NoAxes() + theme(plot.title = element_text( face = "bold.italic"))
dev.off()

################################# Save files to read into scvelo 
library(SeuratDisk)
library(Matrix)
date='20240110'

# save metadata table:
seurat_muscle2$barcode <- colnames(seurat_muscle2)
seurat_muscle2$UMAP_1 <- seurat_muscle2@reductions$umap@cell.embeddings[,1]
seurat_muscle2$UMAP_2 <- seurat_muscle2@reductions$umap@cell.embeddings[,2]
seurat_muscle2$clusterLabels <- seurat_muscle2@active.ident 
write.csv(seurat_muscle2@meta.data, file= paste0(date,"_metadata_muscle",".csv"), quote=F, row.names=F)

# write expression counts matrix
counts_matrix <- GetAssayData(seurat_muscle2, assay='RNA', slot='counts')
writeMM(counts_matrix, file=paste0(date,'_counts_muscle','.mtx'))

# write dimensionality reduction matrix
write.csv(seurat_muscle2@reductions$pca@cell.embeddings, file=paste0(date,'_pca_muscle','.csv'), quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file=paste0(date,'_gene_names_muscle','.csv'),
  quote=F,row.names=F,col.names=F
)

################################################################################################
# Pseudo bulk sc data and perform edgeR or DESEQ differential analysis: ###################
# Similar to: https://github.com/dmichelson/thymic_hnf4/blob/main/scRNAseq/hnf4ctrl_vs_hnf4dtec_scrnaseq.R

seurat_table$clusterLabels <- seurat_table@active.ident #save updated cluster names in metadata
seurat_table@meta.data[["muscle clusters"]] = 
  seurat_table@meta.data[["clusterLabels"]]=='lncRNA-enriched muscle' | 
  seurat_table@meta.data[["clusterLabels"]]=='Late muscle' |
  seurat_table@meta.data[["clusterLabels"]]=='Early muscle'

FeaturePlot(seurat_table, features=c('muscle clusters'), order=TRUE, pt.size = 1) #+ theme(legend.position="none")

dat <- cbind( rowSums(seurat_table@assays$RNA@counts[,seurat_table$orig.ident=="HT2" & seurat_table@meta.data[["muscle clusters"]]==TRUE]),
              rowSums(seurat_table@assays$RNA@counts[,seurat_table$orig.ident=="HT3" & seurat_table@meta.data[["muscle clusters"]]==TRUE]),
              rowSums(seurat_table@assays$RNA@counts[,seurat_table$orig.ident=="HT4" & seurat_table@meta.data[["muscle clusters"]]==TRUE]),
              rowSums(seurat_table@assays$RNA@counts[,seurat_table$orig.ident=="HT5" & seurat_table@meta.data[["muscle clusters"]]==TRUE]),
              rowSums(seurat_table@assays$RNA@counts[,seurat_table$orig.ident=="HT6" & seurat_table@meta.data[["muscle clusters"]]==TRUE]),
              rowSums(seurat_table@assays$RNA@counts[,seurat_table$orig.ident=="HT2" & seurat_table@meta.data[["muscle clusters"]]==FALSE]),
              rowSums(seurat_table@assays$RNA@counts[,seurat_table$orig.ident=="HT3" & seurat_table@meta.data[["muscle clusters"]]==FALSE]),
              rowSums(seurat_table@assays$RNA@counts[,seurat_table$orig.ident=="HT4" & seurat_table@meta.data[["muscle clusters"]]==FALSE]),
              rowSums(seurat_table@assays$RNA@counts[,seurat_table$orig.ident=="HT5" & seurat_table@meta.data[["muscle clusters"]]==FALSE]),
              rowSums(seurat_table@assays$RNA@counts[,seurat_table$orig.ident=="HT6" & seurat_table@meta.data[["muscle clusters"]]==FALSE]))
colnames(dat) <- c("ht2_muscle","ht3_muscle","ht4_muscle","ht5_muscle","ht6_muscle",
                   "ht2_nonmuscle","ht3_nonmuscle","ht4_nonmuscle","ht5_nonmuscle","ht6_nonmuscle")

library(DESeq2)
idx <- substr( colnames(dat), 5, 10 ) == 'muscle'
ctrl_idx <- substr( colnames(dat), 5, 10 ) == 'nonmus'

group <- factor(c(rep(1,times=sum(ctrl_idx)),rep(2,times=sum(idx))))
coldat=DataFrame(group)

dds <- DESeqDataSetFromMatrix(countData = cbind( dat[,ctrl_idx], dat[,idx]), #Deseq2 #adjusted ordering so ctrl_idx treated as control
                              colData = coldat,
                              design = ~ group) #group in coldat$group
dds <- DESeq(dds) #Deseq2

dds.res_df <- as.data.frame(results(dds))
dds.res_df$genes = rownames(dds.res_df)

dds.res_df = na.omit(dds.res_df) 

tflist = read.csv('Lambert et al Cell 2018 Supp data S1_TFs.csv')
tflist$tfs=tflist$Transcription.Factors
dds.res_df = merge(dds.res_df, tflist, by.x='genes',by.y='Transcription.Factors',all.x=T,all.y=F)
dds.res_df$col <- 'black'
dds.res_df[!is.na(dds.res_df$tfs),'col'] <- 'red'


#label specific genes in volcano plot:
for (gene_i in c('MYOG','DYSF','MEF2C','MYOD1','MYF6','TTN','MYH3','DES','DMD','SYNPO2L','RYR1','MYH3','PITX2')){ #
  dds.res_df[dds.res_df$genes==gene_i,'tf_sig'] <- gene_i 
}

pdf( paste0("figures/vol_plot_muscle_pseudobulked-labels_",date,".pdf" ), height=10, width=10 ) 
ggplot(data=dds.res_df, aes(x=log2FoldChange, y=-log10(padj), col=col , label=tf_sig )) + 
  scale_color_manual(values=c("black", "red")) +
  geom_point() + theme_classic()  +  geom_text_repel(max.overlaps = 10, min.segment.length = 0.2,force = 100,seed = 13) + theme(text = element_text(size = 15), legend.position="none")  #
dev.off()

pdf( paste0("figures/vol_plot_muscle_pseudobulked-nolabels_",date,".pdf" ), height=4, width=5 ) 
ggplot(data=dds.res_df %>% arrange(col), aes(x=log2FoldChange, y=-log10(padj), col=col)) + 
  scale_color_manual(values=c("black", "red")) + 
  geom_point() + theme_classic()  + theme(text = element_text(size = 15), legend.position="none") +xlab('log2FC(muscle/non-muscle)') 
dev.off()
#edited x-label

#change x-scale to not be 'log'
pdf( paste0("figures/vol_plot_muscle_pseudobulked3-nolabels_notlogx_",date,".pdf" ), height=4, width=5 ) 
ggplot(data=dds.res_df %>% arrange(col), aes(x=2^log2FoldChange, y=-log10(padj), col=col)) + 
  scale_color_manual(values=c("black", "red")) + #aes(group=col) + 
  geom_point() + theme_classic()  + theme(text = element_text(size = 15), legend.position="none") +xlab('Fold Change (muscle/non-muscle)')  +
  scale_x_continuous(trans='log2', labels = scales::number_format(accuracy = 0.001))
dev.off()

