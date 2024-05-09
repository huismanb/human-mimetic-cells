#####
# Preprocessing of zebrafish data

#####
#start
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

date="20240219"

cmap = c("#A0CBE8","#FABFD2", "#F28E2B","#59A14F","#B07AA1",
          "#8CD17D", "#B6992D","#F1CE63", "#499894", "#86BCB6", "#E15759",
         "#FF9D9A", "#79706E", "#BAB0AC","#D37295", "#4E79A7",  "#ef6f6a", "#D4A6C8", "#9D7660", "#D7B5A6")


#####
#function to process data after removing clusters
## repeat normalization if removed clusters
norm_red <- function(seurat_table, nfeatures = 2000, res = 2, ndims=30){
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
  seurat_table <- FindNeighbors(seurat_table, dims = 1:ndims)
  seurat_table <- FindClusters(seurat_table, resolution = res)
  seurat_table <- RunUMAP(seurat_table, dims = 1:ndims, seed.use = 123)
  
  return(seurat_table)
}

add_more_metadata <- function(seurat_table_all){
  #add to metadata so can split by experiments (merge samples from multiple lanes of 1 experiment)
  seurat_table_all@meta.data$macro.ident = ''
  seurat_table_all@meta.data$macro.ident[seurat_table_all@meta.data$orig.ident == "4wpfThymus_1"] = "4wpfThymus"
  seurat_table_all@meta.data$macro.ident[seurat_table_all@meta.data$orig.ident == "4wpfThymus_3"] = "4wpfThymus"
  seurat_table_all@meta.data$macro.ident[seurat_table_all@meta.data$orig.ident == "4wpfThymus_5"] = "4wpfThymus"
  seurat_table_all@meta.data$macro.ident[seurat_table_all@meta.data$orig.ident == "4wpfThymus_6"] = "4wpfThymus"
  seurat_table_all@meta.data$macro.ident[seurat_table_all@meta.data$orig.ident == "adult_thy_liberase_1"] = "adult_thy_liberase"
  seurat_table_all@meta.data$macro.ident[seurat_table_all@meta.data$orig.ident == "adult_thy_liberase_2"] = "adult_thy_liberase"
  seurat_table_all@meta.data$macro.ident[seurat_table_all@meta.data$orig.ident == "adult_thy_untreated"] = "adult_thy_untreated"
  seurat_table_all@meta.data$macro.ident[seurat_table_all@meta.data$orig.ident == "Fish_22_Thymus_"] = "Fish_22_Thymus"
  seurat_table_all@meta.data$macro.ident[seurat_table_all@meta.data$orig.ident == "Fish_24_Thymus_"] = "Fish_24_Thymus"
  seurat_table_all@meta.data$macro.ident[seurat_table_all@meta.data$orig.ident == "lck-eGFP_Thymus_1_"] = "lck-eGFP_Thymus_1"
  seurat_table_all@meta.data$macro.ident[seurat_table_all@meta.data$orig.ident == "lck-eGFP_Thymus_2_"] = "lck-eGFP_Thymus_2"
  ##
  #add to metadata on age
  seurat_table_all@meta.data$age.ident = ''
  seurat_table_all@meta.data$age.ident[seurat_table_all@meta.data$macro.ident == "4wpfThymus"] = "4wpf"
  seurat_table_all@meta.data$age.ident[seurat_table_all@meta.data$macro.ident == "adult_thy_liberase"] = "adult"
  seurat_table_all@meta.data$age.ident[seurat_table_all@meta.data$macro.ident == "adult_thy_untreated"] = "adult"
  seurat_table_all@meta.data$age.ident[seurat_table_all@meta.data$macro.ident == "Fish_22_Thymus"] = "adult"
  seurat_table_all@meta.data$age.ident[seurat_table_all@meta.data$macro.ident == "Fish_24_Thymus"] = "adult"
  seurat_table_all@meta.data$age.ident[seurat_table_all@meta.data$macro.ident == "lck-eGFP_Thymus_1"] = "adult"
  seurat_table_all@meta.data$age.ident[seurat_table_all@meta.data$macro.ident == "lck-eGFP_Thymus_2"] = "adult"
  ##
  #add to metadata on sort
  seurat_table_all@meta.data$sort.ident = ''
  seurat_table_all@meta.data$sort.ident[seurat_table_all@meta.data$macro.ident == "4wpfThymus"] = "Live"
  seurat_table_all@meta.data$sort.ident[seurat_table_all@meta.data$macro.ident == "adult_thy_liberase"] = "non T or B"
  seurat_table_all@meta.data$sort.ident[seurat_table_all@meta.data$macro.ident == "adult_thy_untreated"] = "non T or B"
  seurat_table_all@meta.data$sort.ident[seurat_table_all@meta.data$macro.ident == "Fish_22_Thymus"] = "Live"
  seurat_table_all@meta.data$sort.ident[seurat_table_all@meta.data$macro.ident == "Fish_24_Thymus"] = "Live"
  seurat_table_all@meta.data$sort.ident[seurat_table_all@meta.data$macro.ident == "lck-eGFP_Thymus_1"] = "Live"
  seurat_table_all@meta.data$sort.ident[seurat_table_all@meta.data$macro.ident == "lck-eGFP_Thymus_2"] = "Live"
  ##
  return(seurat_table_all)
}

#####
## Import the trimmed data with major contaminating clusters removed to proceed with (after cut 1)
seurat_table <- readRDS("Rdata/zebrafish_mec_seurat_table_20231017_afterfilter1-wDoubletFinder.rds")

# With trimmed data in rds file - after cut1:
pdf( paste0("figures/umap_aftercut1_",date,".pdf" ), height=5, width=6 )
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1) + theme(legend.position="none")
dev.off()

pdf( paste0("figures/features_aftercut1_",date,".pdf" ), height=8, width=10 )
FeaturePlot(seurat_table, features = c('epcam',"tcf7","cd4-1",'cd8a','lck','EGFP','cd79b','igl1c3','pax5'), pt.size=0.7, order=T, ncol=3 )
dev.off()

## Cut 2:
keep <- seurat_table@active.ident!=3 & seurat_table@active.ident!=16 & seurat_table@active.ident!=22 & seurat_table@active.ident!=20 & #T
  seurat_table@active.ident!=0 & seurat_table@active.ident!=2 #B
seurat_table <- seurat_table[,keep]

#repeat normalization and reduction after removing contaminating cells
seurat_table <- norm_red(seurat_table)

## After cut 2:
pdf( paste0("figures/umap_aftercut2_",date,".pdf" ), height=5, width=6 )
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1) + theme(legend.position="none")
dev.off()

pdf( paste0("figures/features_aftercut2_",date,".pdf" ), height=8, width=10 )
FeaturePlot(seurat_table, features = c('epcam',"tcf7","cd4-1",'cd8a','lck','EGFP','cd79b','igl1c3','pax5'), pt.size=0.7, order=T, ncol=3 )
dev.off()

## Cut 3:
keep <- seurat_table@active.ident!=11 & seurat_table@active.ident!=13 & seurat_table@active.ident!=30
seurat_table <- seurat_table[,keep]

#repeat normalization and reduction after removing contaminating cells
seurat_table <- norm_red(seurat_table)

# After cut 3:
pdf( paste0("figures/umap_aftercut3_",date,".pdf" ), height=5, width=6 )
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1) + theme(legend.position="none")
dev.off()

pdf( paste0("figures/features_aftercut3_",date,".pdf" ), height=8, width=10 )
FeaturePlot(seurat_table, features = c('epcam',"tcf7","cd4-1",'cd8a','lck','EGFP','cd79b','igl1c3','pax5'), pt.size=0.7, order=T, ncol=3 )
dev.off()

## Cut 4: 
keep <- seurat_table@active.ident!=1 & seurat_table@active.ident!=2 & seurat_table@active.ident!=9 & 
  seurat_table@active.ident!=16 & seurat_table@active.ident!=29 & seurat_table@active.ident!=4 &
  seurat_table@active.ident!=25 & seurat_table@active.ident!=5
seurat_table <- seurat_table[,keep]

#repeat normalization and reduction after removing contaminating cells
seurat_table <- norm_red(seurat_table)

# After cut 4:
pdf( paste0("figures/umap_aftercut4_",date,".pdf" ), height=5, width=6 )
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1) + theme(legend.position="none")
dev.off()

pdf( paste0("figures/features_aftercut4_",date,".pdf" ), height=8, width=10 )
FeaturePlot(seurat_table, features = c('epcam',"tcf7","cd4-1",'cd8a','lck','EGFP','cd79b','igl1c3','pax5'), pt.size=0.7, order=T, ncol=3 )
dev.off()

## Cut 5
keep <- seurat_table@active.ident!=18 #lck/T cell markers
seurat_table <- seurat_table[,keep]
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1) + theme(legend.position="none")

#repeat normalization and reduction after removing contaminating cells
seurat_table <- norm_red(seurat_table)

# After cut 5:
pdf( paste0("figures/umap_",date,".pdf" ), height=5, width=6 )
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1) + theme(legend.position="none")
dev.off()

pdf( paste0("figures/features_aftercut5_",date,".pdf" ), height=8, width=10 )
FeaturePlot(seurat_table, features = c('epcam',"tcf7","cd4-1",'cd8a','lck','EGFP','cd79b','igl1c3','pax5'), pt.size=0.7, order=T, ncol=3 )
dev.off()


###################################################################################################
## Naming clusters
seurat_table <- RenameIdents( object=seurat_table,
                              '22'='cTEC',
                              '1'='aire-stage',
                              '6'='Immature',
                              '4'='Immature',
                              '12'='Immature',
                              '0' = 'Immature',
                              '18'='Muscle',
                              '2'="Periderm", 
                              '13'='Periderm',
                              '11'='Ionocyte (HR rich)', 
                              '10'='Ionocyte (NaR rich)', 
                              '23'='dlx3b+', #dlx3b+ structural
                              '16'='Ear sensory hair cell', 
                              '17'='Ciliated (foxj1a+)', 
                              '7'="Ciliated olfactory sensory neurons", 
                              '8'= 'Neuromast', 
                              '21'='Ear non-sensory cell', 
                              '19'="Microvillous olfactory sensory neurons", 
                              '3'= 'Neuro',
                              '15'= 'Tuft1 (neurosensory)', 
                              '5'='Tuft2 (immune)', 
                              '20'='Tuft (gut-type)', 
                              '9'= 'Metaphocytes', 
                              '14'='spi1a+, spi1b+') #macrophage-like

pdf( paste0("figures/umap_relabeled_",date,".pdf" ), height=10, width=12 ) #height=10, width=12 )
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1, cols=cmap ) + theme(legend.position="none")
dev.off()

pdf( paste0("figures/umap_nolabels_",date,".pdf" ), height=5, width=6 )
DimPlot(seurat_table, reduction="umap", label=F, repel=T, pt.size=1, cols=cmap) + theme(legend.position="none")
dev.off()

########
seurat_table = add_more_metadata(seurat_table)

pdf( paste0("figures/umap_split-macroident_",date,".pdf" ), height=5, width=24 )
DimPlot(seurat_table, reduction="umap", label=F, repel=T, pt.size=1,split.by = 'macro.ident', cols=cmap ) + theme(legend.position="none")
dev.off()

pdf( paste0("figures/umap_split-age_",date,".pdf" ), height=5, width=6 )
DimPlot(seurat_table, reduction="umap", label=F, repel=T, pt.size=1, split.by = 'age.ident', cols=cmap  ) + theme(legend.position="none")
dev.off()

pdf( paste0("figures/umap_split-sort_",date,".pdf" ), height=5, width=6 )
DimPlot(seurat_table, reduction="umap", label=F, repel=T, pt.size=1, split.by = 'sort.ident', cols=cmap  ) + theme(legend.position="none")
dev.off()

#########
   
#Cluster TFs
pdf( paste0("figures/dotplot_tf_",date,".pdf"), height=6, width=15.5 )
DotPlot(seurat_table, features=(c("foxn1","tp63",'aire','gata2a', "gata3", 
                                  'myog','myod1',"klf2b","grhl3",
                                  'grhl1','foxi3b','foxi3a', 
                                  'dlx3b', 'six4b','znf750', 
                                  'pax2a','otofb', 'foxj1a', 'foxj1b', 
                                  'atoh1a','drgx',"sox2", "prox1a", 'sox21a',
                                  'sox10',"msx3",'pax6a','emx2',
                                  'neurod1','zeb2a', "lhx9","myt1la","npas4a","pou2f2a", 
                                  "emx3","pitx3","hoxb2a","hoxb3a","nfat5b", 
                                  "pou2f3",'sox8b',"gata1b",'spi2',
                                  'spic','spi1a', 'spi1b')),
        cols=c("lightgray","#3F007D"), dot.min = 0.05, dot.scale=8 ) + theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust=1)) +scale_y_discrete(limits=rev) +#reverse y-direction
  theme(axis.text.x = element_text(face = "italic"), axis.title.x=element_blank(), axis.title.y=element_blank() )
dev.off()

############################
#with new additions:
heatmap_features_narrow = c('foxn1','psmb11b','psmb9a', 'ccl25a',
                            'aire','psmd10', 'col1a2', 'cxcl8b.1','ccl19b',
                            'myog','ttn.1','ttn.2',
                            'anxa1c','abcb5', "krt17", 
                            'atp6v1aa', 'slc9a3.2', 'ca15a', 
                            'kcnj1a.1','trpv6','igfbp5a','slc8a1b','dlx3b', 
                            'atoh1a', 'otofb','foxj1a', 'foxj1a','capsla', 
                            'or132-5', 'or111-11', 'ompb', 
                            'atoh1a','ntm','doc2b','capga', 'stm', 'otomp', 'otog', 'bricd5',  
                            'emx2','trpc2b', 'olfcd3', 'neurod1','nrg3b', 'pou2f3', 
                            "gng13a","trpm5","dclk1a", 'adgrg11','tnfrsf11a', 
                            'gata1b', 'sox8b','aldh1a3','vil1',
                            'spic','il22ra2','ptprc','mpeg1.1', 
                            'spi1a','spi1b','mpeg1.1') 


breakslist = seq(1,3,by=0.01) 
pdf( paste0("figures/markers_heatmap_selected_purple_narrowlist_",date,".pdf" ), height=10, width=10 )
DoHeatmap(subset(seurat_table, downsample=50), features = unique(heatmap_features_narrow), angle=90, disp.min=0.5, size=4, group.colors=cmap ) + theme(axis.text.y.left = element_text(face = "italic"))+ #NoLegend() + #size changes cluster label size
  scale_fill_gradientn(colors = colorRampPalette((brewer.pal(n=9, name="Purples")))(length(breakslist)))
dev.off()


##############################
# Subcluster subsets of clusters
cmap_pou = c( "#4E79A7",   "#D4A6C8", "#ef6f6a" )
cmap_iono = c("#B6992D","#8CD17D")

## subcluster ionocyte cells 
keep <- seurat_table@active.ident=='Ionocyte (HR rich)' | seurat_table@active.ident=='Ionocyte (NaR rich)'
seurat_table_iono <- seurat_table[,keep]
DimPlot(seurat_table_iono, reduction="umap", label=T, repel=T, pt.size=1 ) + theme(legend.position="none")

seurat_table_iono <- norm_red(seurat_table_iono, 2000, 0.5)

pdf( paste0("figures/umap_ionocyte_",date,".pdf" ), height=5, width=5 )
DimPlot(seurat_table_iono, reduction="umap", label=F, repel=T, pt.size=1) #+ theme(legend.position="none")
dev.off()

seurat_table_iono <- RenameIdents( object=seurat_table_iono,
                                   '0'="NaR-rich",
                                   '1'="NaR-rich",
                                   '2'='HR-rich')

pdf( paste0("figures/umap_ionocyte_merged_",date,".pdf" ), height=4, width=4 )
DimPlot(seurat_table_iono, reduction="umap", label=T, repel=T, pt.size=1, cols = cmap_iono ) + theme(legend.position="none")
dev.off()

pdf( paste0("figures/umap_ionocyte_merged_unlabeled_",date,".pdf" ), height=2.75, width=2.75 )
DimPlot(seurat_table_iono, reduction="umap", label=F, repel=T, pt.size=1, cols=cmap_iono )  & NoAxes() + theme(legend.position="none")
dev.off()

pdf( paste0("figures/clusterfeatures_ionocyte_",date,".pdf" ), height=5.5, width=7.5 )
FeaturePlot(seurat_table_iono, features = c('trpv6', #NaK ionocyte
                                       'slc12a10.2', #NCC
                                       'kcnj1a.1', # KS
                                       'atp6v1aa', 'slc9a3.2', 'ca15a'), #HR  
            pt.size=1, order=T, ncol=3, cols=c("lightgray","#3F007D")) & NoAxes() & theme( plot.title = element_text( face = "bold.italic") ) &NoLegend()
dev.off()

##############################################
## subcluster pou2f3 cells 
keep <- seurat_table@active.ident=='Tuft1 (neurosensory)' | seurat_table@active.ident=='Tuft2 (immune)' | seurat_table@active.ident=='Tuft (gut-type)'
seurat_table_pou <- seurat_table[,keep]
DimPlot(seurat_table_pou, reduction="umap", label=T, repel=T, pt.size=1 ) + theme(legend.position="none")

seurat_table_pou <- norm_red(seurat_table_pou, 2000, 1)

pdf( paste0("figures/umap_pou2f3_",date,".pdf" ), height=5, width=5 )
DimPlot(seurat_table_pou, reduction="umap", label=F, repel=T, pt.size=1, cols=cmap_pou ) #+ theme(legend.position="none")
dev.off()

pdf( paste0("figures/features_pou2f3_",date,".pdf" ), height=5.5, width=12.5 )
FeaturePlot(seurat_table_pou, features = c('alox5a','alox12','il17rc',
                                           'aldh1a3','vil1',
                                           'gnb3a','calb2a','gng13a','trpm5','dclk1a'), pt.size=1, order=T, ncol=5, cols=c("lightgray","#3F007D")) & NoAxes() & theme( plot.title = element_text( face = "bold.italic") ) &NoLegend()
dev.off()

seurat_table_pou <- RenameIdents( object=seurat_table_pou,
                                   '1'="Neurosensory",
                                   '2'='Gut-type',
                                   '0'='Immune-type')

pdf( paste0("figures/umap_pou2f3_relabeled_",date,".pdf" ), height=5, width=5 )
DimPlot(seurat_table_pou, reduction="umap", label=T, repel=T, pt.size=1 ) + theme(legend.position="none")
dev.off()

pdf( paste0("figures/umap_pou2f3_unlabeled_",date,".pdf" ), height=2.75, width=2.75 )
DimPlot(seurat_table_pou, reduction="umap", label=F, repel=T, pt.size=1, cols=cmap_pou )  & NoAxes() + theme(legend.position="none")
dev.off()


##############################################
## subcluster muscle cells 
keep <- seurat_table@active.ident=='Muscle'
seurat_table_muscle <- seurat_table[,keep]
DimPlot(seurat_table_muscle, reduction="umap", label=T, repel=T, pt.size=1 ) + theme(legend.position="none")

seurat_table_muscle <- norm_red(seurat_table_muscle, nfeatures = 2000, res = 1, ndims=50) 

pdf( paste0("figures/umap_muscle_reclustered_",date,".pdf" ), height=5, width=5 )
DimPlot(seurat_table_muscle, reduction="umap", label=T, repel=T, pt.size=1 ) + theme(legend.position="none")
dev.off()

pdf( paste0("figures/umap_muscle_unlabeled-greens_",date,".pdf" ), height=2.75, width=2.75 )
DimPlot(seurat_table_muscle, reduction="umap", label=F, repel=T, pt.size=1, cols=c("#85BF7D","#3D6E36") )  & NoAxes() + theme(legend.position="none")
dev.off()

pdf( paste0("figures/features_muscle_subset3_",date,".pdf" ), height=5.5, width=7.5 )
FeaturePlot(seurat_table_muscle, features = c('epcam','krt5','ccl19b',
                                              'myod1','ttn.1','mylpfa'
), pt.size=1, order=T, cols=c("lightgray","#3F007D"), ncol=3) & NoAxes() & theme( plot.title = element_text( face = "bold.italic") ) &NoLegend()
dev.off()

####
# plotting species side by side
## import human data
human_seurat <- readRDS( "human_seurat.rds" )

## import mouse data
mouse_seurat <- readRDS("mouse_seurat.rds")

#remove non-mimetic cells (eg get rid of CCL19+ mTECs):
keep <- human_seurat@active.ident!='Immature' & human_seurat@active.ident!='Cycling' & human_seurat@active.ident!="Aire-stage"
human_seurat_mimetic <- human_seurat[,keep]

#remove non-mimetic cells:
keep <- mouse_seurat@active.ident!='Aire-stage' & mouse_seurat@active.ident!='TA MEC' & mouse_seurat@active.ident!='perinatal cTEC'  & 
  mouse_seurat@active.ident!='adult cTEC' & mouse_seurat@active.ident!='Immature MEC' 
mouse_seurat_mimetic <- mouse_seurat[,keep]

#remove non-mimetic cells:
keep <- seurat_table@active.ident!='aire-stage' & seurat_table@active.ident!='Immature' & seurat_table@active.ident!='cTEC'
seurat_table_mimetic <- seurat_table[,keep]

#order clusters
seurat_table_mimetic@active.ident <- factor(seurat_table_mimetic@active.ident,
                                    levels=c("Periderm","dlx3b+",
                                             "Ionocyte (HR rich)","Ionocyte (NaR rich)", "Muscle",
                                             "Ciliated (foxj1a+)",
                                             "Ear sensory hair cell","Ciliated olfactory sensory neurons",
                                             "Neuromast","Ear non-sensory cell","Microvillous olfactory sensory neurons","Neuro",
                                             "Tuft1 (neurosensory)","Tuft2 (immune)","Tuft (gut-type)",
                                             "Metaphocytes","spi1a+, spi1b+"))

print(levels(seurat_table_mimetic@active.ident))

#order clusters
print(levels(human_seurat_mimetic@active.ident))
human_seurat_mimetic@active.ident <- factor(human_seurat_mimetic@active.ident,
                                    levels=c("Keratinocyte",
                                             "HMX2+ ionocyte", "Ionocyte", 
                                             "Early muscle","Late muscle","lncRNA-enriched muscle",
                                             "Cochlear hair cell","DRGX+ sensory neuro","SHOX2+ neuro","NKX6-2+ neuro", "FEZF2+ neuro","CUX2+ neuro",
                                             "Tuft",
                                             "Intermediate"))

#order clusters
mouse_seurat_mimetic@active.ident <- factor(mouse_seurat_mimetic@active.ident,
                                    levels=c("Skin, basal","Skin, keratinized","Lung, basal",
                                             "Ionocyte","Muscle",
                                             "Neuroendocrine","Ciliated",
                                             "Tuft1","Tuft2","Mcell", "Gut/Liver","Goblet",
                                             "Ptf1a+ ductal"))

####
p_zf_n = DotPlot(seurat_table_mimetic, features=(c('grhl1','dlx3b','foxi3b','myog','myod1','foxj1a','atoh1a','pou4f1','neurod1','ehf','pou2f3','sox8b','spi1b')),
                 cols=c("lightgray","#3F007D"), dot.min = 0.05, dot.scale=8 ) + theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust=1)) +
  theme(axis.text.x = element_text(face = "italic"), axis.title.x=element_blank(), axis.title.y=element_blank() )+
  scale_y_discrete(limits=rev) #reverse y-direction

p_hu_n = DotPlot(human_seurat_mimetic, features=(c('GRHL1','DLX3','FOXI1','MYOG','MYOD1','FOXJ1','ATOH1','POU4F1','NEUROD1','EHF','POU2F3','SOX8','SPI1')),
                 cols=c("lightgray","#3F007D"), dot.min = 0.05, dot.scale=8 ) + theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust=1)) +
  theme(axis.text.x = element_text(face = "italic"), axis.title.x=element_blank(), axis.title.y=element_blank() )+
  scale_y_discrete(limits=rev) #reverse y-direction #+ theme(plot.margin = unit(c(2,2,2,2), 'cm' )) #add space around margin so can be similar btwn plots

p_ms_n = DotPlot(mouse_seurat_mimetic, features=(c('Grhl1','Dlx3','Foxi1','Myog','Myod1','Foxj1','Atoh1','Pou4f1','Neurod1','Ehf','Pou2f3','Sox8','Spi1')),
                 cols=c("lightgray","#3F007D"), dot.min = 0.05, dot.scale=8 ) + theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust=1)) +
  theme(axis.text.x = element_text(face = "italic"), axis.title.x=element_blank(), axis.title.y=element_blank() )+
  scale_y_discrete(limits=rev) #reverse y-direction

pdf( paste0("figures/dotplot_sharedTFs_3species_narrowed_mimeticonly_",date,".pdf"), height=16, width=11 )
plot_grid(p_hu_n, p_ms_n, p_zf_n, ncol=1, align="v")
dev.off()

