# Analysis of human mTEC bulk RNAseq
#####
library(tidyverse)
library(edgeR) 
library(RColorBrewer)
library(pheatmap)
library(pcaMethods) 
library(ggrepel)
library(scales)
library(ggthemes)
library(DESeq2)
library(gage)
library(pathview)
library(gageData)

date = "20240115"

# Import human data
norm_table <- read.delim( "Batch135_Human_adjust.gct", sep="\t", header=T, skip=2, row.names=1 )
norm_table <- norm_table[,2:ncol(norm_table)]
sample_names <- colnames(norm_table)
count_table <- read.delim( "Genes_count_table.tsv", sep="\t", header=T, row.names=1 ) #unnormalized
#QC
idx <- rowSums( norm_table > 20 ) < 2 
norm_table <- norm_table[!idx,]

# Import mouse data
norm_table_ms <- read.delim( "gene_expression_1_adjust.gct", sep="\t", header=T, skip=2, row.names=1 )
norm_table_ms <- norm_table_ms[,c("MEChi.Thy.1","MEChi.Thy.2","MEChi.Thy.3",
                                  "PDPNp.CD104p.MEClo.Thy.1", "PDPNp.CD104p.MEClo.Thy.2","PDPNp.CD104p.MEClo.Thy.3",
                                  "PDPNn.CD104n.MEClo.Thy.1","PDPNn.CD104n.MEClo.Thy.2","PDPNn.CD104n.MEClo.Thy.3")]
#QC
idx <- rowSums( norm_table_ms > 20 ) < 2 
norm_table_ms <- norm_table_ms[!idx,]

##################################
##pca analysis FUNCTION
pca_func <- function(norm_table, factor_idx, label=''){
  #run PCA
  pca_object <- pca( t(norm_table), method="svd", nPcs=5, scale="uv", center=T )
  
  #examine PC scores
  cmap <- ggthemes_data$tableau$`color-palettes`$regular$`Classic 10`$value
  pdf( paste0("figures/PCA_plot_",label,'_', date, ".pdf" ), height=4.5, width=6 )
  i=1
  temp_plot <- ggplot( as.data.frame(pca_object@scores), aes( x=pca_object@scores[,i],y=pca_object@scores[,i+1] ) ) +
    geom_point( aes( color=factor_idx ), size=4 ) +
    #geom_text_repel( label=anno_names[,1] ) + #turned point labels off
    theme_classic() + #rid of grey background
    theme( legend.position="none" ) +
    xlab( paste0("PC",i, " (",pca_object@R2[i]*100,"%)") ) +
    ylab( paste0("PC",i+1, " (",pca_object@R2[i+1]*100,"%)") ) +
    scale_color_manual(values=c("steelblue","darkgray","deeppink3") )
  #  scale_color_manual(values=cmap ) +
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1) )
  print( temp_plot )
  dev.off()
}

###################################
##### HEATMAP

heatmap_fn <- function(features,file_label,features_label=NULL,ht=13,wt=5.5){
  #######
  #features = list of genes to be plotted
  #file_label = file label extension
  #features_label = which genes to be labelled if a subset to be labelled; default NULL (plots all genes)
  #######
  scaled_norm_table <- log2( norm_table ) 
  scaled_norm_table <- scaled_norm_table - rowMeans(scaled_norm_table)
  
  factor_idx <- factor(c( rep("DN.mTEClo", times=3),rep("DP.mTEClo", times=3),rep("mTEChi", times=5) ))
  annotation_names = as.data.frame(factor_idx)
  rownames( annotation_names ) = colnames( norm_table )
  colnames(annotation_names) <- c("celltype")
  
  anno_colors <- list( celltype=c("steelblue","darkgray","deeppink3"),
                       genes=c("steelblue","darkgray","deeppink3") )
  names(anno_colors$celltype) <- unique(annotation_names$celltype)
  
  breaksList <- seq(-3,3,by=0.1) #color bar scale
  
  features_in <- features[features %in% rownames(scaled_norm_table)] 
  
  #edit to make rownames italics
  newnames <- lapply(
    rownames(scaled_norm_table[features_in,]),
    function(x) bquote(italic(.(x))))
  
  heat = pheatmap( scaled_norm_table[features_in,c(4:6,7:11,1:3)],
            cluster_rows=F, 
            cluster_cols=F, 
            show_annotation_row=F, 
            annotation_col=annotation_names,
            annotation_colors = anno_colors,
            show_rownames = T, 
            labels_row = as.expression(newnames),#edit to make rownames italics
            labels_col = as.vector(annotation_names$genotype),
            main="Human mTECs", 
            annotation_legend = T,
            breaks=breaksList, 
            border_color = NA,
            gaps_col = c(3,8),
            color=colorRampPalette(rev(brewer.pal(name="RdBu",n=11)))(length(breaksList))
  )

  pdf( paste0("figures/heatmap_",file_label,"_", date, ".pdf"), height=ht, width=wt)
  print(heat)
  dev.off()
  
  if (!is.null(features_label)){
    #edit to make rownames italics
    features_label <- lapply(
      features_label,
      function(x) bquote(italic(.(x))))
    heat <- add.flag(heat,kept.labels = features_label,repel.degree = 0)
    pdf( paste0("figures/heatmap_",file_label,"_sparselabels_", date, ".pdf"), height=13, width=5.5)
    #print(heat)
    grid.draw(heat)
    dev.off()
  }
}

########## function "add.flag", to label select rows, is from Z.Lin: https://stackoverflow.com/questions/52599180/partial-row-labels-heatmap-r
library(grid)
add.flag <- function(pheatmap, kept.labels,repel.degree) {
  # repel.degree = number within [0, 1], which controls how much space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  
  heatmap <- pheatmap$gtable
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, new.label$label, "")
  
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,grobs = new.flag,t = 4, l = 4)
  
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  return(heatmap)
}
##########

##################### identify top differential genes in human dataset and plot on a heatmap 
de_deseq2 <- function(i,j){
  #########
  # i = sample 1 name (options = "DN.mTEClo", "DP.mTEClo", "mTEChi")
  # j = other sample name (same options as i); or it can be '1vsRest', in which case it's the samples that aren't i
  #########
  idx <- substr( colnames(count_table), 1, 9 ) == i | substr( colnames(count_table), 1, 6 ) == i
  ctrl_idx <- substr( colnames(count_table), 1, 9 ) == j | substr( colnames(count_table), 1, 6 ) == j
  if (j == '1vsRest') {ctrl_idx <- !idx} #if doing sample i vs both of the other
  
  group <- factor(c(rep(1,times=sum(ctrl_idx)),rep(2,times=sum(idx))))
  coldat=DataFrame(group)
  dds <- DESeqDataSetFromMatrix(countData = cbind( count_table[,ctrl_idx], count_table[,idx]), #Deseq2 #adjusted ordering so ctrl_idx treated as control
                                colData = coldat,
                                design = ~ group) #group in coldat$group
  id_remove <- rowSums( counts(dds) > 20 ) < 2 
  dds <- dds[!id_remove,]
  dds <- DESeq(dds) 
  deseq2.res <- results(dds)
  deseq2.fc=deseq2.res$log2FoldChange
  dds.res_df <- as.data.frame(results(dds))
  return(dds.res_df) #dds)
}

topgenes_deseq2 <- function(dds.res_df){ #take all genes with a given p-value and log2FC cutoff
  upper = dds.res_df[dds.res_df$log2FoldChange>=1,] 
  lower = dds.res_df[dds.res_df$log2FoldChange<=-1,] 
  upper = rownames(upper[upper$padj<0.05,])
  lower = rownames(lower[lower$padj<0.05,])
  return(list(upper, lower))
}

### Plot heatmaps of human and mouse data, with corresponding genes 
# HUMAN
features <- c("COL17A1","ITGB4","ITGB5","PTGDR","TDGF1",
              "AIRE","CD70","CD74","CD80","CD86","HLA-DQA1",'HLA-DQB2',
              "ACTA1",'CFTR',"CHGB","CHRNG",'CKM',"MYOG","POU4F1","SNAP25",'TTN')
heatmap_fn(features,"ms-equiv",features_label=NULL,ht=5.3,wt=5.1)

# MOUSE #skip using function so can directly tailor column naming
features <- c("Col17a1","Itgb4","Itgb5","Ptgdr","Tdgf1",
              "Aire","Cd70","Cd74","Cd80","Cd86","H2-Aa",'H2-Ab1',
              "Acta1",'Cftr',"Chgb","Chrng",'Ckm',"Myog","Neurod1","Pou4f1","Snap25",'Ttn')

factor_idx <- factor(c( rep("MEChi.Thy", times=3),rep("PDPNp.CD104p.MEClo.Thy", times=3),rep("PDPNn.CD104n.MEClo.Thy", times=3) )) #ms
scaled_norm_table <- log2( norm_table_ms)
scaled_norm_table <- scaled_norm_table - rowMeans(scaled_norm_table)

annotation_names = as.data.frame(factor_idx)
rownames( annotation_names ) = colnames( norm_table_ms )
colnames(annotation_names) <- c("celltype")

anno_colors <- list( celltype=c("deeppink3","darkgray","steelblue"),
                     genes=c("deeppink3","darkgray","steelblue") )
names(anno_colors$celltype) <- unique(annotation_names$celltype)

breaksList <- seq(-3,3,by=0.1) #color bar scale

features_in <- features[features %in% rownames(scaled_norm_table)] 

#edit to make rownames italics
newnames <- lapply(
  rownames(scaled_norm_table[features_in,]),
  function(x) bquote(italic(.(x))))

heat = pheatmap( scaled_norm_table[features_in,c(4:6,1:3,7:9)],
                 cluster_rows=F,
                 cluster_cols=F,
                 show_annotation_row=F,
                 annotation_col=annotation_names,
                 annotation_colors = anno_colors,
                 show_rownames = T, 
                 labels_row = as.expression(newnames),#edit to make rownames italics
                 labels_col = as.vector(annotation_names$genotype),
                 main="Mouse mTECs",
                 annotation_legend = T,
                 breaks=breaksList,
                 border_color = NA,
                 gaps_col = c(3,6),
                 color=colorRampPalette(rev(brewer.pal(name="RdBu",n=11)))(length(breaksList))
)
pdf( paste0("figures/heatmap_mouse_", date, ".pdf"), height=5.5, width=5.5)
print(heat)
dev.off()


### Make a heatmap with differential genes between compartments, identified using DESEQ2
# DESEQ2 with different comparisons, to identify genes to plot in heatmap:
dds.res_df_DN <- de_deseq2("DN.mTEClo", "1vsRest")
dds.res_df_hi <- de_deseq2("mTEChi", "1vsRest")
dds.res_df_DP <- de_deseq2("DP.mTEClo", "1vsRest")

##plot differential genes in heatmap, label select genes
features_label = c("PTGDR","TDGF1","C3","IL13RA2","KRT14","COL17A1","KRT5","ITGB4","ITGB5","CCL21",
                   "HLA-DRA","HLA-DRB1","CD74","CD80","CD86","AIRE","CD70","ALOX5","CERKL","BSND","SNAP25","FOXI1",
                   "ACTA1","MYOG","CHGB","DEFA5","NEUROD1","POU4F1","CHRNG",'CKM','CFTR','TTN')

#compartment vs all other samples
features1 = topgenes_deseq2(dds.res_df_DN)
features2 = topgenes_deseq2(dds.res_df_hi)
features3 = topgenes_deseq2(dds.res_df_DP)
combo = unique(c(unlist(features3[1]),unlist(features2[1]),unlist(features1[1]) )) 
heatmap_fn(combo, 'topN_1vsall', features_label)

print(length(unlist(features3[1])))
print(length(unlist(features2[1])))
print(length(unlist(features1[1])))

##Plot PCA on these genes
combo = unique(c(unlist(features3[1]),unlist(features2[1]),unlist(features1[1]),unlist(features3[2]),unlist(features2[2]),unlist(features1[2]) )) #up and down genes
variable <- norm_table[combo,]
variable <- variable %>% drop_na() #remove NA rows
factor_idx <- factor( c(1,1,1,2,2,2,3,3,3,3,3) )
pca_func(variable, factor_idx,'DEG_1vsall_all')
