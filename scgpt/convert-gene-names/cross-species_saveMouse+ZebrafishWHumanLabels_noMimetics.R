#####
# Cross-species analysis - convert gene names
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
library(reticulate)
library(sceasy)
##
plan("sequential")

######################################## Import labelled and processed seurat objects
## import mouse data
mouse_seurat <- readRDS("mouse_seurat.rds")

#remove non-mimetic cells
keep <- mouse_seurat@active.ident!='perinatal cTEC' & mouse_seurat@active.ident!='adult cTEC' &
  mouse_seurat@active.ident!='TA MEC' & mouse_seurat@active.ident!='Immature MEC' & mouse_seurat@active.ident!='Aire-stage'
mouse_seurat <- mouse_seurat[,keep]

## import zebrafish data:
zebra_seurat <- readRDS("zebrafish_seurat.rds")

#remove non-mimetic cells
keep <- zebra_seurat@active.ident!='aire-stage' & zebra_seurat@active.ident!='Immature' & zebra_seurat@active.ident!='cTEC' 
zebra_seurat <- zebra_seurat[,keep]

##################################
#Functions for renaming genes and filtering so can integrate across species

#filter and rename genes for species1 (species being renamed)
#adapted from https://github.com/Papatheodorou-Group/BENGAL/blob/main/bin/concat_by_homology_multiple_species_by_gene_id.R
filter_rename_species1 <- function(dat,one2one,species0,species1){
  
  one2one_now = one2one[ one2one[,paste0(species1,'_homolog_associated_gene_name')] %in% dat@Dimnames[[1]] ,] %>% #keep genes in one2one that are in dat
    distinct( get(paste0(species0, "_homolog_associated_gene_name")), .keep_all = TRUE) 
  
  dat = dat[ # filter out genes that don't have species0 homolog
    dat@Dimnames[[1]] %in% one2one_now[[paste0(species1, "_homolog_associated_gene_name")]],  #filter genes
  ]
  
  one2one_now = one2one_now[match( dat@Dimnames[[1]] , one2one_now[[paste0(species1, '_homolog_associated_gene_name')]]), ] #keep genes in one2one_now that match with dat (+match/order them)
  dat@Dimnames[[1]] = one2one_now[[paste0(species0,'_homolog_associated_gene_name')]] #replace genename with species0 version
  return(dat)
}

#filter genes for species0 (reference species; species1 will take gene names of species0)
filter_species0 <- function(dat, one2one, species0){
  dat = dat[ dat@Dimnames[[1]] %in% one2one[[paste0(species0,"_homolog_associated_gene_name")]], ]  
  return(dat)
}
##################################
# Homology table generated in python
one2one_zf2human = read.csv('/Users/brooke/Dropbox (HMS)/CBDM_Lab/Computational/Mouse-human gene conversion/homology_tbl_hsapiens_drerio-20231117.csv')
one2one_mouse2human = read.csv('/Users/brooke/Dropbox (HMS)/CBDM_Lab/Computational/Mouse-human gene conversion/homology_tbl_hsapiens_mmusculus-20231117.csv')
one2one_zf2mouse = read.csv('/Users/brooke/Dropbox (HMS)/CBDM_Lab/Computational/Mouse-human gene conversion/homology_tbl_mmusculus_drerio-20231117.csv')

##########################
## rename mouse to human genenames
species0 = 'hsapiens'
species1 = 'mmusculus'
one2one = one2one_mouse2human

dat_species1 = mouse_seurat[["RNA"]]@counts
dat_species1_out = filter_rename_species1(dat_species1, one2one, species0, species1)

saveRDS(dat_species1_out, 'Rdata/meclo_seurat_table_mouse-humanlabels_mimeticsonly-20240306.rds')

##########################
## rename zebrafish to human genenames
species0 = 'hsapiens'
species1 = 'drerio'
one2one = one2one_zf2human

dat_species1 = zebra_seurat[["RNA"]]@counts
dat_species1_out = filter_rename_species1(dat_species1, one2one, species0, species1)

saveRDS(dat_species1_out, 'Rdata/zebrafish_mec_aftercut5-humanlabels_mimeticsonly-20240306.rds')
