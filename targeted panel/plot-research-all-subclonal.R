source('../../analysis/twist-panel/resources.R')
library('dplyr')
library("DESeq2")
library('Rsubread')
library('limma')
library('edgeR')
library('reshape2')
library('pheatmap')
library('ggplot2')
library('ggrepel')
library('svglite')
library('RColorBrewer')
library('tidyr')
library('GSEABase')
library('pheatmap')
library('tidyverse')
library(readxl)
library(ComplexHeatmap)
setwd('~/Dropbox/analysis/pdx-panel/')
#short <- c('mutid','ID','gene','name','Consequence','PMCFREQ','PMCAD','PMCDP','CHROM','POS','COSMIC_ids','HGVSp','HGVSc',
#           'PolyPhen','SIFT','CADD_PHRED','gnomAD_AF','CHROM','POS','REF','ALT')

ready <- read_xlsx('~/Dropbox/Monash Staff PDX/All PDX Info/Targeted_Panel_v2/snapshot/pdx-samples-rename-for-thesis-Nov-2020-snapshot.xlsx') %>% filter(research_ready!='X') ##X means remove all - not ready for publication
#order <- read_csv('../pdx-sample-master-list-2020-tree-order.csv') %>% filter(research_ready!='X')
#order <- order[order$name %in% ready$name,]$name
#ready <- ready[ready$name %in% order,]

#info <- read.table('',sep=",",stringsAsFactors = F,header=T)
out1 <- read.table('../twist-panel/cnv-out-2020-09-16.tsv',sep="\t",stringsAsFactors = F, header=T)
out3 <- read.table('../garvan-panel/cnv-out-garvan-2020-09-03.tsv',sep="\t",stringsAsFactors = F, header=T)

out2 <- read_xlsx('all-PDX-variants-compact-2020-09-16-curate.xlsx') %>% filter(curate=='Y' | PMCFREQ2<=0.75)
out2 <- filter(out2,PMCFREQ2>0.75 | as.numeric(gnomAD_AF)==0 | gnomAD_AF==".")
print("looking at all variants")

out1 <- rbind(out1,out3)

gtmp <- read.table('~/Dropbox/Monash Staff PDX/All PDX Info/Targeted_Panel_v2/data/curated-PCA-genes.csv',header=T)
ar <- read.table('~/Dropbox/Monash Staff PDX/All PDX Info/Targeted_Panel_v2/data/prad1000_CNA_Genes.txt',header=T,sep="\t")
names(ar)[8] <- 'OncoKB'

### FILTERS
MIN_FREQ <- 0.4 #minimum frequency to include in plot
THR_GAIN <- 2.8 #theshold for curated copy gain
THR_LOSS <- 1.2 #1 #threshold for curate copy loss
DB_FREQ <- 1 #minimum frequency in Armenia to count as recurrent copy number
MIN_GENE <- 3 #gene must be recurrently altered this many times to be included in plot if in known genes
MIN_GENE2 <- 3 #gene must be recurrently altered this many times to be included in plot if in known genes

out1 <- filter(out1,name %in% ready$name)
out2 <- filter(out2,name %in% ready$name)

out2 <- filter(out2,PMCFREQ2>MIN_FREQ | (curate=='Y' & mutid=='17-7574029-CG-C-TP53'))
out2 <- out2 %>% dplyr::rename(gene=SYMBOL)
snv.genes <- unique(out2$gene)

out1 <- merge(out1,ready,by="name")
out1 <- filter(out1,gene %in% gtmp$gene | gene %in% snv.genes | cnv=='gain-high',copy>=THR_GAIN | copy < THR_LOSS)
#out1 <- filter(out1,copy>=THR_GAIN | copy < THR_LOSS)

ar <- filter(ar,OncoKB=='Yes')
ar$fq <- as.numeric(gsub('<','',gsub("%","",as.character(ar$Freq))))
ar.gain <- ar %>% filter(fq>DB_FREQ,CNA=='AMP') %>% dplyr::distinct(Gene,CNA,fq)
ar.del <- ar %>% filter(fq>DB_FREQ,CNA=='HOMDEL') %>% dplyr::distinct(Gene,CNA,fq)
out1.amp <- out1[grepl('gain',out1$cnv),] %>% filter(!gene %in% ar.del$Gene | gene %in% ar.gain$Gene)
out1.loss <- out1[grepl('loss',out1$cnv),] %>% filter(!gene %in% ar.gain$Gene | gene %in% ar.del$Gene)
out1.homdel <- out1[grepl('homdel',out1$cnv),]
out1 <- rbind(rbind(out1.amp,out1.loss),out1.homdel)

#out2 <- out2 %>% dplyr::rename(gene=SYMBOL)

out1 <- out1 %>% dplyr::select(id=name,gene,type=cnv) 
out2 <- out2 %>% dplyr::select(id=name,gene,type=Consequence) 

tmp <- rbind(out2,out1)
#tmp <- filter(tmp,!gene %in% singeltons$gene |
#                 (gene %in% reoccuring$gene))

tmp <- filter(tmp,type!='neutral')
tmp$id <- as.character(tmp$id)
tmp$gene <- as.character(tmp$gene)
tmp$type <- as.character(tmp$type)
tmp[grepl('^loss$',tmp$type),]$type <- 'LOSS'
tmp[grepl('^gain-high',tmp$type),]$type <- 'AMP;AMPHIGH'
tmp[grepl('^gain',tmp$type),]$type <- 'AMP'
tmp[grepl('homdel',tmp$type),]$type <- 'HOMDEL'
tmp[grepl('frameshift',tmp$type),]$type <- 'STP'
tmp[grepl('stop',tmp$type),]$type <- 'STP'
tmp[grepl('start_lost',tmp$type),]$type <- 'STP'
tmp[grepl('missense',tmp$type),]$type <- 'MUT'
tmp[grepl('splice_acceptor_variant',tmp$type),]$type <- 'MUT'
tmp[grepl('splice_donor_variant',tmp$type),]$type <- 'MUT'
tmp[grepl('splice_region_variant',tmp$type),]$type <- 'MUT'
tmp[grepl("inframe_insertion",tmp$type),]$type <- 'MUT'
tmp[grepl("inframe_deletion",tmp$type),]$type <- 'MUT'

#manually add 294 germline
#manually add germline
tmp <- rbind(tmp,c('CDG-031-294R-PDX','BRCA2','germline'))
#manually add 27-2 finding
tmp <- rbind(tmp,c('44-CA27LN2-Cx5','AR','LOSS'))

print('leftovers')
print(unique(tmp[!tmp$type %in% c('AMP','AMP;AMPHIGH','HOMDEL','LOSS','STP','germline'),]$type))

#count genes to get subset by frequency

tmp2 <- merge(tmp,phen,by.x="id",by.y="name")
tmp2 <- tmp2 %>% dplyr::distinct(gene,pdx) %>% dplyr::count(gene) %>% arrange(desc(n))
tmp2$perc <- tmp2$n/30*100

#tmp <- filter(tmp,type!='LOSS') #no hemizygous
data <- as.data.frame(dcast(tmp,gene~id,value.var="type"))
data$gene <- as.character(data$gene)
data[is.na(data)] <- ''
row.names(data) <- data$gene
data <- dplyr::select(data,-gene)
data[!is.na(data)] <- ''
for (i in 1:nrow(tmp)) {
  id <- tmp[i,]$id
  g <- tmp[i,]$gene
  type <- tmp[i,]$type
  if (data[g,id]!='') {data[g,id] <- paste0(data[g,id],type,';')}
  if (data[g,id]=='') {data[g,id] <- paste0(type,';')}
}

alter_fun_list = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "red", col = NA))
  },
  AMPHIGH = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.3, gp = gpar(fill = "orange", col = NA))
  },
  LOSS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "steelblue1", col = NA))
  },
  HOMDEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "blue", col = NA))
  },
  STP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.35, gp = gpar(fill = "#444444", col = NA))
  },
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.3, gp = gpar(fill = "#008000", col = NA))
  },
  germline = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.25, gp = gpar(fill = "purple", col = NA))
  }
)
col = c("AMP" = "firebrick2", "AMPHIGH" = "orange","HOMDEL" = "blue", "LOSS" = "steelblue1", "STP"="#444444","MUT" = "#008000", "germline"='purple')
#col = c("AMP" = "red", "HOMDEL" = "blue", "MUT" = "#008000")

hmap <-  list(title = "Alterations", at = c("AMP","AMPHIGH","HOMDEL","LOSS","STP","MUT","germline"), 
              labels = c("Amplification","Amplification 8+","Deep deletion","Partial Loss","Stop Gain / Frameshift","Mutation","Germline"))

ndx <- read_xlsx('~/Dropbox/Monash Staff PDX/All PDX Info/Targeted_Panel_v2/snapshot/pdx-samples-rename-for-thesis-Nov-2020-snapshot.xlsx')
map.ndx <- setNames(ndx$name3,ndx$name)
names(data) <- map.ndx[names(data)]
order <- ndx[ndx$name3 %in% names(data),]$name3
###Plot all by pathway

gtmp <- read.table('~/Dropbox/analysis/twist-panel/curated-PCA-genes.csv',header=T)
twist <- read.table('/Users/bakshiandrew/Dropbox/Monash\ Staff\ PDX/All\ PDX\ Info/Targeted_Panel_v2/gene-lists/twist-genelist-2020.txt',stringsAsFactors = F,header=T)

path_ar <- c('AR','NCOR1','NCOR2','FOXA1','SPOP','SPEN') #AR pathway - remove SPEN?
path_dna <- c('BRCA2','ATM','BRCA1','FANCA','RAD51C','RAD51','MLH1','MSH2') #DNA repair
path_dna2 <- c('ATM', 'BRCA1', 'BRIP1', 'CHEK2', 'EPCAM', 'HOXB13', 'MLH1', 'MSH2', 'MSH6', 'NBN', 'PALB2', 'PMS2', 'RAD51C', 'RAD51D', 'TP53','FANCA','CDK12')
path_wnt <- c('APC','CTNNB1','RNF43','AXIN1','AXIN2','ZNRF3','RSPO2') #WNT pathway
path_pi3k <- c('PTEN','PIK3CA','PIK3CB','PIK3R1','AKT1')
path_ubi <- c('SF3B1','U2AF1','GEMIN5','TCERG1','PRPF8') #Ubiquitin and Splicing pathways

#lookup names in order

#path_dna2[path_dna2 %in% twist$gene]
#path_dna2[!path_dna2 %in% twist$gene]

hmap <-  list(title = "Alterations", at = c("AMP","AMPHIGH","HOMDEL","LOSS","STP","MUT","germline"), 
              labels = c("Amplification","Amplification 8+","Deep deletion","Partial Loss","Truncation/Frameshift","Mutation","Germline"))

#ndx <- read_xlsx('../pdx-sample-master-list-2020-tree-order.xlsx') %>% filter(research_ready!='X'
#map.ndx <- setNames(ndx$name2,ndx$name)
#names(data) <- map.ndx[names(data)]
#

co_dna <- c('TP53','PTEN','MLH1','MSH2') #DNA repair

pdf(paste0('plots/oncoprint-co-dna-',MIN_FREQ,'-',Sys.Date(),'.pdf'),height=6,width=11)
oncoPrint(data[co_dna,],
          get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun_list, col = col, column_title = "Androgen Receptor",
          top_annotation = NULL, right_annotation = NULL,
          show_column_names = TRUE, remove_empty_rows = TRUE, remove_empty_columns = TRUE,
          heatmap_legend_param = hmap)
dev.off()


#pdf(paste0('plots-extended/oncoprint-androgen-receptor-',MIN_FREQ,'-',Sys.Date(),'.pdf'),height=4,width=11)
png(paste0('plots/oncoprint-androgen-receptor-',MIN_FREQ,'-',Sys.Date(),'.png'),height=250,width=600)
oncoPrint(data[path_ar,],#[,grepl('201|435',names(data))],
          get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun_list, col = col, column_title = "",
          top_annotation = NULL, right_annotation = NULL,
          show_column_names = TRUE, remove_empty_rows = TRUE, remove_empty_columns = TRUE,
          column_names_gp=gpar(fontsize=10),
          heatmap_legend_param = hmap)
dev.off()

#pdf(paste0('plots-extended/oncoprint-dna-repair-',MIN_FREQ,'-',Sys.Date(),'.pdf'),height=4,width=11)
png(paste0('plots/oncoprint-dna-repair-',MIN_FREQ,'-',Sys.Date(),'.png'),height=250,width=600)
oncoPrint(data[path_dna,],#[,grepl('201|435',names(data))],
          get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun_list, col = col, column_title = "",
          show_column_names = TRUE, remove_empty_rows = TRUE, remove_empty_columns = TRUE,
          top_annotation = NULL, right_annotation = NULL,
          column_names_gp=gpar(fontsize=9),
          heatmap_legend_param = hmap)
dev.off()


#pdf(paste0('plots-extended/oncoprint-wnt-',MIN_FREQ,'-',Sys.Date(),'.pdf'),height=3.5,width=10)
png(paste0('plots/oncoprint-wnt-',MIN_FREQ,'-',Sys.Date(),'.png'),height=250,width=600)
oncoPrint(data[path_wnt,],#[,grepl('201|435',names(data))],
          get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun_list, col = col, column_title = "",
          top_annotation = NULL, right_annotation = NULL,
          show_column_names = TRUE, remove_empty_rows = TRUE, remove_empty_columns = TRUE,
          column_names_gp=gpar(fontsize=9),
          heatmap_legend_param = hmap)
dev.off()

#pdf(paste0('plots-extended/oncoprint-pi3k-',MIN_FREQ,'-',Sys.Date(),'.pdf'),height=4,width=11)
png(paste0('plots/oncoprint-pi3k-',MIN_FREQ,'-',Sys.Date(),'.png'),height=250,width=600)
oncoPrint(data[path_pi3k,],#[,grepl('201|435',names(data))],
          get_type = function(x) strsplit(x, ";")[[1]],
          top_annotation = NULL, right_annotation = NULL,
          alter_fun = alter_fun_list, col = col, column_title = "",
          show_column_names = TRUE, remove_empty_rows = TRUE, remove_empty_columns = TRUE,
          column_names_gp=gpar(fontsize=8),
          heatmap_legend_param = hmap)
dev.off()

#pdf(paste0('plots-extended/oncoprint-ubi-',MIN_FREQ,'-',Sys.Date(),'.pdf'),height=5,width=10)
#oncoPrint(data[path_dna2,],#[,grepl('201|435',names(data))],
 #         get_type = function(x) strsplit(x, ";")[[1]],
#          alter_fun = alter_fun_list, col = col, column_title = "DNA2",
#          show_column_names = TRUE, remove_empty_rows = TRUE, remove_empty_columns = TRUE,
#          column_names_gp=gpar(fontsize=9),
#          heatmap_legend_param = hmap)
#dev.off()

