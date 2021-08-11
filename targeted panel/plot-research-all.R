source('resources.R')
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
#source('make_phenotype_info.R')

# Read metadata
# Data labelled with 'X' in research ready should be eXcluded
#ready <- read_xlsx('pdx-samples-rename-for-paper-Nov-2020.xlsx') %>% filter(research_ready!='X')
ready <- read_xlsx('~/Dropbox/Monash Staff PDX/All PDX Info/Targeted_Panel_v2/snapshot/pdx-samples-rename-for-thesis-Nov-2020-snapshot.xlsx') %>% filter(research_ready!='X')

# Read variants
snv <- read_xlsx('all-PDX-variants-compact-2020-09-16-curate.xlsx') %>% filter(curate=='Y')

# Read cnvs
cnv.tw <- read.table('../twist-panel/cnv-out-2020-09-16.tsv',sep="\t",stringsAsFactors = F, header=T)
cnv.gr <- read.table('../garvan-panel/cnv-out-garvan-2020-09-03.tsv',sep="\t",stringsAsFactors = F, header=T)

#####################################################
## write snvs for publication
#out2.all <- read.table('../chapter1/all-PDX-variants-2020-09-16.tsv',header=T,sep="\t",quote="\"",stringsAsFactors = F)
#out2.remove <- read_xlsx('../chapter1/all-PDX-variants-compact-2020-09-16-curate.xlsx') %>% filter(curate=='N')
#rename <- read_xlsx('~/Dropbox/Monash Staff PDX/All PDX Info/Targeted_Panel_v2/snapshot/pdx-samples-rename-for-thesis-Nov-2020-snapshot.xlsx') %>% dplyr::select(name,name3)
#keepthis <- out2.all %>% filter(name %in% ready$name) %>% filter(!mutid %in% out2.remove$mutid)
##out2.all %>% filter(name %in% ready$name) %>% distinct(name)
#keepthis[keepthis$mutid %in% out2.all$mutid,]$curate <- 'Y'
#keepthis <- merge(keepthis,rename,by="name")
#write.table(keepthis,'all-PDX-variants-2020-09-16-curate-remove-unwanted-samples.tsv',
#            quote=F,row.names=F,sep="\t",col.names=T)
#keep2 <- keepthis %>% dplyr::select(name3,SYMBOL,HGVSc,dbSNP_ids,COSMIC_ids,POS,REF,ALT,DEPTH=PMCDP,ALT_FREQ=PMCFREQ,
#                                    Consequence,EXON,HGVSp,PolyPhen,SIFT,CADD_PHRED,gnomAD_AF,clinvar_clnsig,
#                                    clinvar_id,fathmm.MKL_coding_pred,MetaSVM_pred,MetaLR_pred,curate)
#write.table(keep2,'all-PDX-variants-2020-09-16-curate-remove-unwanted-samples-short.tsv',
#            quote=F,row.names=F,sep="\t",col.names=T)
######################################################

out1 <- rbind(cnv.tw,cnv.gr)

#######################################################
# write copy number variants for paper - all variants
out1.keep <- merge(out1,rename,by="name") %>% dplyr::select(name3,gene,chromosome,start,end,log2,copy,cnv)
write.table(out1.keep,'cnv-out-2020-09-16-rename-all.tsv',quote=F,row.names=F,col.names=T,sep="\t")
#######################################################

cur.genes <- read.table('~/Dropbox/Monash Staff PDX/All PDX Info/Targeted_Panel_v2/data/curated-PCA-genes.csv',header=T)
cna.genes <- read.table('~/Dropbox/Monash Staff PDX/All PDX Info/Targeted_Panel_v2/data/prad1000_CNA_Genes.txt',header=T,sep="\t")
names(cna.genes)[8] <- 'OncoKB'

### FILTERS
MIN_FREQ <- 0.75  #minimum frequency to include in plot
#THR_GAIN <- 2.8   #theshold for curated copy gain
#THR_LOSS <- 1.2   #1 #threshold for curate copy loss
DB_FREQ <- 1      #minimum frequency in Armenia to count as recurrent copy number

out1 <- filter(out1,name %in% ready$name)
out2 <- filter(snv,name %in% ready$name)

out2 <- filter(out2,PMCFREQ2>MIN_FREQ | (curate=='Y' & mutid=='17-7574029-CG-C-TP53'))
out2 <- out2 %>% dplyr::rename(gene=SYMBOL)
snv.genes <- unique(out2$gene)

out1 <- merge(out1,ready,by="name")
out1 <- filter(out1,gene %in% cur.genes$gene | gene %in% snv.genes | cnv=='gain-high')#,copy>=THR_GAIN | copy < THR_LOSS)
#out1 <- filter(out1,copy>=THR_GAIN | copy < THR_LOSS)


cna.genes <- filter(cna.genes,OncoKB=='Yes')
cna.genes$fq <- as.numeric(gsub('<','',gsub("%","",as.character(cna.genes$Freq))))
cna.genes.gain <- cna.genes %>% filter(fq>DB_FREQ,CNA=='AMP') %>% dplyr::distinct(Gene,CNA,fq)
cna.genes.del <- cna.genes %>% filter(fq>DB_FREQ,CNA=='HOMDEL') %>% dplyr::distinct(Gene,CNA,fq)
out1.amp <- out1[grepl('gain',out1$cnv),] %>% filter(!gene %in% cna.genes.del$Gene | gene %in% cna.genes.gain$Gene)
out1.loss <- out1[grepl('loss',out1$cnv),] %>% filter(!gene %in% cna.genes.gain$Gene | gene %in% cna.genes.del$Gene)
out1.homdel <- out1[grepl('homdel',out1$cnv),]
out1 <- rbind(rbind(out1.amp,out1.loss),out1.homdel)

###### write copy number variants for paper - curated
out1.keep <- out1 %>% dplyr::select(name3,gene,chromosome,start,end,log2,copy,cnv)
write.table(out1.keep,'cnv-out-2020-09-16-rename.tsv',quote=F,row.names=F,col.names=T,sep="\t")
#######################################################

out1 <- out1 %>% dplyr::select(id=name,gene,type=cnv) 
out2 <- out2 %>% dplyr::select(id=name,gene,type=Consequence) 

tmp <- rbind(out2,out1)
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
tmp[grepl('missense',tmp$type),]$type <- 'MUT'
tmp[grepl('splice_acceptor_variant',tmp$type),]$type <- 'MUT'
tmp[grepl('splice_donor_variant',tmp$type),]$type <- 'MUT'
tmp[grepl("inframe_insertion,splice_region_variant",tmp$type),]$type <- 'MUT'

#manually add germline, somatic variants 
tmp <- rbind(tmp,c('CDG-031-294R-PDX','BRCA2','germline'))
tmp <- rbind(tmp,c('44-CA27LN2-Cx5','AR','LOSS'))

#count genes to get subset by frequency
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

hmap <-  list(title = "Alterations", at = c("AMP","AMPHIGH","HOMDEL","LOSS","STP","MUT","germline"), 
              labels = c("Amplification","Amplification 8+","Deep deletion","Partial Loss",
                         "Stop Gain / Frameshift","Mutation","Germline"))


ndx <- read_xlsx('~/Dropbox/Monash Staff PDX/All PDX Info/Targeted_Panel_v2/snapshot/pdx-samples-rename-for-thesis-Nov-2020-snapshot.xlsx') %>% filter(research_ready!='X')
map.ndx <- setNames(ndx$name3,ndx$name)
names(data) <- map.ndx[names(data)]
order <- ndx[ndx$name3 %in% names(data),]$name3

#hacky way to add NE
#backup <- data['SUFU',]
#data['SUFU',]='MUT;'
#
#pdf(paste0('plots/oncoprint-subset-',MIN_FREQ,'-',Sys.Date(),'.pdf'),height=13,width=14)
#oncoPrint(data[c('SUFU','TP53','RB1','PTEN','MYC','AR'),],
#          get_type = function(x) strsplit(x, ";")[[1]],
#          alter_fun = alter_fun_list, col = col, 
#          column_title = "",
#          show_column_names = TRUE,
#          remove_empty_rows = FALSE,
#          column_order = order,
#          remove_empty_columns = FALSE,
#          heatmap_legend_param = hmap)
#dev.off()

pdf(paste0('plots/oncoprint-all-',MIN_FREQ,'-',Sys.Date(),'.pdf'),height=13,width=14)
oncoPrint(data,
          get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun_list, col = col, 
          column_title = "PDX",
          show_column_names = TRUE,
          remove_empty_rows = TRUE,
          column_order = order,
          remove_empty_columns = FALSE,
          heatmap_legend_param = hmap)
dev.off()

png(paste0('plots/oncoprint-all-',MIN_FREQ,'-',Sys.Date(),'.png'),height=1000,width=1200)
oncoPrint(data,
          get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun_list, col = col, 
          column_title = "PDX",
          show_column_names = TRUE,
          remove_empty_rows = TRUE,
          column_order = order,
          remove_empty_columns = FALSE,
          heatmap_legend_param = hmap)
dev.off()

# Re-order genes
#ndx <- read_xlsx('pdx-samples-rename-for-paper-Nov-2020.xlsx') %>% filter(research_ready!='X')
#map.ndx <- setNames(ndx$name2,ndx$name)
#names(data) <- map.ndx[names(data)]
#order <- ndx[ndx$name2 %in% names(data),]$name2
#
#pdf(paste0('plots/oncoprint-all-rename-',MIN_FREQ,'-',Sys.Date(),'.pdf'),height=13,width=14)
#oncoPrint(data,
#          get_type = function(x) strsplit(x, ";")[[1]],
#          alter_fun = alter_fun_list, col = col, 
#          column_title = "PDX",
#          show_column_names = TRUE,
#          remove_empty_rows = TRUE,
#          column_order = order,
#          remove_empty_columns = FALSE,
#          heatmap_legend_param = hmap)
#dev.off()
#
#png(paste0('plots/oncoprint-all-rename-',MIN_FREQ,'-',Sys.Date(),'.png'),height=1000,width=1200)
#oncoPrint(data,
#          get_type = function(x) strsplit(x, ";")[[1]],
#          alter_fun = alter_fun_list, col = col, 
#          column_title = "PDX",
#          show_column_names = TRUE,
#          remove_empty_rows = TRUE,
#          column_order = order,
#          remove_empty_columns = FALSE,
#          heatmap_legend_param = hmap)
#dev.off()
#
###Plot all by pathway
cur.genes <- read.table('~/Dropbox/analysis/twist-panel/curated-PCA-genes.csv',header=T)
twist <- read.table('/Users/bakshiandrew/Dropbox/Monash\ Staff\ PDX/All\ PDX\ Info/Targeted_Panel_v2/gene-lists/twist-genelist-2020.txt',stringsAsFactors = F,header=T)

path_ar <- c('AR','NCOR1','NCOR2','FOXA1','SPOP','SPEN') #AR pathway - remove SPEN?
path_dna <- c('BRCA2','ATM','BRCA1','FANCA','RAD51C','RAD51','MLH1','MSH2') #DNA repair
path_dna2 <- c('ATM', 'BRCA1', 'BRIP1', 'CHEK2', 'EPCAM', 'HOXB13', 'MLH1', 'MSH2', 'MSH6', 'NBN', 'PALB2', 'PMS2', 'RAD51C', 'RAD51D', 'TP53','FANCA','CDK12')
path_wnt <- c('APC','CTNNB1','RNF43','AXIN1','AXIN2','ZNRF3','RSPO2') #WNT pathway
path_pi3k <- c('PTEN','PIK3CA','PIK3CB','PIK3R1','AKT1')
path_ubi <- c('SF3B1','U2AF1','GEMIN5','TCERG1','PRPF8') #Ubiquitin and Splicing pathways


#lookup names in order

path_dna2[path_dna2 %in% twist$gene]
path_dna2[!path_dna2 %in% twist$gene]

hmap <-  list(title = "Alterations", at = c("AMP","AMPHIGH","HOMDEL","LOSS","STP","MUT","germline"), 
              labels = c("Amplification","Amplification 8+","Deep deletion","Partial Loss","Stop Gain / Frameshift","Mutation","Germline"))




pdf(paste0('plots/oncoprint-androgen-receptor-',MIN_FREQ,'-',Sys.Date(),'.pdf'),height=6,width=11)
oncoPrint(data[path_ar,],
          get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun_list, col = col, column_title = "Androgen Receptor",
          top_annotation = NULL, right_annotation = NULL,
          show_column_names = TRUE, remove_empty_rows = TRUE, remove_empty_columns = TRUE,
          heatmap_legend_param = hmap)
dev.off()

pdf(paste0('plots/oncoprint-dna-repair-',MIN_FREQ,'-',Sys.Date(),'.pdf'),height=6,width=11)
oncoPrint(data[path_dna,],
          get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun_list, col = col, column_title = "DNA repair",
          show_column_names = TRUE, remove_empty_rows = TRUE, remove_empty_columns = TRUE,
          top_annotation = NULL, right_annotation = NULL,
          heatmap_legend_param = hmap)
dev.off()


pdf(paste0('plots/oncoprint-wnt-',MIN_FREQ,'-',Sys.Date(),'.pdf'),height=5,width=10)
oncoPrint(data[path_wnt,],
          get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun_list, col = col, column_title = "WNT",
          top_annotation = NULL, right_annotation = NULL,
          show_column_names = TRUE, remove_empty_rows = TRUE, remove_empty_columns = TRUE,
          heatmap_legend_param = hmap)
dev.off()

pdf(paste0('plots/oncoprint-pi3k-',MIN_FREQ,'-',Sys.Date(),'.pdf'),height=5,width=11)
oncoPrint(data[path_pi3k,],
          get_type = function(x) strsplit(x, ";")[[1]],
          top_annotation = NULL, right_annotation = NULL,
          alter_fun = alter_fun_list, col = col, column_title = "PI3K",
          show_column_names = TRUE, remove_empty_rows = TRUE, remove_empty_columns = TRUE,
          heatmap_legend_param = hmap)
dev.off()

