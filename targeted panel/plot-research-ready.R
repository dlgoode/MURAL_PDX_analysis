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
ready <- read_xlsx('~/Dropbox/Monash Staff PDX/All PDX Info/Targeted_Panel_v2/snapshot/pdx-samples-rename-for-paper-Nov-2020-snapshot.xlsx') %>% 
  filter(research_ready=='Y')

# Read variants
snv <- read_xlsx('all-PDX-variants-compact-2020-09-16-curate.xlsx') %>% filter(curate=='Y')
# Read cnvs
cnv.tw <- read.table('../twist-panel/cnv-out-2020-09-16.tsv',sep="\t",stringsAsFactors = F, header=T)
cnv.gr <- read.table('../garvan-panel/cnv-out-garvan-2020-09-03.tsv',sep="\t",stringsAsFactors = F, header=T)


out1 <- rbind(cnv.tw,cnv.gr)

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
#tmp <- rbind(tmp,c('CDG-031-294R-PDX','BRCA2','germline'))
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


ndx <- read_xlsx('~/Dropbox/Monash Staff PDX/All PDX Info/Targeted_Panel_v2/snapshot/pdx-samples-rename-for-paper-Nov-2020-snapshot.xlsx') %>% filter(research_ready=='Y')
map.ndx <- setNames(ndx$name3,ndx$name)
names(data) <- map.ndx[names(data)]
#order <- ndx[ndx$name3 %in% names(data),]$name3
order <- c("167.1R", "224R", "287R", "305R", "27.1A Cx", "27.2A Cx", "167.2M Cx", "201.1A Cx",
"201.2A Cx", "224R Cx", "305R Cx", "373M Cx", "374M Cx", "407M Cx", "410.51A Cx", "426M Cx", "435.1A Cx")

pdf(paste0('plots/oncoprint-ready-',MIN_FREQ,'-',Sys.Date(),'.pdf'),height=12,width=8)
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

png(paste0('plots/oncoprint-ready-',MIN_FREQ,'-',Sys.Date(),'.png'),height=900,width=600)
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
