
#locations
working_dir_fname <- '~/Dropbox/analysis/twist-panel/'
oncokb_fname <- '/Users/bakshiandrew/Dropbox/Monash\ Staff\ PDX/All\ PDX\ Info/Targeted_Panel_v2/data/OncoKB_cancerGeneList_26_June_2020.tsv'
curated_genes_fname <- '/Users/bakshiandrew/Dropbox/Monash\ Staff\ PDX/All\ PDX\ Info/Targeted_Panel_v2/data/prad1000_CNA_Genes.txt'

setwd(working_dir)
source('resources.R')
library('reshape2')
library('ggrepel')
library('svglite')
library('RColorBrewer')
library('tidyverse')
setwd(working_dir_fname)

#Get set of important genes
oncokb <- read.table(oncokb_fname,sep="\t",header=T,quote="\"")
genes <- read.table(curated_genes_fname,header=T,sep="\t")
names(genes)[8] <- 'OncoKB'
genes$fq <- as.numeric(gsub('<','',gsub("%","",as.character(genes$Freq))))
genes <- filter(genes,Gene %in% oncokb$Hugo.Symbol | fq > 5)
gtmp <- read.table('~/Dropbox/analysis/twist-panel/curated-PCA-genes.csv',header=T)
gset <- unique(c(as.character(gtmp$gene),as.character(genes$Gene)))

cnr.dat <- read.multi(loc='cnvkit/',pat ="cnr.gm",sep = "\t")
cnr.dat$copy <- (2 ^ (cnr.dat$log2)) * 2
cnr.dat[cnr.dat$chromosome=='X',]$copy <- (2 ^ (cnr.dat[cnr.dat$chromosome=='X',]$log2)) 

cnv.dat <- cnr.dat %>% filter(gene %in% gset)
cnv.dat$cnv <- 'neutral'
cnv.dat[cnv.dat$chromosome!='X' & cnv.dat$copy >= 2.8,]$cnv <- 'gain'
cnv.dat[cnv.dat$chromosome!='X' & cnv.dat$copy >= 8,]$cnv <- 'gain-high'
cnv.dat[cnv.dat$chromosome!='X' & cnv.dat$copy <= 1.2,]$cnv <- 'loss'
cnv.dat[cnv.dat$chromosome!='X' & cnv.dat$copy <= 0.2,]$cnv <- 'homdel'

cnv.dat[cnv.dat$chromosome=='X' & cnv.dat$copy >= 2,]$cnv <- 'gain'
cnv.dat[cnv.dat$chromosome=='X' & cnv.dat$copy >= 8,]$cnv <- 'gain-high'
cnv.dat[cnv.dat$chromosome=='X' & cnv.dat$copy <= 0.2,]$cnv <- 'homdel'

#cnv.dat <- filter(cnv.dat,cnv!='neutral')

cnv.dat$name <- gsub('_merged_sorted.markdups.realign.cnr.gm.txt','',cnv.dat$name)
cnv.exn <- cnv.dat
cnv.exn$id <- paste0(cnv.exn$gene,"-",cnv.exn$cnv,"-",cnv.exn$name)
cnv.exn$id2 <- paste0(cnv.exn$gene,"-",cnv.exn$name)

cns.dat <- read.multi(loc='cnvkit/',pat ="cns.gm",sep = "\t")
cns.dat$copy <- (2 ^ (cns.dat$log2)) * 2
cns.dat[cns.dat$chromosome=='X',]$copy <- (2 ^ (cns.dat[cns.dat$chromosome=='X',]$log2)) 

cnv.dat <- cns.dat %>% filter(gene %in% gset)
cnv.dat$cnv <- 'neutral'
cnv.dat[cnv.dat$chromosome!='X' & cnv.dat$copy >= 2.8,]$cnv <- 'gain'
cnv.dat[cnv.dat$chromosome!='X' & cnv.dat$copy >= 8,]$cnv <- 'gain-high'
cnv.dat[cnv.dat$chromosome!='X' & cnv.dat$copy <= 1.2,]$cnv <- 'loss'
cnv.dat[cnv.dat$chromosome!='X' & cnv.dat$copy <= 0.2,]$cnv <- 'homdel'

cnv.dat[cnv.dat$chromosome=='X' & cnv.dat$copy >= 2,]$cnv <- 'gain'
cnv.dat[cnv.dat$chromosome=='X' & cnv.dat$copy >= 8,]$cnv <- 'gain-high'
cnv.dat[cnv.dat$chromosome=='X' & cnv.dat$copy <= 0.2,]$cnv <- 'homdel'

#cnv.dat <- filter(cnv.dat,cnv!='neutral')

cnv.dat$name <- gsub('_merged_sorted.markdups.realign.cns.gm.txt','',cnv.dat$name)
cnv.dat$id <- paste0(cnv.dat$gene,"-",cnv.dat$cnv,"-",cnv.dat$name)
cnv.dat$id2 <- paste0(cnv.dat$gene,"-",cnv.dat$name)
#View(filter(cnv.dat,cnv=='gain-high'))
tmp.exn <- cnv.exn[!cnv.exn$id %in% cnv.dat$id,] %>% filter(cnv=='homdel')
tmp.exn$info <- 'exon'
cnv.dat$info <- ''
cnv.out <- rbind(dplyr::select(cnv.dat,-segment_weight,-segment_probes),
                 tmp.exn)

write.table(cnv.out,paste0('cnv-out-',Sys.Date(),'.tsv'),quote=F,row.names=F,col.names=T,sep="\t")


