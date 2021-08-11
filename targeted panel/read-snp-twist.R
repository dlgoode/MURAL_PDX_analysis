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
setwd('~/Dropbox/analysis/twist-panel/')
source('resources.R')

#read this as template for which columns to select
group1 <- read.multi.snv(loc='variants/',pat="fieldsReduced_transcriptIndex.tsv")
names(group1) <- gsub("X1.27.1.cx4.","",names(group1))
group1$name <- gsub('_merged_sorted.markdups.realign.recal_HAP.bamstats_fieldsReduced_transcriptIndex.tsv','',group1$name)

snp <- read.multi.snv(loc='variants/',pat ="bamstats_transcriptIndex.tsv")
back <- snp
snp$name <- gsub('_merged_sorted.markdups.realign.recal_HAP.bamstats_transcriptIndex.tsv','',snp$name)
names(snp) <- gsub("X41.435A.31.Cx2.","",names(snp))
snp <- snp[,names(snp) %in% names(group1)] %>% dplyr::select(-DP.1)

snp <- rbind(snp,group1)

snp$PMCFREQ2 <- snp$PMCFREQ
snp[snp$PMCFREQ2==".",]$PMCFREQ2 <- '0'
snp$PMCFREQ2 <- as.numeric(snp$PMCFREQ2)

snp$PMCAD2 <- snp$PMCAD
snp[snp$PMCAD2==".",]$PMCAD2 <- '0'
snp$PMCAD2 <- as.numeric(snp$PMCAD2)

snp$CADD_PHRED2 <- snp$CADD_PHRED
snp[snp$CADD_PHRED2==".",]$CADD_PHRED2 <- '0'
snp$CADD_PHRED2 <- as.numeric(snp$CADD_PHRED2)

snp$gnomAD_AF2 <- snp$gnomAD_AF
snp[snp$gnomAD_AF2==".",]$gnomAD_AF2 <- '0'
snp$gnomAD_AF2 <- as.numeric(snp$gnomAD_AF2)

dat <- snp %>% filter(PMCFREQ2>0.25,
                      PMCAD2>=10,
                      CANONICAL=='YES',
                      IMPACT %in% c('HIGH','MODERATE'))

dat$mutid <- paste(dat$CHROM,dat$POS,dat$REF,dat$ALT,dat$SYMBOL,sep="-")
info <- read.table('twist-panel-details.csv',sep=",",stringsAsFactors = F,header=T)
back.dat <- dat
dat <- merge(dat,info,by="name")
write.table(dat,paste0('../twist-panel/prefilter-variants-twist-',Sys.Date(),'.tsv'),quote=T,row.names=F,col.names=T,sep="\t")
#dat <- read.table('../twist-panel/prefilter-variants-twist-2020-08-08.tsv',quote="\"",header=T,sep="\t",stringsAsFactors = F)
out <- filter(dat,HGVSc!=".") %>% dplyr::distinct(HGVSc)
write.table(out,paste0('pred-HGVSc-twist-',Sys.Date(),'.tsv'),sep=",",quote=F,row.names=F,col.names=F)
