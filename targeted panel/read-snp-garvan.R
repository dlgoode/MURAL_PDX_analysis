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
setwd('~/Dropbox/analysis/garvan-panel/')
#source('resources.R')
source('../twist-panel/resources.R')

snp <- read.multi.snv(loc='variants/',pat =".tsv")

back <- snp
names(snp) <- gsub('CDG.025.156R.PDX.','',names(snp))
snp$name <- gsub('_merged_sorted.markdups.realign.recal_HAP.bamstats_fieldsReduced_transcriptIndex.tsv','',snp$name)
names(snp) <- gsub(":","",names(snp))

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

phen <- read.table('phen_targeted_panel.csv',sep=",",stringsAsFactors = F,header = T) %>% dplyr::select(name=PDX,pdx=PDX.1,subline)
#dat <- merge(dat,phen,by="name")

dat.gar <- dat
cnvkit.index <- read.table('~/Dropbox/analysis/garvan-panel/index_rename_cnvkit_extended.csv',header=T,sep=",")
map.ndx <- setNames(cnvkit.index$NewID,cnvkit.index$CNVKIT)
dat$name <- map.ndx[dat$name]

dat$mutid <- paste(dat$CHROM,dat$POS,dat$REF,dat$ALT,dat$SYMBOL,sep="-")
dat <- merge(dat,phen,by="name")

paste0('shortlist-variants-garvan-',Sys.Date(),'.tsv')

write.table(dat,paste0('../garvan-panel/prefilter-variants-garvan-',Sys.Date(),'.tsv'),quote=T,row.names=F,col.names=T,sep="\t")
#dat <- read.table('../garvan-panel/prefilter-variants-garvan-2020-08-09.tsv',quote="\"",header=T,sep="\t",stringsAsFactors = F)
out <- filter(dat,HGVSc!=".") %>% dplyr::distinct(HGVSc)
write.table(out,paste0('pred-HGVSc-garvan-',Sys.Date(),'.tsv'),sep=",",quote=F,row.names=F,col.names=F)
