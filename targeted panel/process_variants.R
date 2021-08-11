# *****************************************************************************
# Process Variants
# -----------------
# Andrew Bakshi 2021
# Read in pre-filtered variants, and filter according to parameters described
# in PDX paper manuscript, remove very likely germline, write out likely germline
# for review (then remove these too unless something stands out that should be included)
# then write to files

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
setwd('~/Dropbox/analysis/chapter1/')
source('lib_variants.R')

# *****************************************************************************
# INPUT

tmp1 <- read.table('../garvan-panel/prefilter-variants-garvan-2020-08-09.tsv',quote="\"",header=T,sep="\t",stringsAsFactors = F)
tmp1$batch <- 'g'
tmp2 <- read.table('../twist-panel/prefilter-variants-twist-2020-09-15.tsv',quote="\"",header=T,sep="\t",stringsAsFactors = F)
tmp <- rbind(tmp1,tmp2)

source('../chapter1/readvep.R')
tmp.vep <- merge(tmp,vep,by="HGVSc",all.x=T)

# *****************************************************************************
# FILTERING

tmp.keep <- tmp.vep %>% filter(PMCFREQ2>0.25,gnomAD_AF2<0.001 | gnomAD_AF==".") %>%
                        filter(IMPACT=='HIGH' |
                                IMPACT=='MODERATE' &
                                (
                                 MetaSVM_pred=='D' | 
                                 MetaLR_pred=='D' |
                                 grepl('pathogenic',clinvar_clnsig) |
                                 grepl('Pathogenic',clinvar_clnsig) |
                                 clinvar_clnsig=='Pathogenic' |
                                    (
                                     grepl('^deleterious',SIFT) & 
                                     (grepl('^damaging',PolyPhen) | grepl('probably_damaging',PolyPhen)) & 
                                     `fathmm-MKL_coding_pred`=='D'
                                     ) |
                                 CADD_PHRED2 >= 30 | CADD_PHRED=="."
                                )
                              )

#manually exclude those that had sample swap or QC issuies
exclude <- c('8-373M-cx6','2-27-2-cx9','14-435A-31-cx5','41-435A-31-Cx2','twist-2-43-373M-Cx11','45-CA27LN2-Cx14')
tmp.keep <- filter(tmp.keep,!name %in% exclude)

#*****************************************************************************
# REVIEW VERY LIKELY & POSSIBLE GERMLINE
# Manually review these - as of latest data, all of these seem ok to remove,
# hence we just remove anything in tmp.dup

# ALMOST CETAIN GERMLINE 
# snps that occur in 10+ samples are so common they are likely germline snps
tmp.germline <- tmp.keep %>% filter(PMCFREQ2>0.25) %>% group_by(mutid) %>% dplyr::distinct(mutid,pdx,IMPACT) %>% dplyr::count(mutid,IMPACT) %>% filter(n>=10)

# LIKELY GERMLINE 
#snps that occur in less than 10, but greater than 3 individual PDXs (not sublines) are probably duplicates,
#but they need to be manually reviewed
tmp.dup <- tmp.keep %>% filter(PMCFREQ2>0.25) %>% group_by(mutid) %>% dplyr::distinct(mutid,pdx,IMPACT) %>% dplyr::count(mutid,IMPACT) %>% filter(n>=3 & n<10)
tmp.out <- tmp.keep %>% filter(mutid %in% tmp.dup$mutid) %>% arrange(mutid)

# Manually review this file, if you want to change it,
# read it back in and overwrite tmp.dup with it
write.table(dplyr::select(tmp.out,mutid,short),
            paste0('manual-review-if-artifact',Sys.Date(),'.tsv'),
            quote=T,row.names=F,col.names=T,sep="\t")

# list of column names to keep
short <- c('SYMBOL','name','Consequence','PMCFREQ2','PMCAD2','PMCDP','CHROM','POS','COSMIC_ids','HGVSp','HGVSc',
           'PolyPhen','SIFT','CADD_PHRED','gnomAD_AF','ID','CHROM','POS','REF','ALT','clinvar_clnsig','MetaLR_pred','MetaSVM_pred','fathmm-MKL_coding_pred')

tmp.keep <- tmp.keep %>% filter(!mutid %in% tmp.germline$mutid,!mutid %in% tmp.dup$mutid) 

#these columns manually filled during manual review (ie. in excel)
tmp.keep$curate <- ""
tmp.keep$note <- ""

# *****************************************************************************
# OUTPUT
# write out data to various files, some with more columns

write.table(dplyr::select(tmp.keep,mutid,short,curate,note) %>% arrange(mutid),
            paste0('keep-all-variants',Sys.Date(),'.tsv'),
            quote=T,row.names=F,col.names=T,sep="\t")

write.table(tmp.keep,
            paste0('all-PDX-variants-',Sys.Date(),'.tsv'),
            quote=T,row.names=F,col.names=T,sep="\t")

write.table(dplyr::select(tmp.keep,short,subline,IMPACT,EXON,INTRON,mutid),
            paste0('all-PDX-variants-compact-',Sys.Date(),'.tsv'),
            quote=T,row.names=F,col.names=T,sep="\t")

ready <- read_xlsx('pdx-samples-rename-for-paper-Nov-2020.xlsx') %>% filter(research_ready=='Y')

write.table(filter(tmp.keep,name %in% ready$name) %>% dplyr::select(short,subline,IMPACT,EXON,INTRON,mutid,'clinvar_clnsig','MetaLR_pred','MetaSVM_pred','fathmm-MKL_coding_pred','REVEL_score','curate','note'),
            paste0('all-PDX-variants-research-ready-compact-',Sys.Date(),'.tsv'),
            quote=T,row.names=F,col.names=T,sep="\t")

