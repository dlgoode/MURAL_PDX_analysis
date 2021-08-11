source('~/Dropbox/Monash Staff PDX/All PDX Info/Targeted_Panel_v2/src/lib/resources.R')
source('~/Dropbox/Monash Staff PDX/All PDX Info/Targeted_Panel_v2/src/lib/useful.R')
library('ggtree')
library(dplyr)
library(ggplot2)
library(tidyverse)
library(singscore)
library(msigdbr)
library(edgeR)
library(gplots)
library(ape)
library(RColorBrewer)
library(readxl)
library(heatmap3)
setwd('~/Dropbox/analysis/pdx-quantseq/')

#use this one - AB 2020

data <- read.table('~/Dropbox/Monash Staff PDX/All PDX Info/QuantSeq/QuantSeq-MURAL-PDX-counts.txt',header=T,sep="\t",stringsAsFactors = F)
sample_info <- read.table('~/Dropbox/Monash Staff PDX/All PDX Info/Targeted_Panel_v2/data/sample-info-july2019.txt',header=T,sep="\t",stringsAsFactors = F)
sample_info <- filter(sample_info,!grepl('Original',sample_type))
counts.keep <- filterRNASEQ(data,threshold=10,cpm_min=1)
dge <- DGEList(counts=counts.keep)

log.cpm <- cpm(dge,log=T)
log.cpm <- log.cpm[,colnames(log.cpm) %in% sample_info$sample]
dat.cpm <- add_symbol(log.cpm,gene_names)
row.names(dat.cpm) <- dat.cpm$Symbol
dat.cpm <- dplyr::select(dat.cpm,-Symbol,-GeneID)

dat.cpm <- dat.cpm[,phen$sample]

phen <- read_xlsx('~/Dropbox/Monash Staff PDX/All PDX Info/Targeted_Panel_v2/snapshot/pdx-samples-for-paper-quantseq-Nov-2020.xlsx')
names(phen)[1] <- 'sample'
names(phen)[7] <- 'paperlist'
phen$source <- 'metastatic'
phen[phen$Sample.Site=='Prostate',]$source <- 'primary'
phen$type <- phen$pathology

pca <- prcomp(t(dat.cpm))
pca.tmp <- as.data.frame(pca$x)
pca.tmp$sample <- row.names(pca$x)
pca.tmp <- merge(pca.tmp,phen,by="sample")
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
group.colors <- c('mediumpurple4','goldenrod2','tomato3')
pca.tmp$shape <- paste0(pca.tmp$source,"-",pca.tmp$paperlist)

print("plotting")
name<-paste0('PDX-',Sys.Date())
ggplot(pca.tmp) + 
  geom_point(data=pca.tmp,aes(x=PC1,y=PC2,color=as.factor(pathology),shape=as.factor(source),size=4),alpha=0.95) +
  scale_color_manual(values=group.colors) +
  xlab(paste0("PC1: var explained: ",format(pca.var.per[1],digits=3),"%")) + 
  ylab(paste0("PC2: var explained: ",format(pca.var.per[2],digits=3),"%")) +
  theme_classic(base_size = 16) + theme(legend.position="right") +
  theme(legend.title=element_blank()) +
  scale_alpha(guide = 'none') + scale_size(guide = 'none')
ggsave(paste0('plots/PCA-',name,'-label.pdf'),width=8,height=6.5)
ggsave(paste0('plots/PCA-',name,'-label.png'),width=8,height=6.5)

ggplot(pca.tmp) + 
  geom_point(data=pca.tmp,aes(x=PC1,y=PC2,color=as.factor(pathology),shape=as.factor(source),size=4),alpha=0.95) +
  scale_color_manual(values=group.colors) +
  geom_label_repel(aes(x=PC1,y=PC2,label=sample),alpha=0.85, 
  box.padding = 0.25, size=4, show.legend=F,
  point.padding = 0.25, segment.color='black') +
  xlab(paste0("PC1: var explained: ",format(pca.var.per[1],digits=3),"%")) + 
  ylab(paste0("PC2: var explained: ",format(pca.var.per[2],digits=3),"%")) +
  theme_classic(base_size = 16) + theme(legend.position="right") +
  theme(legend.title=element_blank()) +
  scale_alpha(guide = 'none') + scale_size(guide = 'none')
ggsave(paste0('plots/PCA-',name,'-detailed.pdf'),width=8*2,height=6.5*2)
ggsave(paste0('plots/PCA-',name,'-detailed.png'),width=8*2,height=6.5*2)


print("scree plot")
png(paste0('plots/skree-plot',name,'.png'),height=300,width=400)
names(pca.var.per) <- 1:20
barplot(pca.var.per[1:20]) + title('Variance explained by first 20 PCs')
dev.off()

print("plotting")
name<-'PDX-all-PC3PC4'
ggplot(pca.tmp) + 
  geom_point(data=pca.tmp,aes(x=PC3,y=PC4,color=as.factor(pathology),shape=as.factor(source),size=4),alpha=0.95) +
  scale_color_manual(values=group.colors) +
  xlab(paste0("PC3: var explained: ",format(pca.var.per[3],digits=3),"%")) + 
  ylab(paste0("PC4: var explained: ",format(pca.var.per[4],digits=3),"%")) +
  theme_classic(base_size = 16) + theme(legend.position="right") +
  theme(legend.title=element_blank()) +
  scale_alpha(guide = 'none') + scale_size(guide = 'none')
ggsave(paste0('plots/PCA-',name,'-label.pdf'),width=8,height=6.5)
ggsave(paste0('plots/PCA-',name,'-label.png'),width=8,height=6.5)


options(digits=3)

for (i in 1:4) {
  loading_scores <- pca$rotation[,i]
  gene_scores <- loading_scores ## get the magnitudes
  gene_score_ranked <- sort(abs(gene_scores), decreasing=TRUE)
  top_10_genes <- names(gene_score_ranked[1:10])
  filter(gene_names,Symbol %in% top_10_genes)
  gscores <- as.data.frame(gene_scores)
  gscores$Symbol <- row.names(gscores)
  int.load <- filter(gene_names,Symbol %in% names(boxplot(gene_scores,plot=FALSE)$out))
  int.load <- merge(int.load,gscores,by="Symbol")
  arrange(int.load,desc(abs(int.load$gene_scores))) %>% head(20)
  int.load$gene_scores <- as.numeric(formatC(int.load$gene_scores, digits = 4, format = "f"))
  write.table(arrange(int.load,desc(abs(gene_scores))) %>% head(20),paste0('data/PC',i,'-loadings.txt'),quote=F,row.names=F,col.names=T,sep="\t")
}

options(digits=8)


library(gplots)

mypalette <- brewer.pal(9,"YlGnBu")
morecols <- colorRampPalette(mypalette)
var_genes <- apply(dat.cpm[,phen$sample],1,var)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:100]
dat.cpm.var <- dat.cpm[select_var,phen$sample]
dim(dat.cpm.var)
group.colors <- c('mediumpurple4','goldenrod2','tomato3')
col.cell <- col.cell <- group.colors[as.factor(phen$type)]
png(filename = paste('plots/heatmap-ALL-50.png',sep="/"),
    width = 1040, height = 824, units = "px", pointsize = 12, bg = "white")
heatmap.2(as.matrix(dat.cpm.var[1:50,]),col=rev(morecols(20)),trace="none", main="Top 50 most variable genes across all samples",
          ColSideColors=col.cell,scale="row", keysize = 1, cexRow = 1, cexCol=1,margins=c(14,10))
dev.off()

phen.back <- phen
phen <- filter(phen,paperlist=='Y')

library(GSEABase)
m_df = msigdbr(species = "Homo sapiens", category = "H")
gsets <- unique(m_df$gs_name)
gmt <- read.msig(gsets)
gmt <- GeneSetCollection(c(gmt))
singscore.out <- run.singscore(dat.cpm[,phen$sample],gmt)
singdata <- t(singscore.out %>% dplyr::select(-sample))
singtidy <- gather(singscore.out,key='set',value='score',-sample) #%>% filter(set %in% )

rownames(singdata) <- gsub("HALLMARK_","",rownames(singdata))

mypalette <- brewer.pal(9,"RdYlBu")
morecols <- colorRampPalette(mypalette)
col.cell <- group.colors[as.factor(phen[phen$sample!='M167M-R2',]$type)]


png(filename = paste('plots/HALLMARK-MSIGDB50-scaled.png',sep="/"),
    width = 1040, height = 824, units = "px", pointsize = 12, bg = "white")
heatmap.2(as.matrix(singdata[,colnames(singdata) != "M167M-R2"]),col=rev(morecols(25)),trace="none", main="MSigDB Hallmark 50 Gene Sets",
          ColSideColors=col.cell,scale="row", keysize = 1, cexRow = 1, cexCol=1,margins=c(8,20),
          Rowv=F,Colv=T,dendrogram="column")
dev.off()

png(filename = paste('plots/HALLMARK-MSIGDB50.png',sep="/"),
    width = 1040, height = 824, units = "px", pointsize = 12, bg = "white")
heatmap.2(as.matrix(singdata[,colnames(singdata) != "M167M-R2"]),col=rev(morecols(25)),trace="none", main="MSigDB Hallmark 50 Gene Sets",
          ColSideColors=col.cell,scale="none", keysize = 1, cexRow = 1, cexCol=1,margins=c(8,20),
          Rowv=T,Colv=T,dendrogram="column")
dev.off()

kegg <- c('KEGG_BASE_EXCISION_REPAIR','KEGG_HOMOLOGOUS_RECOMBINATION',
          'KEGG_MISMATCH_REPAIR','KEGG_NON_HOMOLOGOUS_END_JOINING','KEGG_NUCLEOTIDE_EXCISION_REPAIR')
gmt <- read.msig(kegg)
gmt <- GeneSetCollection(c(gmt))
singscore.out <- run.singscore(dat.cpm[,phen$sample],gmt)
singdata <- t(singscore.out %>% dplyr::select(-sample))
singtidy <- gather(singscore.out,key='set',value='score',-sample) #%>% filter(set %in% )
singtidy <- merge(singtidy,phen,by="sample")

rownames(singdata) <- gsub("KEGG_","",rownames(singdata))

source.colors <- c("mediumseagreen","orange2")
ggplot(singtidy) + geom_boxplot(aes(set,score,fill=as.factor(source)),width=0.75) + 
  coord_flip() + scale_fill_manual(values=source.colors) + theme_minimal() +
  theme(legend.title=element_blank()) + xlab("") + ylab("") #+ ggtitle('')
ggsave('plots/KEGG-DNA-REPAIR-boxplot.png',width=10,height=6)

#ggplot(singtidy[singtidy$name3 %in% phen[phen$host=='Castrate',]$name3,]) + 
#  geom_boxplot(aes(set,score,fill=as.factor(host)),width=0.75) + 
#  coord_flip() + scale_fill_manual(values=c("grey","black")) + theme_minimal() +
#  theme(legend.title=element_blank()) + xlab("") + ylab("") #+ ggtitle('')
#ggsave('plots/KEGG-DNA-REPAIR-boxplot-cx-in.png',width=10,height=6)


singtidy$MMR <- 0
singtidy[singtidy$name3 %in% c('201.1A','201.1A Cx'),]$MMR <- '201.1A'
singtidy[singtidy$name3 %in% c('272R'),]$MMR <- '272R'
singtidy[singtidy$name3 %in% c('287R'),]$MMR <- '287R'
singtidy$HRD <- 0
singtidy[singtidy$name3 %in% c('294R'),]$HRD <- '294R'
singtidy[singtidy$name3 %in% c('435.10A','435.1A','435.23A','435.7A'),]$HRD <- '435M'
names(singtidy)
phen$MMR <- 0
phen[phen$name3 %in% c('201.1A','201.1A Cx','272R','287R'),]$MMR <- 1
phen$HRD <- 0
phen[phen$name3 %in% c('294R','435.10A','435.1A','435.23A','435.7A'),]$HRD <- 1

t.test(singtidy[singtidy$MMR=='272R',]$score,singtidy[singtidy$MMR==0,]$score,alternative="less")
t.test(singtidy[singtidy$HRD=='294R',]$score,singtidy[singtidy$HRD==0,]$score,alternative="less")

library(ggrepel)
singtidy$label <- paste0(singtidy$name3,"_",singtidy$host)
singtidy$state <- 'IN'
singtidy[singtidy$state.x]
singtidy$set <- gsub("KEGG_","",singtidy$set) 
ggplot(singtidy %>% filter(set=='MISMATCH_REPAIR')) + geom_boxplot(aes(scale(score),set,color=as.factor(MMR),fill=as.factor(MMR)),width=0.4,alpha=0.25) +
  theme_classic() + xlab('normalised score')
#  geom_label_repel(data=filter(singtidy,MMR==1,set=='MISMATCH_REPAIR'),aes(scale(score),set,label=label),alpha=.8, 
#                   box.padding = 0.25, size=3, show.legend=F,
#                   point.padding = 0.25, segment.color='black') + ylab("")

ggplot(singtidy %>% filter(set=='HOMOLOGOUS_RECOMBINATION')) + geom_boxplot(aes(scale(score),set,color=as.factor(HRD),fill=as.factor(HRD)),width=0.4,alpha=0.25) +
  theme_classic() + xlab('normalised score')
#  geom_label_repel(data=filter(singtidy,HRD==1,set=='HOMOLOGOUS_RECOMBINATION'),aes(score,set,label=label),alpha=.8, 
#                   box.padding = 0.25, size=3, show.legend=F,
#                   point.padding = 0.25, segment.color='black') + ylab("")


c.source <- c("mediumseagreen","orange1")[as.factor(phen$source)]
c.erg <- c("grey","black")[as.factor(phen$host)]
c.type <- group.colors[as.factor(phen$type)]
mmr.type <- c("white","purple")[as.factor(phen$MMR)]
hr.type <-  c("white","blue")[as.factor(phen$HRD)]

clab <- matrix(cbind(c.source,c.erg,c.type),nrow=dim(phen)[1],ncol=3)
colnames(clab) <- c('source','state','pathology')
mypalette <- brewer.pal(9,"RdYlBu")
morecols <- colorRampPalette(mypalette)
png(filename = paste('plots/KEGG_DNA_REPAIR.png',sep="/"),
    width = 1240, height = 600, units = "px", pointsize = 10, bg = "white")
heatmap3(as.matrix(singdata),col=rev(morecols(25)),  distfun = dist, 
         main="", showColDendro = T,cexRow = 1,ColSideWidth=5,
         ColSideColors=clab,scale="row",method="complete",margins=c(7,0))
dev.off()


#neuroendocrine related signatures
gmt.1 <- getGmt("data/AR.NE.gene.sets.symbols.gmt")
gmt <- GeneSetCollection(c(gmt.1))
singscore.out <- run.singscore(dat.cpm,gmt)
singdata <- t(singscore.out %>% dplyr::select(-sample))
singtidy <- gather(singscore.out,key='set',value='score',-sample)
col.cell <- group.colors[as.factor(phen$type)]

#neuroendocrine related signatures
tmp.cpm <- dat.cpm[,phen$sample]
gmt.1 <- getGmt("data/AR_NED_sym.gmt")
gsets <- c('KAMMINGA_EZH2_TARGETS','BENPORATH_SOX2_TARGETS','HALLMARK_ANDROGEN_RESPONSE',
           'NELSON_RESPONSE_TO_ANDROGEN_UP','KEGG_MAPK_SIGNALING_PATHWAY','REACTOME_SIGNALING_BY_FGFR','HALLMARK_MYC_TARGETS_V2')
gmt <- read.msig(gsets)
gmt <- GeneSetCollection(c(gmt,gmt.1))
singscore.out <- run.singscore(tmp.cpm,gmt)
singdata <- t(singscore.out %>% dplyr::select(-sample))
singtidy <- gather(singscore.out,key='set',value='score',-sample)
col.cell <- group.colors[as.factor(phen$type)]
png(filename = paste('plots/MULTI_PRAD_SIG.png',sep="/"),
    width = 740, height = 624, units = "px", pointsize = 10, bg = "white")
heatmap.2(as.matrix(singdata),col=rev(morecols(25)),trace="none", main="Prostate Cancer Signatures",
          ColSideColors=col.cell,scale="row", keysize = 0.5, cexRow = 2, cexCol=0.5,
          Rowv=F,dendrogram = 'col',margins=c(10,40),labCol = FALSE)
dev.off()

gsets <- c('HALLMARK_WNT_BETA_CATENIN_SIGNALING','REACTOME_BETA_CATENIN_INDEPENDENT_WNT_SIGNALING',
           'KEGG_WNT_SIGNALING_PATHWAY','REACTOME_REPRESSION_OF_WNT_TARGET_GENES','REACTOME_TCF_DEPENDENT_SIGNALING_IN_RESPONSE_TO_WNT')
gmt <- read.msig(gsets)
gmt <- GeneSetCollection(c(gmt))
singscore.out <- run.singscore(dat.cpm[phen$sample],gmt)
singdata <- t(singscore.out %>% dplyr::select(-sample))
singtidy <- gather(singscore.out,key='set',value='score',-sample)

c.source <- c("mediumseagreen","orange1")[as.factor(phen$source)]
c.erg <- c("grey","black")[as.factor(phen$host)]
c.type <- group.colors[as.factor(phen$type)]

clab <- matrix(cbind(c.source,c.erg,c.type),nrow=dim(phen)[1],ncol=3)
colnames(clab) <- c('source','state','pathology')
mypalette <- brewer.pal(9,"RdYlBu")
morecols <- colorRampPalette(mypalette)
heatmap3(as.matrix(singdata),col=rev(morecols(25)),  distfun = dist, 
         main="", showColDendro = T,cexRow = 1,ColSideWidth=5,
         ColSideColors=clab,scale="none",method="complete",margins=c(7,0))

heatmap3(as.matrix(singdata[,phen[phen$paperlist=='Y',]$sample]),col=rev(morecols(25)),  distfun = dist, 
         main="", showColDendro = T,cexRow = 1,ColSideWidth=5,
         ColSideColors=clab,scale="row",method="complete",margins=c(7,0))

singtidy <- merge(singtidy,phen,by="sample")
wnt <- c()
for (i in unique(singtidy$set)) {
  tmp <- singtidy[singtidy$set==i,]
  print(t.test(tmp[tmp$pathology=='AR+/NE-',]$score,tmp[tmp$pathology=='AR-/NE+',]$score))
  wnt <- rbind(wnt,c(i,t.test(tmp[tmp$pathology=='AR+/NE-',]$score,tmp[tmp$pathology=='AR-/NE+',]$score)$p.value))
}
wnt <- as.data.frame(wnt)
names(wnt) <- c('set','p')
p.adjust(method="fdr",wnt$p)
#source.colors <- c("mediumseagreen","orange2")
#ggplot(singtidy) + geom_boxplot(aes(set,scale(score),fill=as.factor(source)),width=0.75) + 
#  coord_flip() + scale_fill_manual(values=source.colors) + theme_minimal() +
#  theme(legend.title=element_blank()) + xlab("") + ylab("")
#ggsave('plots/WNT-boxplot-source.png',width=10,height=6)
#
#ggplot(singtidy) + geom_boxplot(aes(set,scale(score),fill=as.factor(type)),width=0.75) + 
#  coord_flip() + scale_fill_manual(values=group.colors) + theme_minimal() +
#  theme(legend.title=element_blank()) + xlab("") + ylab("")
#ggsave('plots/WNT-boxplot-source.png',width=10,height=6)
#
#ggplot(singtidy) + geom_boxplot(aes(pathology,score),width=0.75) + 
#  coord_flip() + scale_fill_manual(values=group.colors) + theme_minimal() +
#  theme(legend.title=element_blank()) + xlab("") + ylab("")




####### hormone receptor pathways

gsets <- c(
  #'BIOCARTA_IGF1_PATHWAY',
  'GO_RESPONSE_TO_PROGESTERONE',
  'GO_RESPONSE_TO_ESTROGEN',
  'GO_RESPONSE_TO_CORTICOSTEROID',
  'GO_PROGESTERONE_METABOLIC_PROCESS',
  'GO_PROGESTERONE_RECEPTOR_SIGNALING_PATHWAY',
  #'GO_INSULIN_BINDING',
  #'GO_EPINEPHRINE_BINDING',
  'GO_GLUCOCORTICOID_SECRETION',
  #'GO_GONADOTROPIN_SECRETION',
  #'GO_ENDOCRINE_PROCESS',
  'GO_CELLULAR_RESPONSE_TO_HORMONE_STIMULUS',
  'GO_ESTROGEN_RECEPTOR_ACTIVITY')

gmt <- read.msig(gsets)
gmt <- GeneSetCollection(c(gmt))
singscore.out <- run.singscore(dat.cpm[phen[phen$paperlist=='Y',]$sample],gmt)
singdata <- t(singscore.out %>% dplyr::select(-sample))
singtidy <- gather(singscore.out,key='set',value='score',-sample)

png(filename = paste('plots/HORMONE-SIG.png',sep="/"),
    width = 1240, height = 600, units = "px", pointsize = 10, bg = "white")
#heatmap.2(as.matrix(singdata),col=rev(morecols(25)),  distfun = dist, ColSideColors=c.source,scale="row")
heatmap3(as.matrix(singdata),col=rev(morecols(25)),  distfun = dist, 
         main="", showColDendro = T,cexRow = 1,
         ColSideColors=clab,scale="none",method="complete",margins=c(7,1))
dev.off()

### hormone genes
png(filename = paste('plots/HORMONE-gene-unscaled.png',sep="/"),
    width = 1240, height = 600, units = "px", pointsize = 12, bg = "white")
heatmap3(as.matrix(dat.cpm[c('AR','PGR','NR3C2','NR3C1','ESR1','ESR2'),
                           phen[phen$paperlist=='Y',]$sample]),cexRow=3,cexCol=3,cex=0.5,
         ColSideColors=clab,col=rev(morecols(30)),
         scale="none",margins=c(5,5))
dev.off()

png(filename = paste('plots/HORMONE-gene-unscaled-set2.png',sep="/"),
    width = 1240, height = 600, units = "px", pointsize = 10, bg = "white")
heatmap3(as.matrix(dat.cpm[c('PGR','NR3C2','ESR1','ESR2'),
                           phen[phen$paperlist=='Y',]$sample]),cexRow=2,cex=0.5,
         ColSideColors=clab,col=rev(morecols(30)),
         scale="none",margins=c(5,5))
dev.off()

gsets <- c('BENPORATH_ES_1',
           'BENPORATH_ES_2',
           'BENPORATH_PRC2_TARGETS',
           'BENPORATH_NOS_TARGETS')

gmt <- read.msig(gsets)
gmt <- GeneSetCollection(c(gmt))
singscore.out <- run.singscore(dat.cpm[phen$sample],gmt)
singdata <- t(singscore.out %>% dplyr::select(-sample))
singtidy <- gather(singscore.out,key='set',value='score',-sample)
singtidy <- merge(singtidy,phen,by="sample")
#ggplot(singtidy %>% filter(set=='YEMELYANOV_GR_TARGETS_DN')) + geom_boxplot(aes(scale(score),set,color=as.factor(pathology),fill=as.factor(pathology)),width=0.4,alpha=0.25) +
#  theme_classic() + xlab('normalised score')

png(filename = paste('plots/STEM-SIG.png',sep="/"),
    width = 1240, height = 600, units = "px", pointsize = 10, bg = "white")
#heatmap.2(as.matrix(singdata),col=rev(morecols(25)),  distfun = dist, ColSideColors=c.source,scale="row")
heatmap3(as.matrix(singdata[,phen[phen$paperlist=='Y',]$sample]),col=rev(morecols(25)),  distfun = dist, 
         main="", showColDendro = T,cexRow = 1,ColSideWidth=5,
         ColSideColors=clab,scale="row",method="complete",margins=c(7,0))
dev.off()

###PAM50
################
x <- read.table('~/Dropbox/Projects__//data/PAM50_pdx/bioclassifier_example/pdx-cohort_pam50scores.txt',header=T,sep="\t")
x <- x[complete.cases(x),] %>% filter(Confidence>0.9)
mat <- as.data.frame(x %>% dplyr::select(X,Basal,Her2,LumA,LumB,Normal))
row.names(mat) <- mat$X
mat <- dplyr::select(mat,-X)
pam_phen <- phen.all[phen.all$sample %in% row.names(mat),]
mat <- mat[row.names(mat) %in% pam_phen$sample,]
#pam_phen <- pam_phen[pam_phen$paperlist=='Y',]

tmp <- t(mat)
tmp <- tmp[c('Basal','LumA','LumB'),]
tmp <- tmp[,colnames(tmp) %in% phen$sample]
#col.cell <- group.colors[as.factor(pam_phen[pam_phen$paperlist=='Y',]$type)]

c.source <- c("mediumseagreen","orange1")[as.factor(pam_phen$source)]
c.erg <- c("grey","black")[as.factor(pam_phen$state.x)]
c.type <- group.colors[as.factor(pam_phen$type)]

clab <- matrix(cbind(c.source,c.erg,c.type),nrow=dim(pam_phen)[1],ncol=3)
colnames(clab) <- c('source','state','pathology')
mypalette <- brewer.pal(9,"RdYlBu")
morecols <- colorRampPalette(mypalette)

png(filename = paste('plots/PAM50-PDX.png',sep="/"),width = 600, height = 340, units = "px", pointsize = 12, bg = "white",  res = NA )
heatmap.2(as.matrix(tmp),col=rev(morecols(50)),trace="none", main="PAM50",ColSideColors=col.cell,cexCol=1,cexRow=1.5,scale="none",margins=c(10,10))
dev.off()


#FOXA1 and AR - relationship
tmp <- as.data.frame(t(dat.cpm[c('FOXA1','AR'),]))
tmp$sample <- row.names(tmp)
ggplot(tmp[tmp$sample %in% phen[phen$paperlist=='Y',]$sample,]) + geom_point(aes(FOXA1,AR))




###### BELTRAN NE 

gsets <- c("HALLMARK_ANDROGEN_RESPONSE","HALLMARK_DNA_REPAIR","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
           "HALLMARK_MYC_TARGETS_V2","HALLMARK_TGF_BETA_SIGNALING",
           "HALLMARK_WNT_BETA_CATENIN_SIGNALING","KAMMINGA_EZH2_TARGETS","KEGG_MAPK_SIGNALING_PATHWAY",
           "REACTOME_SIGNALING_BY_FGFR","KEGG_HOMOLOGOUS_RECOMBINATION","KEGG_MISMATCH_REPAIR")
#gmt.1 <- getGmt("data/AR_NED_sym.gmt")
gmt <- read.msig(gsets)

gmt.1 <- getGmt("../chapter2/data/AR.NE.gene.sets.symbols.gmt")
gmt <- GeneSetCollection(c(gmt,gmt.1))
singscore.out <- run.singscore(dat.cpm,gmt)
singdata <- t(singscore.out %>% dplyr::select(-sample))
singtidy <- gather(singscore.out,key='set',value='score',-sample)

rownames(singdata) <- gsub("Beltran.up.NE.29","BELTRAN NEUROENDOCRINE UP",rownames(singdata))
ndx <- rownames(singdata) %in% c(gsets,"BELTRAN NEUROENDOCRINE UP")
rownames(singdata) <- gsub("HALLMARK_","HLMRK_",rownames(singdata))
rownames(singdata) <- gsub("EPITHELIAL_MESENCHYMAL_TRANSITION","EMT",rownames(singdata))
rownames(singdata) <- gsub("_"," ",rownames(singdata))
png(filename = paste('plots/EXTENDED-SIGNATURES-subset.pdf',sep="/"),
    width = 860, height = 720, units = "px", pointsize = 12, bg = "white")
heatmap.2(as.matrix(singdata[ndx,]),col=rev(morecols(50)),trace="none", main="MSigDB Extended Sets",
          ColSideColors=col.cell,scale="row", keysize = 1, cexRow = 1, cexCol=1,margins=c(10,20),#labCol = TRUE,
          Rowv=T,Colv=T,dendrogram="col",lhei = c(1.5,4))
dev.off()


singtidy <- dplyr::select(singscore.out,'sample','HALLMARK_ANDROGEN_RESPONSE','Beltran.up.NE.29')
tmp.ne <- merge(singtidy,phen,by="sample")
tmp.ne$shape <- paste0(tmp.ne$source,"-",tmp.ne$paperlist)
tmp.ne$shape2 <- paste0(tmp.ne$source,"-",tmp.ne$paperlist)
scale1 <- c(16,17)
ggplot(tmp.ne) + geom_point(aes(HALLMARK_ANDROGEN_RESPONSE,Beltran.up.NE.29,color=as.factor(type),shape=as.factor(shape)),size=2,alpha=0.75,stroke=1) +
  geom_label_repel(data=tmp.ne,aes(x=HALLMARK_ANDROGEN_RESPONSE,y=Beltran.up.NE.29,label=sample),alpha=0.85, 
                   box.padding = 0.25, size=2, show.legend=F,
                   point.padding = 0.25, segment.color='black') +
  theme_classic() + scale_color_manual(values=group.colors) + ylab("BELTRAN NEUROENDOCRINE UP") + xlab("HALLMARK ANDROGEN RESPONSE") + 
  scale_shape_manual(values=scale1) 
ggsave('plots/AR-V-NE-select-lines-label.png',width=12,height=9)
ggsave('plots/AR-V-NE-select-lines-label.pdf',width=12,height=9)

singtidy <- dplyr::select(singscore.out,'sample','HALLMARK_ANDROGEN_RESPONSE','Beltran.up.NE.29')
tmp.ne <- merge(singtidy,phen,by="sample")
tmp.ne$shape <- paste0(tmp.ne$source,"-",tmp.ne$state.y)
tmp.ne$shape2 <- paste0(tmp.ne$source,"-",tmp.ne$state.y)
scale1 <- c(1,16,2,17)
ggplot(tmp.ne) + geom_point(aes(HALLMARK_ANDROGEN_RESPONSE,Beltran.up.NE.29,color=as.factor(type),shape=as.factor(shape)),size=2,alpha=0.75,stroke=1) +
  #geom_label_repel(data=tmp.ne,aes(x=HALLMARK_ANDROGEN_RESPONSE,y=Beltran.up.NE.29,label=sample),alpha=0.85, 
  #                 box.padding = 0.25, size=2, show.legend=F,
  #                 point.padding = 0.25, segment.color='black') +
  theme_minimal() + scale_color_manual(values=group.colors) + ylab("BELTRAN NEUROENDOCRINE UP") + xlab("HALLMARK ANDROGEN RESPONSE") + 
  scale_shape_manual(values=scale1) 
ggsave('plots/AR-V-NE-select-lines.png',width=5,height=4)
ggsave('plots/AR-V-NE-select-lines.pdf',width=5,height=4)

singtidy <- dplyr::select(singscore.out,'sample','HALLMARK_ANDROGEN_RESPONSE','Beltran.up.NE.29')
tmp.ne <- merge(singtidy,phen,by="sample")
tmp.ne$shape <- paste0(tmp.ne$source,"-",tmp.ne$paperlist)
tmp.ne$shape2 <- paste0(tmp.ne$source,"-",tmp.ne$paperlist)
scale1 <- c(16,17)
ggplot(tmp.ne) + geom_point(aes(HALLMARK_ANDROGEN_RESPONSE,Beltran.up.NE.29,color=as.factor(type),shape=as.factor(shape)),size=3,alpha=0.75,stroke=1) +
  geom_label_repel(aes(HALLMARK_ANDROGEN_RESPONSE,Beltran.up.NE.29,label=sample),alpha=0.85, 
                   box.padding = 0.25, size=3, show.legend=F,
                   point.padding = 0.25, segment.color='black') +
  theme_minimal() + scale_color_manual(values=group.colors) + ylab("BELTRAN NEUROENDOCRINE UP") + xlab("HALLMARK ANDROGEN RESPONSE") + 
  scale_shape_manual(values=scale1) 
ggsave('plots/AR-V-NE-select-lines-labelled.pdf',width=10,height=8)
