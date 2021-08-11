##preload data

library(org.Hs.eg.db)
columns(org.Hs.eg.db)
gene_names <- AnnotationDbi::select(org.Hs.eg.db,keys=(keys(org.Hs.eg.db)),columns=c("SYMBOL"))
names(gene_names) <- c("GeneID","Symbol")

#gene_names <- read.table('data/Homo_sapiens.gene_info.simple.txt',header=T,sep="\t")

run.singscore <- function(counts.log,gene.set) {
  rankData <- rankGenes(counts.log)
  score_dat <- c()
  tmp <- c()
  for (gs in gene.set) {
    up.set <- geneIds(gs)[geneIds(gs) %in% row.names(counts.log)]
    tmp_score <- simpleScore(rankData, upSet = up.set)
    tmp_score$set <- setName(gs)
    tmp_score$sample <- row.names(tmp_score)
    tmp <- rbind(tmp,tmp_score)
    tmp_score <- tmp_score[,c("TotalScore","sample")]
    names(tmp_score) <- c(setName(gs),"sample")
    if (!is.null(score_dat)) score_dat <- cbind(score_dat,tmp_score[setName(gs)])
    if (is.null(score_dat)) score_dat <- tmp_score[,c(2,1)] 
  }
  return(score_dat)
}

get.int.load  <- function(pca,n=1) {
  loading_scores <- pca$rotation[,n]
  gene_scores <- abs(loading_scores) ## get the magnitudes
  gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
  top_10_genes <- names(gene_score_ranked[1:10])
  filter(gene_names,GeneID %in% top_10_genes)
  gscores <- as.data.frame(gene_scores)
  gscores$GeneID <- row.names(gscores)
  int.load <- filter(gene_names,GeneID %in% names(boxplot(gene_scores,plot=FALSE)$out))
  int.load <- merge(int.load,gscores,by="GeneID")
  return(arrange(int.load,-gene_scores))
}

ggmds <- function(dge,sample_info,name='all',repel=F) {
  out <- plotMDS(dge)
  out.mds <- as.data.frame(cbind(out[6]$x,out[7]$y))
  out.mds$sample <- row.names(out.mds)
  out.mds <- merge(out.mds,sample_info,by="sample")
  g <- ggplot(dat=out.mds) + geom_label(aes(x=V1,y=V2,color=type,label=sample)) +
    theme_classic(base_size = 12) + theme(legend.position="right") +
    theme(legend.title=element_blank())  
  if (repel) {
    g <- g + geom_label_repel(aes(x=V1,y=V2,label=sample), box.padding = 0.35, 
                     fontface = 'bold', color = 'black', size=3, show.legend=F,
                     point.padding = 0.5, segment.color='grey50') 
    }
    ggsave(paste0('plots/mds-',name,'-label.png'),width=8,height=6)
    return(out.mds)
}

ggpca_line <- function(dge,samples,sample_info,name) {
  log.cpm <- cpm(dge,log=T)
  log.cpm <- log.cpm[,colnames(log.cpm) %in% samples]
  print("running pca")
  pca <- prcomp(t(log.cpm))
  pca.tmp <- as.data.frame(pca$x)
  pca.tmp$sample <- row.names(pca$x)
  pca.tmp <- merge(pca.tmp,sample_info,by="sample")
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
  
  pca.tmp$PC2 <- 1
  print("plotting")
  ggplot(dat=pca.tmp) +
    geom_point(aes(x=PC1,y=PC2,color=as.factor(type))) +
    xlab(paste0("PC1: var explained: ",format(pca.var.per[1],digits=3),"%")) + 
    ylab(paste0("PC2: var explained: ",format(pca.var.per[2],digits=3),"%")) +
    theme_classic(base_size = 12) + theme(legend.position="right") +
    theme(legend.title=element_blank()) 
  ggsave(paste0('plots/PC1-',name,'label.png'),width=8,height=6)
  return(pca)
  
}

ggpca <- function(dge,samples,sample_info,name,width=10,height=4) {
  log.cpm <- cpm(dge,log=T)
  log.cpm <- log.cpm[,colnames(log.cpm) %in% samples]
  print("running pca")
  pca <- prcomp(t(log.cpm))
  pca.tmp <- as.data.frame(pca$x)
  pca.tmp$sample <- row.names(pca$x)
  pca.tmp <- merge(pca.tmp,sample_info,by="sample")
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
  
  print("plotting")
  ggplot(dat=pca.tmp) +
    #geom_label(aes(PC1,PC2,label=sample,color=as.factor(type),alpha=0.5),show.legend=F) +
    geom_label_repel(aes(x=PC1,y=PC2,label=sample), 
                     box.padding = 0.25, size=4, show.legend=F,
                     point.padding = 0.25, segment.color='black') +
    geom_point(aes(x=PC1,y=PC2,color=as.factor(type),size=4)) +
    xlab(paste0("PC1: var explained: ",format(pca.var.per[1],digits=3),"%")) + 
    ylab(paste0("PC2: var explained: ",format(pca.var.per[2],digits=3),"%")) +
    theme_classic(base_size = 16) + theme(legend.position="right") +
    theme(legend.title=element_blank()) 
  ggsave(paste0('plots/PCA-',name,'label.png'),width=width,height=height)
  return(pca)
  
}

##useful functions
plot.cpm <- function(dge,file,title) {
  counts.log <- cpm(dge,log=TRUE)
  png(file,width=1000,height=400)
  boxplot(counts.log, xlab="", ylab="Log2 counts per million",las=2)
  abline(h=median(counts.log),col="blue")
  title(title)
  dev.off()
}

plot.libsize <- function(dge,filename) {
  png(filename,width=1400,height=400)
  barplot(dge$samples$lib.size,names=colnames(dge),las=2)
  title("Library sizes")
  dev.off() 
}

filterRNASEQ <- function(data,threshold,cpm_min=1) {
  counts <- reshape2::dcast(data, GeneID  ~ sample,value.var="counts")
  row.names(counts) <- counts$GeneID
  counts <- counts[,2:dim(counts)[2]]
  cpm <- reshape2::dcast(data, GeneID  ~ sample,value.var="cpm")
  row.names(cpm) <- cpm$GeneID
  cpm <- cpm[,2:dim(cpm)[2]]
  #cpm_min <- 1
  thresh <- cpm > cpm_min
  table(rowSums(thresh))
  keep <- rowSums(thresh) >= threshold
  counts.keep <- counts[keep,]
  print("Reads kept and filtered")
  print(summary(keep))
  return(counts.keep)
}



##featurecount multiple sam files
featureCounts.multi <- function(path,pattern) {
  data <- c() 
  rna.files <- list.files(path,pattern)
  for (file in rna.files) {
    fc <- featureCounts(paste0(path,'/',file),annot.ext=ann,isPairedEnd=FALSE)
    out <- cbind(fc$counts,fc$annotation[,c("GeneID","Length")])
    names(out)[1] <- "counts"
    out$sample <- file
    out$cpm <- as.vector(cpm(DGEList(counts=out$counts, genes=out$annotation[,c("GeneID","Length")])))
    out <- merge(out,gene_names,by="GeneID",all.x=T)
    data <- rbind(data,out)
  }
  return(data)
}

add_symbol <- function(data,gene_names){
  data <- as.data.frame(data)
  data$GeneID <- row.names(data)
  names(data) <- gsub("\\.","",names(data))
  return(merge(data,gene_names,by="GeneID",sort=F))
}

ggplot_de <- function(data,n=0,pthreshold=1e-6,title='DE plot') {
  if (n>0) pthreshold <- data[n,]$adjPVal 
  g <- ggplot(data) + 
    geom_point(aes(x=logFC,y=-log10(adjPVal))) +
    geom_point(data=filter(data,adjPVal<pthreshold),aes(x=logFC,y=-log10(adjPVal)),color="red") +
    geom_label_repel(data=filter(data,adjPVal<pthreshold),aes(x=logFC,y=-log10(adjPVal),label=Symbol), 
                     box.padding = 0.25, size=3, show.legend=F,
                     point.padding = 0.25, segment.color='black') +
    theme_classic(base_size = 12) + theme(legend.position="right") +
    theme(legend.title=element_blank()) +
    ggtitle(title)
  return(g)
}

ggplot_cpm <- function(data,sam=c(),sym=c(),title="Select genes",cols=1) {
  subset.rna <-  data.rna %>% 
    #dplyr::select(sample,Symbol,cpm,Disease.stage,responsive) %>% 
    dplyr::select(sample,Symbol,cpm,type) %>% 
    filter(Symbol %in% sym) %>% 
    filter(sample %in% sam)
  g <- ggplot(subset.rna) + 
    #geom_bar(aes(sample,cpm,fill=as.factor(responsive)),stat="identity") + 
    geom_bar(aes(sample,cpm,fill=as.factor(type)),stat="identity") + 
    facet_wrap(~Symbol,scales="free_y",ncol=cols) +
    theme(axis.text.x=element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("") +
    guides(fill=guide_legend(title="Treatment Sensitive")) 
  return(g)
}

hypergeom_msigdb <- function(mset,data,nlength) {
    mdbset <- mset %>% 
              dplyr::select(gs_name,gene_symbol) %>%
              rename(gene_symbol='Symbol') %>% 
              as.data.frame()
    out <- c()
    for (gs in unique(mset$gs_name)) {
      mdbset <- mset %>% 
        dplyr::filter(gs_name==gs) %>% 
        dplyr::select(gs_name,gene_symbol) %>%
        rename(gene_symbol='Symbol') %>%
        as.data.frame()
      q<-length(filter(data,Symbol %in% mdbset$Symbol)$Symbol)
      m<-length(mdbset$Symbol)
      n<-nlength-m
      k<-length(data$Symbol)
      out <- rbind(out,as.data.frame(cbind(gs,phyper(q-1,m,n,k,lower.tail=F))))
    }
    out$V2 <- as.numeric(as.character(out$V2))
    out$p.adjust <- p.adjust(p=out$V2,method="BH")
    return(out)
}



  
read.multi <- function(loc="./",pat=NULL,sep="\t",header=T,quote="\"",stringsAsFactors=F) {
  files <- as.data.frame(list.files(pattern=pat,path=loc))
  data <- c()
  
  for (f in t(files)) {
    pass = list.files(path=paste0(loc,f))#,pattern="*_calls.vcf$")
    fname <- paste0(loc,f)
    cl <- read.table(fname,sep=sep,header=header,quote=quote,stringsAsFactors=stringsAsFactors)
    if (nrow(cl)>0) { 
    cl$name <- f
    data <- rbind(data,cl)
    }
  }
  return(data)
}

read.msig <- function(gsets) {
  mset = msigdbr(species = "Homo sapiens")
  gmt.set <- c()
  for (gs in gsets) {
    tmp <- mset %>% filter(gs_name==gs)
    tmp <- GeneSet(tmp$human_gene_symbol,setIdentifier = gs)
    setName(tmp) <- gs
    gmt.set <- c(gmt.set,tmp)
  }
  return(gmt.set)
}

#take list of genes, return gene set object with one gene for each
single.genes.to.sets <- function(genes) {
  gmt.set <- c()
  for (g in genes) {
    gmt.set <- c(gmt.set,GeneSet(g,setIdentifier = g,setName=g))
  }
  return(gmt.set)
}

