#library('VariantAnnotation')
#library('StructuralVariantAnnotation')
library('dplyr')
library('readr')
library('GenomicRanges')
library('ggplot2')
library('stringr')
library('tidyr')
library('magrittr')
#library('skimr',lib="~/lib")


## note to self - really useful
#tmp %>% filter(grepl("APOBEC", Symbol))

chr.offset <- function()
{
    tot <- read.table('~/res/chr-hg19.txt',header=T)

    tot$total <- as.numeric(tot$end)
    for (i in 2:nrow(tot)) {
            tot[i,]$start <- tot[i-1,]$end + 1 
        tot[i,]$end <- tot[i,]$start + tot[i,]$end
    }   
    return(tot)
}


#read many functions, put them into genomic ranges for overlapping
range.as.data.frame <- function(range) 
{
    out <- as.data.frame(range)
    names(out) <- gsub("mcols.","",names(out))
    names(out)[1] <- "chr"
    return(out)
}

read.cnvkit <- function(loc=NULL,files=NULL,as.range=TRUE,xnum=F)
{
    cnvkit <- c()
    if (is.null(loc)) loc <- "/researchers/david.goode/Analysis/170224_prostate_wgs/cnvkit/"
    if (is.null(files)) {
        files <- as.data.frame(list.files(pattern="*call.cns$",path=loc)) 
    } else { files <- as.data.frame(list.files(pattern=files,path=loc)) }
    names(files) <- "ID"
    for (i in files$ID) {
        #fname <- (list.files(pattern="_merged.cns",path=paste0(loc,i)))[1]
        x <- read.table(paste0(loc,"/",i),header=T,check.names=F,stringsAsFactors=F)
        x$ID <- i
        cnvkit <- rbind(cnvkit,x)
    }
    if (xnum) cnvkit$chromosome <- gsub("Y",24,gsub("X",23,cnvkit$chromosome))
    if (as.range) cnvkit <- GRanges(cnvkit$chromosome,IRanges(cnvkit$start,cnvkit$end),mcols=cnvkit[,c(5:8,10)])
    return(cnvkit)
}



read.multi <- function(loc="./",pat=NULL,sep="\t",header=T,quote="",stringsAsFactors=F) {

    #if (is.null(loc)) loc <- "."
    #if (is.null(pat)) pat <- "" else print("try s <- read.somatic.snps(pat='--FR07888586$')")
    files <- as.data.frame(list.files(pattern=pat,path=loc))
    data <- c()

    for (f in t(files)) {
        pass = list.files(path=paste0(loc,f))#,pattern="*_calls.vcf$")
        fname <- paste0(loc,f)
        cl <- read.table(fname,sep=sep,header=header,quote=quote,stringsAsFactors=stringsAsFactors)
        #names(cl) <- c("CHROM","POS","ID","REF","INFO","QUAL","TYPE","FORMAT")
        cl$name <- f
        data <- rbind(data,cl)
    }
    return(data)
}


read.gridss <- function(loc=NULL,files=NULL,pat=NULL,incfail=FALSE,as.range=T) {
    if (is.null(loc)) loc <- "/researchers/david.goode/Analysis/170224_prostate_wgs/Gridss/"
    if (is.null(pat)) pat <- print("try s <- read.gridss(pat='--FR07888586$')") 
    if (is.null(files)) files <- as.data.frame(list.files(pattern=pat,path=loc))

    data <- c()

    print('reading ')
    print(list.files(path=loc,pat=pat))
    for (i in list.files(path='../Gridss/',pat=pat)) {
        tmp <- read.table(paste0('../Gridss/',i,'/SV/',i,'_gridss.vcf'),header=T,sep="\t")
        tmp$name <- i
        data <- rbind(data,tmp)
    }

    max_struct_qual <- max(data$QUAL)
    data$scale_qual <- data$QUAL / max_struct_qual

    return(data)

}

read.clove <- function(loc=NULL,files=NULL,pat=NULL,incfail=FALSE,as.range=T) {

    if (is.null(loc)) loc <- "/researchers/david.goode/Analysis/170224_prostate_wgs/Clove/"
    if (is.null(pat)) pat <- "^FR*" else print("try s <- read.clove(pat='--FR07888586$')")
    if (is.null(files)) files <- as.data.frame(list.files(pattern=pat,path=loc))
    clove <- c()

    for (f in t(files)) {
        print(f)
        if (incfail) {
            pass = list.files(path=paste0(loc,f,"/"),pattern="*_calls.vcf$")
        } else { 
            pass = list.files(path=paste0(loc,f,"/"),pattern="*.pass.vcf$")
            print(pass)
        }
        fname <- paste0(loc,f,"/",pass[1])
        cl <- read.table(fname)
        names(cl) <- c("CHROM","POS","ID","REF","INFO","QUAL","TYPE","FORMAT")
        cl$name <- f
        clove <- rbind(clove,cl)
    }
    #if (as.range) {
    #    clove <- GRanges(clove$CHROM,IRanges(clove$start,clove$end), mcols=clove[c("freq","cna","depth","name")])
    #}
    return(clove)
}

read.battenberg <- function(loc=NULL,files=NULL,as.range=T,pat=NULL) {

    results <- c()
    if (is.null(loc)) loc <- "/researchers/david.goode/Analysis/170224_prostate_wgs/Battenberg/"
    if (is.null(pat)) pat <- "^FR*" else print("try s <- read.somatic.snps(pat='--FR07888586$')")
    if (is.null(files)) files <- as.data.frame(list.files(pattern=pat,path=loc))
    print(files)
    names(files) <- "ID"
    for (i in files$ID) {
        ii <-  as.character(str_split_fixed(i,"-",2)[1])
        if (file.exists(paste0(loc,i,"/",ii,"_battenberg_cn.vcf.gz")))
        {
            print(paste0(loc,i,"/",ii,"_battenberg_cn.vcf.gz"))
            x <- read.table(gzfile(paste0(loc,i,"/",ii,"_battenberg_cn.vcf.gz")),stringsAsFactors=F)
            #x <- read.table(gzfile("FR07888618_battenberg_cn.vcf.gz"),stringsAsFactors=F)
            names(x) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","NORMAL","TUMOUR")
            x$CHROM <- gsub("X",23,x$CHROM)

            cnv <- as_tibble(str_split_fixed(x$TUMOUR,":",n=7))
            names(cnv) <- c("Genotype","Total copy number","Minor allele copy number",
                             "Fraction Cells first state","Total copy number second state",
                             "Minor allele copy number second state","Fraction Cells second state")
            cnv[,2:4] <- lapply(cnv[,2:4], function(x) as.numeric(as.character(x)))
            cnv[,5:7] <- lapply(cnv[,5:7], function(x) as.numeric(as.character(x)))

            pos <- as_tibble(str_split_fixed(x$INFO,"=",n=3))
            pos$V3 <- as.integer(pos$V3)
            names(pos)[3] <- "end"
            x <- x[,c(1:5,8)]
            x <- cbind(x,pos[,3])
            x <- cbind(x,cnv)
            x$id <- i
        }
        results <- rbind(results,x)
    }
    if (as.range) results <- GRanges(results$CHROM,IRanges(results$POS,results$end),mcols=results[,3:15])
    return(results)
}

read.vcf.seg <- function(vcf,chr,start,end) {
   #need to lock file or create random filename then delete it!
    ucode <- runif(1, min = 0, max = 10) 
    vcf_tmp <- paste(vcf,".tmp.",ucode,sep="")
    system(paste("head -n 1 ",vcf," > ",vcf_tmp,sep=""))
    system(paste("cat ",vcf," | awk ' $1==",chr," && $2 > ",start-1, " && $2 < ",end-1,"' >> ",vcf_tmp,sep=""))
    #vcf <- read.table(vcf_tmp,header=T,sep="\t",stringsAsFactors=F)
    vcf <- read_table2(vcf_tmp)
    system(paste0('rm ',vcf_tmp))
    return(vcf)
}

bp2mbp <- function(bp) {
    return(bp /1000000)
}
mbp2bp <- function(mbp) {
    return(mbp * 1000000)
}

read.somatic.snps <- function(loc=NULL,files=NULL,as.range=TRUE,var='mut',pat=NULL,filter.canonical=T,all.cols=F)
{
    som <- c()
    if (is.null(loc)) loc <- "/researchers/david.goode/Analysis/170224_prostate_wgs/Stats/"
    if (is.null(pat)) pat <- "^CPCG*" else print("try s <- read.somatic.snps(pat='--FR07888586$')")
    if (is.null(files)) files <- as.data.frame(list.files(pattern=pat,path=loc))
    print(files)
    names(files) <- "ID"
    for (i in files$ID) {
        if (file.exists(paste0(loc,i,"/",i,".",var,".tsv"))) 
        {
            x <- read.table(paste0(loc,i,"/",i,".",var,".tsv"),header=T,
                            check.names=F,sep="\t",quote="\"",stringsAsFactors=F)
            print(paste0(loc,i,"/",i,".",var,".tsv"))
        }

        #MUTECT2 fix
        x[x$TYPE==".",]$TYPE <- x[x$TYPE==".",]$Variant_Type
        x$TYPE <- gsub("deletion","Deletion",x$TYPE)
        x$TYPE <- gsub("insertion","Insertion",x$TYPE)
        #x$ID <- i
        #print(head(x)) #include depth here
        if (filter.canonical) x <- filter(x,CANONICAL=='YES')
        if (all.cols==T) {
            som <- rbind(som,x)
        } else { 
            x <- x[c("CHROM","POS","Gene","Consequence_Rank","Identified","ID","REF",
                     "ALT","TYPE","SYMBOL","sample_name","tumour_rd","tumour_ad","Variant_Type","Consequence")]
            som <- rbind(som,x)
        }
    }
    if (as.range) som <- GRanges(som$CHROM,IRanges(som$POS-1,som$POS),mcols=som[,c("Gene","SYMBOL","Consequence_Rank","Identified","ID","REF","ALT","TYPE","sample_name","tumour_rd","tumour_ad","Variant_Type","Consequence")])
    return(som)
}

read.freec.pvalues <- function(loc=NULL,files=NULL,as.range=TRUE,subclonal=TRUE)
{
    info <- c()
    if (is.null(loc)) loc <- "/researchers/david.goode/Analysis/170224_prostate_wgs/phylo-v2/pvals"
    if (is.null(files)) {
        files <- as.data.frame(list.files(pattern="FR*",path=loc)) 
    } else { files <- as.data.frame(list.files(pattern=files,path=loc)) }
    names(files) <- "ID"
    print("reading")
    print(files)
    for (i in files$ID) {
        #fname <- (list.files(pattern="_CNVs",path=paste0(loc,i,"/")))
        #    for (f in fname) {
                x <- read.table(paste0(loc,'/',i),header=T,check.names=F,sep='\t')
                #x <- read.table(paste0(loc,i,"/",f),header=F,check.names=F)
                x$ID <- i
                info <- rbind(info,x)
        #}#
    }
    #chr    start   end copy number status  WilcoxonRankSumTestPvalue   KolmogorovSmirnovPvalue
    names(info) <- c("chrom","start","end","cnv","type","WilcoxonRankSumTestPvalue","KolmogorovSmirnovPvalue","ID")
    info$chrom <- gsub("X",23,info$chrom)
    info$chrom <- gsub("Y",24,info$chrom)
    info <- cbind(info,str_split_fixed(info$ID,".sim.",n=2))
    names(info)[9:10] <- c("f1","sim")
    if (as.range) info <- GRanges(info$chrom,
                                  IRanges(info$start,info$end),
                                  mcols=info[c("cnv","type","ID","sim","WilcoxonRankSumTestPvalue","KolmogorovSmirnovPvalue")])
    return(info)
}



read.freec <- function(loc=NULL,files=NULL,as.range=TRUE,subclonal=TRUE)
{
    info <- c()
    if (is.null(loc)) loc <- "/researchers/david.goode/Analysis/170224_prostate_wgs/phylo-v2/cna"
    if (is.null(files)) {
        files <- as.data.frame(list.files(pattern="FR*",path=loc)) 
    } else { files <- as.data.frame(list.files(pattern=files,path=loc)) }
    names(files) <- "ID"
    print("reading")
    print(files)
    for (i in files$ID) {
        #fname <- (list.files(pattern="_CNVs",path=paste0(loc,i,"/")))
        #    for (f in fname) {
                x <- read.table(paste0(loc,'/',i),header=F,check.names=F)
                #x <- read.table(paste0(loc,i,"/",f),header=F,check.names=F)
                x$ID <- i
                info <- rbind(info,x)
        #}#
    }
    names(info) <- c("chrom","start","end","cnv","type","ID")
    info$chrom <- gsub("X",23,info$chrom)
    info$chrom <- gsub("Y",24,info$chrom)
    info <- cbind(info,str_split_fixed(str_split_fixed(info$ID,".sim.",n=2)[,2],'[.]',n=4))
    #info <- cbind(info,str_split_fixed(info$ID,".sim.",n=2))
    names(info)[7:9] <- c("f1","sim","win")
    if (as.range) info <- GRanges(info$chrom,IRanges(info$start,info$end),mcols=info[c("cnv","type","ID","f1","sim","win")])
    return(info)
}


read.facets <- function(loc=NULL,files=NULL,as.range=TRUE,subclonal=TRUE)
{
#eg usage facets <- read.facets(files="^FR*")
    facet <- c()
    if (is.null(loc)) loc <- "/researchers/david.goode/Analysis/170224_prostate_wgs/Facets/"
    if (is.null(files)) {
        files <- as.data.frame(list.files(pattern="^CPCG*",path=loc)) 
    } else { files <- as.data.frame(list.files(pattern=files,path=loc)) }
    names(files) <- "ID"
    for (i in files$ID) {
        if (subclonal) fname <- (list.files(pattern="*.subclonal.WGS.tsv",path=paste0(loc,i)))[1]
        if (!subclonal) fname <- (list.files(pattern="*.clonal.WGS.tsv",path=paste0(loc,i)))[1]
        x <- read.table(paste0(loc,i,"/",fname),header=T,check.names=F,stringsAsFactors=F)
        if (is.null(x$clonal.cluster)) x$clonal.cluster <- -1 #-1 indicates missing field
        x$ID <- i
        facet <- rbind(facet,x)
    }
    if (as.range) facet <- GRanges(facet$chrom,IRanges(facet$start,facet$end),mcols=facet[,c(3,4,12,13,14,15,16)])
    return(facet)
}

read.blacklist <- function(loc=NULL,as.range=TRUE)
{
    if (is.null(loc)) loc <- "~/res/Scheinin_et_al_Blacklist.txt"
    print(paste0('reading ',loc))
    blacklist <- read.table(loc)
    names(blacklist) <- c("chr","start","end","width")
    blacklist$chr <- gsub("X",23,blacklist$chr)
    blacklist <- filter(blacklist,chr %in% c(1:23))
    blacklist$chr <- as.numeric(blacklist$chr)
    if (as.range) blacklist <- GRanges(blacklist$chr,IRanges(blacklist$start,blacklist$end))
    return(blacklist)
}

read.genelist <- function(loc=NULL,as.range=TRUE)
{
    if (is.null(loc)) loc <- '~/res/glist-hg19'
    hg19 <- read.table(loc,header=F)
    names(hg19) <- c("chr","start","end","gene")
    hg19$chr <- gsub("X",23,hg19$chr)
    genelist <- filter(hg19,chr %in% 1:23)
    if (as.range) genelist <- GRanges(genelist$chr,IRanges(genelist$start,genelist$end),mcols=genelist[c("gene")])
    return(genelist)

}

read.gaps <- function(loc=NULL,as.range=TRUE)
{
    if (is.null(loc)) loc <- '~/res/hg19-gaps-ucsc.txt'
    gaps <- read.table(loc,header=F)
    gaps <- gaps[,c(2,3,4,8,9)]
    names(gaps) <- c("chr","start","end","type","detail")
    gaps$chr <- gsub('chr','',gaps$chr)
    if (as.range) gaps <- GRanges(gaps$chr,IRanges(gaps$start,gaps$end),mcols=gaps[c("type","detail")])
    return(gaps)
}

read.exon <- function(loc=NULL, as.range=TRUE) {
    if (is.null(loc)) loc <- '~/res/UCSC_exons_hg19.bed'
    exon <- read.table(loc,header=T)
    print(paste0('reading exon file ',loc))
    names(exon) <- c("chr","start","end","exid","ii","strand") 
    exon$chr <- gsub("chr","",exon$chr)
    exon <- filter(exon,chr %in% c(1:23))
    if (as.range) exon <- GRanges(exon$chr,IRanges(exon$start,exon$end),mcols=exon[c("exid")])
    return(exon)
}

#read sureselect file
read.bed <- function(loc=NULL, as.range=TRUE) {
    if (is.null(loc)) loc <- '/data/reference/bed_files/AgilentSureSelect_SUPER/SureSelect_SUPER_0628271_Covered_mod_MERGED.bed'
    bed <- read.table(loc,header=T)
    print(paste0('reading bed file ',loc))
    bed <- bed[,1:3]
    names(bed) <- c("chr","start","end")
    bed$chr <- gsub("chr","",bed$chr)
    bed <- filter(bed,chr %in% c(1:23))
    if (as.range) bed <- GRanges(bed$chr,IRanges(bed$start,bed$end))
    return(bed)
}



facet.to.range <- function(facet) 
{
    facet <- GRanges(facet$chrom,IRanges(facet$start,facet$end),mcols=facet[,c(3,4,12,13,14,15,16)])
#    names(data) <- gsub("chrom|CHR|chromosome","chr",names(data))
    return(facet)
}

facet.from.range <- function(facet) 
{
    facet <- as.data.frame(facet)
    names(facet) <- gsub("mcols.","",names(facet))
    names(facet) <- gsub("seqnames","chrom",names(facet))
    return(facet)
}


segment.germline.freq <- function(germline,WINSIZE=500000,plot=TRUE) {
    all_seg_means_by_sample <- c()
    for (chr in c(1:23)) {
        #read blacklist by chr
        print(chr)
        grm_chr <- filter(germline,CHROM==chr)
        seg_means_by_sample <- c()
        for (id in unique(grm_chr$Sample)) {
            grm_chr_sam <- filter(grm_chr, Sample==id)
            grm_hom <- filter(grm_chr_sam, freq==1)
            grm_het <- filter(grm_chr_sam, freq<1)
            grm_het[grm_het$freq > 0.5,]$freq <- 1 - grm_het[grm_het$freq > 0.5,]$freq
            grm_all <- rbind(grm_het,grm_hom) #adding backing homozygous reads for depth
            #grm_all <- rbind(grm_het)
            MAXPOS <- range(grm_all$POS)[2]
            MINPOS <- range(grm_all$POS)[1]

            seg_means <- c()
            prv <- 0
            for (pos in seq(MINPOS,MAXPOS,WINSIZE)) {
                seg_all <- grm_all[grm_all$POS >= prv & grm_all$POS < pos,]
                seg <- grm_het[grm_het$POS >= prv & grm_het$POS < pos,]
                seg_mean <- c(mean(seg$freq))
                seg_depth <- c(mean(as.numeric(seg_all$tumour_rd + seg_all$tumour_ad)))
                if (nrow(seg)<5) {seg_mean <- 0} ###5 is a rounding error... just a guess
                if (is.finite(seg_depth)) {
                    if (seg_depth > 100) {seg_depth <- 100}
                }
                if (nrow(seg)<5) {seg_mean <- 0} ###10 is a rounding error... just a guess
                tmp_mean <- c(seg_mean,(pos-(WINSIZE/2)),prv,pos,seg_depth) #WINSIZE/2 is midpoint
                seg_means <- rbind(seg_means,tmp_mean)
                prv <- pos
            }

            seg_means <- as.data.frame(seg_means)
            names(seg_means) <- c("freq","POS","start","end","depth")
            seg_means$name <- id
            seg_means$chr <- chr
            seg_means <- seg_means[is.finite(seg_means$depth),]
            seg_means$sc_depth <- seg_means$depth / 100 #/max(seg_means_by_sample$depth)
            seg_means_by_sample <- rbind(seg_means_by_sample,seg_means)
            
        }

        facets <- filter(read.facets(as.range=F),chrom==chr)
        facets$cna <- 0
        facets[facets$tcn.em>2,]$cna <- 1
        facets[facets$tcn.em<2,]$cna <- -1
        blacklist <- read.blacklist(as.range=F)
        blacklist <- blacklist[blacklist$chr==chr,]

        if (plot==T) plot.cna(grm_chr,facets,blacklist,seg_means_by_sample,chr)

        all_seg_means_by_sample <- rbind(all_seg_means_by_sample,seg_means_by_sample)
        }
        write.table(all_seg_means_by_sample,'~/data/CPCG_all_seg_means_by_sample.txt',quote=F,row.names=F,col.names=T)

}

segment.to.cna <- function(segment,as.range=T,read_seg=F) {
    if (read_seg==T) segment <- as.data.frame(read_table2('~/data/CPCG_all_seg_means_by_sample.txt'))
    all <- segment
    all$cna <- 0
    all[all$freq < 0.37,]$cna <- 1
    all[all$freq < 0.275,]$cna <- 2
    all[all$freq < 0.01,]$cna <- -1

    all$chr <- as.numeric(all$chr)
    #merge into contiguous blocks of cna change
    newset <- c() 
    for (sample in unique(all$name)) {
        tmp <- filter(all, name==sample)
        tmp <- tmp[order(tmp$chr,tmp$start),]
        prv <- tmp[1,]
        for (row in 2:nrow(tmp)) {
            cur <- tmp[row,]
            if ((prv$cna==cur$cna) & (prv$chr==cur$chr)) {
                prv$end <- cur$end
            } else {
                newset <- rbind(newset,prv)
                prv <- cur 
            }   
        }   
        if (prv$cna==cur$cna) { newset <- rbind(newset,prv) }
    }
    segment <- newset
    if (as.range==T) {
        segment <- GRanges(segment$chr,IRanges(segment$start,segment$end), mcols=segment[c("freq","cna","depth","name")])
    }
    return(segment)
}


plot.cna <- function(data,facets,blacklist,seg_means_by_sample,chr) {
        seg_means_by_sample$freq <- seg_means_by_sample$freq * 2
        names(seg_means_by_sample) <- gsub("name","Sample",names(seg_means_by_sample))

        loc <- '~/data/pairs.grm.txt'
        y <- read.table(loc,header=T,sep=" ",stringsAsFactors=F)
        y <- filter(y, type != "Gross")
        names(y)[1] <- "Sample"
        names(y)[2] <- "ID"
        facets <- merge(facets,y,by="ID")
        print('doing a facets')
        print(unique(seg_means_by_sample$Sample))
        facets <- facets[facets$Sample %in% unique(seg_means_by_sample$Sample),]
        blacklist_chr <- blacklist[blacklist$chr==chr,]

        names(data) <- gsub("ID","UID",names(data))
        data <- merge(data,y,by="Sample")
        print(head(data))
        print(head(seg_means_by_sample))

        ggplot(data, aes(x=POS,y=freq)) + 
        geom_point(aes(alpha=0.2,color=factor(type))) + 
        geom_point(data=seg_means_by_sample,
                     aes(x=POS,y=freq,alpha=1),color="#99AADD") + 
        geom_segment(data=filter(blacklist_chr,chr==chr),
                     aes(x = start, xend = end, y=0.1, yend=0.1)) +
        geom_segment(data=facets,aes(x = start, xend = end, y=0.2+(clonal.cluster/10), 
                     yend=0.2+(clonal.cluster/10),color=factor(cna),
                     linetype=factor(clonal.cluster)),size=1.5) +
        geom_segment(data=seg_means_by_sample, 
                     aes(x=start,xend=end,y=freq,yend=freq,
                     alpha=1),color="#CC2288") + 
        ggtitle(paste0(chr," - VAF - germline")) + facet_wrap(~Sample,ncol=1) 

        print('about to plot')
        #plot_height <- length(unique(data$Sample)) * 3.6
        plot_height <- 80
        ggsave(paste0('cna-blist-chr',chr,'.png'),width=17,height=plot_height,limitsize=FALSE)

}

read.all.germline.vcf <- function(loc=NULL) 
{
    data <- c()
    if (is.null(loc)) loc <- '~/data/pairs.grm.txt'
    dat <- read.table(loc,header=T,sep=" ",stringsAsFactors=F)
    i <- 1
    for (i in 1:nrow(dat)) {
        print(dat[i,]$germline)
        grmloc <- paste0(dat[i,]$germline,"/Vcf/",dat[i,]$germline,"_gcv_HAP.bamstats.vcf")
        somloc <- paste0(dat[i,]$somatic,"/Vcf/","*_combined.bedFiltered.avgDpAd.bamstats.vcf")
        system(paste0(" cat ",somloc, "| grep -v ^## > ", "tmp/",dat[i,]$somatic,".vcf"))
        system(paste0(" cat ",grmloc, "| grep -v ^## > ", "tmp/",dat[i,]$germline,".vcf"))
        som <- read_table2(paste0("tmp/",dat[i,]$somatic,".vcf"))
        som <- as.data.frame(som)
        grm <- read_table2(paste0("tmp/",dat[i,]$germline,".vcf"))
        grm <- as.data.frame(grm)
        names(som) <- gsub("#","",names(grm))
        names(grm) <- gsub("#","",names(grm))
        som$ID <- paste(som$CHROM,som$POS,som$REF,som$ALT,sep="-")
        grm$ID <- paste(grm$CHROM,grm$POS,grm$REF,grm$ALT,sep="-")
        grm <- grm[!grm$ID %in% som$ID,]
        rm(som)
        grm$Sample <- dat[i,]$germline
        grm <- grm[,c(1:7,10,11)]
        vcfdet <- as.data.frame(str_split_fixed(grm[,8],":",n=9),stringsAsFactors=F)
        vcfdet <- as.data.frame((str_split_fixed(vcfdet[,2],",",2)))
        names(vcfdet) <- c("tumour_rd","tumour_ad")
        grm$tumour_rd <- vcfdet$tumour_rd
        grm$tumour_ad <- vcfdet$tumour_ad
        grm <- grm[c(1:7,9,10,11)]
        data <- rbind(data, grm)
    }
    data$tumour_rd <- as.numeric(as.character(data$tumour_rd))
    data$tumour_ad <- as.numeric(as.character(data$tumour_ad))
    data$freq <- data$tumour_ad / (data$tumour_rd + data$tumour_ad)
    write.table(data,'~/data/CPCG_germline_all.vcf',quote=F,row.names=F,col.names=T,sep="\t")
    data <- filter(data, tumour_rd + tumour_ad > 20)
    data <- data[nchar(data$REF)==1,] #no indel
    data <- data[nchar(data$ALT)==1,] #no indel
    write.table(data,'~/data/CPCG_germline_snp.vcf',quote=F,row.names=F,col.names=T,sep="\t")

}

loh.type <- function(snps,rsid) {

    if (mean(filter(snps,ID==rsid)$tumour_rd)==0) return("TOTAL LOH")
    if (mean(filter(snps,type=='RP',ID==rsid)$tumour_rd)==0) return("RP LOH")
    if (mean(filter(snps,type=='MET',ID==rsid)$tumour_rd)==0) return("MET LOH")
    if (mean(filter(snps,type=='MET',gen==5,ID==rsid)$tumour_rd)==0) return("MET G5 LOH")
    if (mean(filter(snps,type=='RP',name=='FR07888602',ID==rsid)$tumour_rd)==0) return("RP G5-I LOH")
    if (mean(filter(snps,type=='RP',name=='FR07888610',ID==rsid)$tumour_rd)==0) return("RP G5-H LOH")
    return("nearbye")
}

make.bed.file <- function(grm,snp) {
    data <- grm
    bp <- filter(data, ID==snp)[1,]
    data <- filter(data, POS >= bp$POS-200, POS <= bp$POS+200, CHROM==bp$CHROM)
    data$chromStart <- data$POS - 1
    data$chromEnd <- data$POS 
    out <- data[c("CHROM","chromStart","chromEnd","ID")]
    return(out)
}

make.bed.file.pos <- function(grm,bp,chr) {
    data <- grm
    data <- filter(data, POS >= bp-200, POS <= bp+200, CHROM==chr)
    if (nrow(data) > 0) {
        data$chromStart <- data$POS - 1
        data$chromEnd <- data$POS 
        out <- data[c("CHROM","chromStart","chromEnd","ID")]
        out$select <- 'FALSE'
        out$type <- 'nearbye'
        return(out)
    }
    return(c())
}



#quit()
#grm_range <- GRanges(data$CHROM,IRanges(data$POS-1,data$POS),mcols=data[,c("REF","ALT","QUAL","tumour_rd","tumour_ad")])
#
##read.germline()
##filter out germline snps that are in somatic already
#
#gaps <- read.gaps()
#glist <- read.genelist()
#facet <- read.facets()
#snps <- read.somatic.snps()
#
#
####testing
#loh <- facet.to.range(filter(facet,tcn.em==1,lcn.em==0))
#gain <- facet.to.range(filter(facet,tcn.em>2))
#
#snps <- snps[snps$mcols.Consequence_Rank <4,]
#for (s in unique(as.list(snps$mcols.sample_name))) {
#    tmp <- snps[snps$mcols.sample_name==s,]
#    #tmploh <- loh[loh$mcols.ID==s,]
#    tmpgain <- gain[gain$mcols.ID==s,]
#    #print(tmp[tmp %over% tmploh,])
#    print(tmp[tmp %over% tmpgain,])
#}
#
#out <- as.data.frame(tmp %>% group_by(sample_name,chromStart,chromEnd,ID,select) %>% summarise(type = paste(type, collapse=", ")))
#
#loh <- facet.to.range(filter(facet,tcn.em==1,lcn.em==0))
#loh <- loh[loh$mcols.ID=="CPCG0103-F1--CPCG0103-B1-F1"]
#tmp <- snps[snps$mcols.sample_name=="CPCG0103-F1--CPCG0103-B1-F1",]
#tmp <- tmp[tmp$mcols.Consequence_Rank<4,]
#tmp[tmp %over% loh,]
#
#for (i in unique(facet$ID)) {
#  tmp <- filter(facet,ID==i)
#  tmp <- tmp[complete.cases(tmp),]
#  tmp <- filter(tmp,tcn.em>2 | tcn.em==2 & lcn.em==0 | tcn.em==1)
#  print(paste0("PGA for ",i))
#  print(sum(tmp$end - tmp$start) / 3095677412)
#}
#
#for (i in unique(facet$ID)) {
#  tmp <- filter(facet,ID==i)
#  tmp <- tmp[complete.cases(tmp),]
#  #tmp <- filter(tmp, lcn.em==0)
#  tmp <- filter(tmp, tcn.em >2)
#  print(paste0("PGA for ",i," ",sum(as.numeric(tmp$end - tmp$start)) / 3095677412))
#  #print(sum(as.numeric(tmp$end - tmp$start)))
#  #print(sum(as.numeric(tmp$end - tmp$start)) / 3095677412)
#}
#
#
#
#for (i in unique(facet$ID)) {
#  tmp <- filter(facet,ID==i)
#  tmp <- tmp[complete.cases(tmp),]
#  #tmp <- filter(tmp, lcn.em==0)
#  tmp <- filter(tmp, tcn.em >2)
#  print(paste0("PGA for ",i," ",sum(as.numeric(tmp$end - tmp$start)) / 3095677412))
#  #print(sum(as.numeric(tmp$end - tmp$start)))
#  #print(sum(as.numeric(tmp$end - tmp$start)) / 3095677412)
#}
#
#
#
#
##check out gene: TPBG


#take list of msigbdb names and return gene set object for each
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
