params <-
list(isSlides = "no")

## ----include=FALSE------------------------------------------------------------
suppressPackageStartupMessages(require(knitr))
knitr::opts_chunk$set(echo = TRUE, tidy = T)


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides != "yes"){
  cat("# RNAseq (part 4)

---
"    
  )
  
}



## ----setup, include=FALSE-----------------------------------------------------
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicAlignments)
library(DESeq2)
library(tximport)
library(org.Mm.eg.db)
library(goseq)
library(DEXSeq)
library(limma)



## ----eval=FALSE,echo=FALSE----------------------------------------------------
## 
## library(BiocParallel)
## library(BatchJobs)
## library(ngsPipeR)
## 
## loadConfig()
## ## register SLURM cluster instructions from the template file
## 
## #setConfig(debug = TRUE, fs.timeout=20)
## 
## #funs <- makeClusterFunctionsSLURM("simple.tmpl")
## funs <- makeClusterFunctionsSLURM("uberSimple.tmpl")
## 
## param <- BatchJobsParam(resources=list(ncpus=1,ntasks=1,walltime=2000),
##                         cluster.functions=funs)
## param$log <- T
## register(param)
## 
## library(GenomicFeatures)
## nonOverlappingExons <- disjointExons(TxDb.Hsapiens.UCSC.hg19.knownGene)
## names(nonOverlappingExons) <- paste(mcols(nonOverlappingExons)$gene_id,
##                                     mcols(nonOverlappingExons)$exonic_part,
##                                     sep="_")
## 
## 
## # myBams <- c("/rugpfs/fs0/ruit/scratch/tcarroll/treg/rnaseq/BAMs/Sorted_T_reg_act_1.bam","/rugpfs/fs0/ruit/scratch/tcarroll/treg/rnaseq/BAMs/Sorted_T_reg_1.bam","/rugpfs/fs0/ruit/scratch/tcarroll/treg/rnaseq/BAMs/Sorted_T_reg_act_2.bam","/rugpfs/fs0/ruit/scratch/tcarroll/treg/rnaseq/BAMs/Sorted_T_reg_2.bam","/rugpfs/fs0/ruit/scratch/tcarroll/treg/rnaseq/BAMs/Sorted_T_reg_act_3.bam")
## myBams <- dir("/rugpfs/fs0/brc/scratch/tcarroll/Test/ptbp1/subreadAlign_salmonDE_txdbAnno/BAMs/",pattern="Sorted.*.bam$",full.names = TRUE)
## senescence_ExonCounts <- summarizeOverlaps(nonOverlappingExons,
##                                     myBams,
##                                     ignore.strand = TRUE,
##                                     inter.feature=FALSE)
## 
## senescence_ExonCounts <- senescence_ExonCounts[,grep("shPTBP1_53|vec_4OHT",colnames(senescence_ExonCounts))]
## 
## colData(senescence_ExonCounts)$condition <- c("shPTBP1_53", "shPTBP1_53", "shPTBP1_53","senescence","senescence","senescence")
## rownames(senescence_ExonCounts) <- NULL
## 
## table(rowSums(assay(senescence_ExonCounts)) > 10)
## senescence_ExonCounts <- senescence_ExonCounts[rowSums(assay(senescence_ExonCounts)) > 10,]
## table(apply(assay(senescence_ExonCounts),1,function(x)sum(x > 10)) > 2)
## senescence_ExonCounts <- senescence_ExonCounts[apply(assay(senescence_ExonCounts),1,function(x)sum(x > 10)) > 2]
## 
## 
## 
## dxd <- DEXSeqDataSetFromSE(senescence_ExonCounts,design= ~ sample + exon + condition:exon)
## dxd <- estimateSizeFactors(dxd)
## dxd <- estimateDispersions(dxd)
## dxd <- testForDEU(dxd, reducedModel = ~ sample + exon)
## 
## dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")
## 
## dxr1 = DEXSeqResults( dxd )
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## myBams <- c("/rugpfs/fs0/ruit/scratch/tcarroll/Tissues/rnaseq/BAMs/Sorted_Heart_1.bam","/rugpfs/fs0/ruit/scratch/tcarroll/Tissues/rnaseq/BAMs/Sorted_Heart_2.bam","/rugpfs/fs0/ruit/scratch/tcarroll/Tissues/rnaseq/BAMs/Sorted_Liver_1.bam","/rugpfs/fs0/ruit/scratch/tcarroll/Tissues/rnaseq/BAMs/Sorted_Liver_2.bam","/rugpfs/fs0/ruit/scratch/tcarroll/Tissues/rnaseq/BAMs/Sorted_Kidney_1.bam","/rugpfs/fs0/ruit/scratch/tcarroll/Tissues/rnaseq/BAMs/Sorted_Kidney_2.bam")
## 
## 
## tissueExonCounts <- summarizeOverlaps(nonOverlappingExons,
##                                     myBams,
##                                     ignore.strand = TRUE,
##                                     inter.feature=FALSE)
## 
## tissueExonCounts <- tissueExonCounts
## colData(tissueExonCounts)$tissue <- c("Heart", "Heart", "Liver","Liver")
## rownames(tissueExonCounts) <- NULL
## 
## table(rowSums(assay(tissueExonCounts)) > 10)
## tissueExonCounts <- tissueExonCounts[rowSums(assay(tissueExonCounts)) > 10,]
## table(apply(assay(tissueExonCounts),1,function(x)sum(x > 10)) > 2)
## tissueExonCounts <- tissueExonCounts[apply(assay(tissueExonCounts),1,function(x)sum(x > 10)) > 2]
## 
## 
## dxd <- DEXSeqDataSetFromSE(tissueExonCounts,design= ~ sample + exon + tissue:exon)
## dxd <- estimateSizeFactors(dxd)
## dxd <- estimateDispersions(dxd)
## dxd <- testForDEU(dxd, reducedModel = ~ sample + exon)
## 
## dxd = estimateExonFoldChanges( dxd, fitExpToVar="tissue")
## 
## dxr1 = DEXSeqResults(dxd)
## 
## 
## 
## 
## 
## load("../../../exonCounts.RData")
## library(limma)
## library(edgeR)
## tregExonCounts <- tregExonCounts[rowSums(assay(tregExonCounts)) > 10,]
## 
## cefcef <- names(table(as.vector(unlist(rowData(tregExonCounts)$gene_id)))[table(as.vector(unlist(rowData(tregExonCounts)$gene_id))) > 1])
## 
## #tregExonCounts <- tregExonCounts[as.vector(unlist(rowData(tregExonCounts)$gene_id)) %in% cefcef,]
## dge <- DGEList(counts=assay(tregExonCounts))
## dge$genes <- data.frame(GeneID=as.vector(unlist(rowData(tregExonCounts)$gene_id)),ExonID=paste0(as.vector(unlist(rowData(tregExonCounts)$gene_id)),"_",as.vector(unlist(rowData(tregExonCounts)$exonic_part))))
## dge <- calcNormFactors(dge)
## 
## 
## f <- factor(c("Act", "Naive", "Act","Naive","Act"))
## design <- model.matrix(~0+f)
## colnames(design) <- c("Act","Naive")
## 
## v <- voom(dge, design, plot=TRUE)
## fit <- lmFit(v, design)
## ex <- diffSplice(fit, geneid="GeneID",exonid="ExonID")
## fdefe <- topSplice(ex,test = "F",number = 50)
## 
## function (fit, geneid, exonid = NULL, robust = FALSE, verbose = TRUE)
## {
##   fit <- fit
##   geneid <- "GeneID"
##    exonid = "ExonID"
##    robust = FALSE
##    verbose = TRUE
##     exon.genes <- fit$genes
##     if (is.null(exon.genes))
##         exon.genes <- data.frame(ExonID = 1:nrow(fit))
##     if (length(geneid) == 1) {
##         genecolname <- as.character(geneid)
##         geneid <- exon.genes[[genecolname]]
##     }else {
##         exon.genes$GeneID <- geneid
##         genecolname <- "GeneID"
##     }
##     if (is.null(exonid)) {
##         exoncolname <- NULL
##     }else {
##         if (length(exonid) == 1) {
##             exoncolname <- as.character(exonid)
##             exonid <- exon.genes[[exoncolname]]
##         }else {
##             exon.genes$ExonID <- exonid
##             exoncolname <- "ExonID"
##         }
##     }
##     if (anyNA(geneid)) {
##         isna <- which(is.na(geneid))
##         geneid[isna] <- paste0("NA", 1:length(isna))
##     }
##     if (is.null(exonid))
##         o <- order(geneid)
##     else o <- order(geneid, exonid)
##     geneid <- geneid[o]
##     exon.genes <- exon.genes[o, , drop = FALSE]
##     exon.coefficients <- fit$coefficients[o, , drop = FALSE]
##     exon.stdev.unscaled <- fit$stdev.unscaled[o, , drop = FALSE]
##     exon.df.residual <- fit$df.residual[o]
##     exon.s2 <- fit$sigma[o]^2
##     exon.stat <- cbind(1, exon.df.residual, exon.s2)
##     gene.sum <- rowsum(exon.stat, geneid, reorder = FALSE)
##     gene.nexons <- gene.sum[, 1]
##     gene.df.residual <- gene.sum[, 2]
##     gene.s2 <- gene.sum[, 3]/gene.sum[, 1]
##     if (verbose) {
##         cat("Total number of exons: ", length(geneid), "\\n")
##         cat("Total number of genes: ", length(gene.nexons), "\\n")
##         cat("Number of genes with 1 exon: ", sum(gene.nexons ==
##             1), "\\n")
##         cat("Mean number of exons in a gene: ", round(mean(gene.nexons),
##             0), "\\n")
##         cat("Max number of exons in a gene: ", max(gene.nexons),
##             "\\n")
##     }
##     squeeze <- squeezeVar(var = gene.s2, df = gene.df.residual,
##         robust = robust)
##     gene.keep <- gene.nexons > 1
##     ngenes <- sum(gene.keep)
##     if (ngenes == 0)
##         stop("No genes with more than one exon")
##     exon.keep <- rep(gene.keep, gene.nexons)
##     geneid <- geneid[exon.keep]
##     exon.genes <- exon.genes[exon.keep, , drop = FALSE]
##     exon.coefficients <- exon.coefficients[exon.keep, , drop = FALSE]
##     exon.stdev.unscaled <- exon.stdev.unscaled[exon.keep, , drop = FALSE]
##     exon.df.residual <- exon.df.residual[exon.keep]
##     gene.nexons <- gene.nexons[gene.keep]
##     gene.df.test <- gene.nexons - 1
##     gene.df.residual <- gene.df.residual[gene.keep]
##     if (robust)
##         squeeze$df.prior <- squeeze$df.prior[gene.keep]
##     gene.df.total <- gene.df.residual + squeeze$df.prior
##     gene.df.total <- pmin(gene.df.total, sum(gene.df.residual))
##     gene.s2.post <- squeeze$var.post[gene.keep]
##     u2 <- 1/exon.stdev.unscaled^2
##     u2.rowsum <- rowsum(u2, geneid, reorder = FALSE)
##     gene.betabar <- rowsum(exon.coefficients * u2, geneid, reorder = FALSE)/u2.rowsum
##     g <- rep(1:ngenes, times = gene.nexons)
##     exon.coefficients <- exon.coefficients - gene.betabar[g,
##         , drop = FALSE]
##     exon.t <- exon.coefficients/exon.stdev.unscaled/sqrt(gene.s2.post[g])
##     gene.F <- rowsum(exon.t^2, geneid, reorder = FALSE)/gene.df.test
##     exon.1mleverage <- 1 - (u2/u2.rowsum[g, , drop = FALSE])
##     exon.coefficients <- exon.coefficients/exon.1mleverage
##     exon.t <- exon.t/sqrt(exon.1mleverage)
##     exon.p.value <- 2 * pt(abs(exon.t), df = gene.df.total[g],
##         lower.tail = FALSE)
##     gene.F.p.value <- pf(gene.F, df1 = gene.df.test, df2 = gene.df.total,
##         lower.tail = FALSE)
##     out <- new("MArrayLM", list())
##     out$genes <- exon.genes
##     out$genecolname <- genecolname
##     out$exoncolname <- exoncolname
##     out$coefficients <- exon.coefficients
##     out$t <- exon.t
##     out$p.value <- exon.p.value
##     out$gene.df.prior <- squeeze$df.prior
##     out$gene.df.residual <- gene.df.residual
##     out$gene.df.total <- gene.df.total
##     out$gene.s2 <- gene.s2[gene.keep]
##     out$gene.s2.post <- gene.s2.post
##     out$gene.F <- gene.F
##     out$gene.F.p.value <- gene.F.p.value
##     gene.lastexon <- cumsum(gene.nexons)
##     gene.firstexon <- gene.lastexon - gene.nexons + 1
##     no <- logical(nrow(exon.genes))
##     isdup <- vapply(exon.genes, duplicated, no)[-gene.firstexon,
##         , drop = FALSE]
##     isgenelevel <- apply(isdup, 2, all)
##     out$gene.genes <- exon.genes[gene.lastexon, isgenelevel,
##         drop = FALSE]
##     out$gene.genes$NExons <- gene.nexons
##     out$gene.firstexon <- gene.firstexon
##     out$gene.lastexon <- gene.lastexon
##     penalty <- rep_len(1L, length(g))
##     penalty[gene.lastexon] <- 1L - gene.nexons
##     penalty <- cumsum(penalty)[-gene.lastexon]
##     penalty <- penalty/rep(gene.nexons - 1L, gene.nexons - 1L)
##     g2 <- g[-gene.lastexon]
##     out$gene.simes.p.value <- gene.F.p.value
##     for (j in 1:ncol(fit)) {
##         o <- order(g, exon.p.value[, j])
##         p.adj <- pmin(exon.p.value[o, j][-gene.lastexon]/penalty,
##             1)
##         o <- order(g2, p.adj)
##         out$gene.simes.p.value[, j] <- p.adj[o][gene.firstexon -
##             0L:(ngenes - 1L)]
##     }
##     out
## }


## ----eval=FALSE,echo=TRUE,tidy=FALSE------------------------------------------
## library(Rsamtools)
## bamFilesToCount <- c("Sorted_shPTBP1_53_rep1.bam","Sorted_shPTBP1_53_rep2.bam","Sorted_shPTBP1_53_rep3.bam",
##                      "Sorted_vec_4OHT_rep1.bam","Sorted_vec_4OHT_rep2.bam","Sorted_vec_4OHT_rep3.bam" )
## myBams <- BamFileList(bamFilesToCount,yieldSize = 10000)
## 


## ----eval=TRUE,echo=TRUE------------------------------------------------------
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
nonOverlappingExons <- disjointExons(TxDb.Hsapiens.UCSC.hg19.knownGene)
names(nonOverlappingExons) <- paste(mcols(nonOverlappingExons)$gene_id,
                                    mcols(nonOverlappingExons)$exonic_part,
                                    sep="_")
nonOverlappingExons[1:2,]


## ----eval=FALSE,echo=TRUE-----------------------------------------------------
## library(GenomicAlignments)
## senescence_ExonCounts <- summarizeOverlaps(nonOverlappingExons,
##                                     myBam,
##                                     ignore.strand = TRUE,
##                                     inter.feature=FALSE)
## senescence_ExonCounts


## ----gC1z,eval=TRUE,echo=FALSE,cache=FALSE------------------------------------
load("data/senescence_ExonCounts.RData")
senescence_ExonCounts


## ----gC2,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gC1z"----------------------
assay(senescence_ExonCounts)[1:2,]


## ----gC3,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gC1z"----------------------
rowRanges(senescence_ExonCounts)[1:5,]


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Differential transcript usage

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Differential transcript usage

---
"    
  )
  
}


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# DEXSeq

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# DEXSeq

---
"    
  )
  
}


## ----de1z,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gC1z"---------------------
metaData <- data.frame(condition=c("shPTBP1_53", "shPTBP1_53", "shPTBP1_53","senescence","senescence","senescence"),
                       row.names = colnames(senescence_ExonCounts))
metaData


## ----de22z,eval=TRUE,echo=TRUE,cache=TRUE,dependson="de1z"--------------------
countMatrix <- assay(senescence_ExonCounts)
countGRanges <- rowRanges(senescence_ExonCounts)
geneIDs <- as.vector(unlist(countGRanges$gene_id))
exonIDs <- rownames(countMatrix)


## ----include=FALSE------------------------------------------------------------
rownames(senescence_ExonCounts) <- NULL


## ----de222z,eval=TRUE,echo=TRUE,cache=TRUE,dependson="de22z",warning=FALSE,message=FALSE----
library(DEXSeq)
ddx <- DEXSeqDataSet(countMatrix,metaData,
                     design= ~ sample + exon + condition:exon,
                     featureID = exonIDs,
                     groupID = geneIDs)
ddx


## ----de222zTT,include=FALSE,cache=TRUE,dependson="de222z"---------------------
mcols(senescence_ExonCounts)$gene_id <- unlist(mcols(senescence_ExonCounts)$gene_id)


## ----de2223z,eval=TRUE,echo=TRUE,cache=TRUE,dependson="de222zTT",warning=FALSE,message=FALSE----
colData(senescence_ExonCounts)$condition <- factor(c("shPTBP1_53", "shPTBP1_53", "shPTBP1_53",
                                       "senescence","senescence","senescence"))
ddx <- DEXSeqDataSetFromSE(senescence_ExonCounts,
                     design= ~ sample + exon + condition:exon)
ddx


## ----de223l,eval=TRUE,echo=TRUE,cache=TRUE,dependson="de2223z",warning=FALSE,message=FALSE----
ToFilter <- apply(assay(senescence_ExonCounts),1,function(x) sum(x > 10)) >= 2
senescence_ExonCounts <- senescence_ExonCounts[ToFilter,]
table(ToFilter)
ddx <- DEXSeqDataSetFromSE(senescence_ExonCounts,
                     design= ~ sample + exon + condition:exon)
ddx


## ----de22k4,eval=FALSE--------------------------------------------------------
## ddx <- estimateSizeFactors(ddx)
## ddx <- estimateDispersions(ddx)
## ddx <- testForDEU(ddx, reducedModel = ~ sample + exon)
## ddx = estimateExonFoldChanges( ddx, fitExpToVar="condition")


## ----include=FALSE------------------------------------------------------------
load("data/senescence_Exon_Res.RData")
load("data/senescence_Exon_DEXSeqDataset.RData")
ddx<-dxd



## ----da1s,eval=FALSE,echo=TRUE,cache=TRUE-------------------------------------
## dxr1 <- DEXSeqResults(ddx)
## dxr1 <- dxr1[order(dxr1$pvalue),]


## ----da1aa,eval=TRUE,echo=FALSE,cache=TRUE------------------------------------
dxr1 <- dxr1[order(dxr1$pvalue),]


## ----de224,eval=FALSE---------------------------------------------------------
## dxr1 <- DEXSeq(ddx)
## dxr1 <- dxr1[order(dxr1$pvalue),]


## ----da2,eval=TRUE,echo=TRUE,cache=TRUE,dependson="da2"-----------------------
dxr1[1,]


## ----da3,eval=TRUE,echo=TRUE,cache=TRUE,dependson="da2"-----------------------
as.data.frame(dxr1)[1,]


## ----da4,eval=TRUE,echo=TRUE,cache=TRUE,dependson="da2",fig.height=4,fig.width=7----
plotDEXSeq(dxr1,"57142")


## ----da5c,eval=TRUE,echo=TRUE,cache=TRUE,dependson="da2",fig.height=5,fig.width=7----
plotDEXSeq(dxr1,"57142",displayTranscripts=TRUE)


## ----da5,eval=TRUE,echo=TRUE,cache=TRUE,dependson="da2",fig.height=9,fig.width=12----
plotDEXSeq(dxr1,"7168",displayTranscripts=TRUE)


## ----da6,eval=FALSE-----------------------------------------------------------
## DEXSeqHTML(dxr1,"57142")


## ----da7,eval=TRUE,echo=TRUE,cache=TRUE,dependson="da2",fig.height=5,fig.width=7----
geneQ <- perGeneQValue(dxr1)
geneQ[1:10]


## ----da8,eval=TRUE,echo=TRUE,cache=TRUE,dependson="da7",fig.height=5,fig.width=7----
gQFrame <- as.data.frame(geneQ,row.names = names(geneQ))
dxDF <- as.data.frame(dxr1)
dxDF <-merge(gQFrame,dxDF,by.x=0,by.y=1,all=TRUE)
dxDF <- dxDF[order(dxDF$geneQ,dxDF$pvalue),]
dxDF[1:2,]


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# RNAseq with multiple groups

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# RNAseq with multiple groups

---
"    
  )
  
}


## ----tissue,eval=TRUE,echo=TRUE,cache=TRUE,fig.height=5,fig.width=7-----------
load("data/gC_TissueFull.RData")
geneCounts_Tissue


## ----tissue2,eval=TRUE,echo=TRUE,cache=TRUE,fig.height=5,fig.width=7,dependson="tissue"----
colData(geneCounts_Tissue)$Tissue <- c("Heart","Heart",
                                       "Kidney","Kidney",
                                       "Liver","Liver")


## ----tissue3,eval=TRUE,echo=TRUE,cache=TRUE,fig.height=5,fig.width=7,dependson="tissue2",warning=FALSE,message=FALSE----
dds <- DESeqDataSet(geneCounts_Tissue,design = ~Tissue)
dds <- DESeq(dds)


## ----tissue4,eval=TRUE,echo=TRUE,cache=TRUE,fig.height=5,fig.width=7,dependson="tissue3"----
heartVsLiver <- results(dds,c("Tissue","Heart","Liver"))
heartVskidney <- results(dds,c("Tissue","Heart","Kidney"))
heartVskidney


## ----tissue5,eval=TRUE,echo=TRUE,cache=TRUE,fig.height=5,fig.width=7,dependson="tissue4"----
heartVsLiverDF <- as.data.frame(heartVsLiver)
heartVskidneyDF <- as.data.frame(heartVskidney)

heartVsLiverDF <- heartVsLiverDF[!is.na(heartVsLiverDF$padj),]
heartVskidneyDF <- heartVskidneyDF[!is.na(heartVskidneyDF$padj),]

heartVsLiverDF <- heartVsLiverDF[,c("log2FoldChange","padj")]
heartVskidneyDF <- heartVskidneyDF[,c("log2FoldChange","padj")]


## ----tissue6,eval=TRUE,echo=TRUE,cache=TRUE,fig.height=5,fig.width=7,dependson="tissue5"----
colnames(heartVskidneyDF) <- paste0("HeartVsKidney","_",colnames(heartVskidneyDF))
colnames(heartVsLiverDF) <- paste0("HeartVsLiver","_",colnames(heartVsLiverDF))
fullTable <- merge(heartVsLiverDF,heartVskidneyDF,by=0)
fullTable[1:2,]


## ----tissue7,eval=TRUE,echo=TRUE,cache=TRUE,fig.height=5,fig.width=7,dependson="tissue8"----
upInHeart <- fullTable$HeartVsLiver_log2FoldChange > 0 &
             fullTable$HeartVsKidney_log2FoldChange > 0 &
             fullTable$HeartVsLiver_padj < 0.05 &
             fullTable$HeartVsKidney_padj < 0.05
upInHeartTable <- fullTable[upInHeart,]
upInHeartTable[1:2,]


## ----tissue9,eval=TRUE,echo=TRUE,cache=TRUE,fig.height=5,fig.width=7,dependson="tissue8"----
forVenn <- data.frame(UpvsLiver=fullTable$HeartVsLiver_log2FoldChange > 0 &
                       fullTable$HeartVsLiver_padj < 0.05,
                     UpvsKidney=fullTable$HeartVsKidney_log2FoldChange > 0 &
                       fullTable$HeartVsKidney_padj < 0.05)
forVenn[1:3,]


## ----tissue10,eval=TRUE,echo=TRUE,cache=TRUE,dependson="tissue9"--------------
library(limma)
vennDiagram(forVenn)


## ----tissue11,eval=TRUE,echo=TRUE,cache=TRUE,dependson="tissue10",message=FALSE,warning=FALSE----
dds2 <- DESeq(dds,test="LRT",reduced=~1)
acrossGroups <- results(dds2)
acrossGroups <- acrossGroups[order(acrossGroups$pvalue),]
acrossGroups[1:3,]


## ----tissue12,eval=TRUE,echo=TRUE,cache=TRUE,dependson="tissue11",fig.width=6,fig.height=5----
plotCounts(dds2,gene="17888",intgroup = "Tissue")


## ----tissue13,eval=TRUE,echo=TRUE,cache=TRUE,dependson="tissue12",fig.width=6,fig.height=5----
acrossGroups <- as.data.frame(acrossGroups[,c("pvalue","padj")])
colnames(acrossGroups) <- paste0("AllGroups","_",colnames(acrossGroups))
allGroups <- merge(fullTable,acrossGroups,by.x=1,by.y=0,all=TRUE)
allGroups <- allGroups[order(allGroups$AllGroups_pvalue),]
allGroups[1:2,]

