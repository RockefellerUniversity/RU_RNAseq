params <-
list(isSlides = "no")

## ----include=FALSE------------------------------------------------------------
suppressPackageStartupMessages(require(knitr))
knitr::opts_chunk$set(echo = TRUE, tidy = T)


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides != "yes"){
  cat("# RNAseq (part 3)

---
"    
  )
  
}



## ----setup, include=FALSE-----------------------------------------------------
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicAlignments)
library(DESeq2)
library(org.Mm.eg.db)
library(goseq)
library(limma)



## ----eval=FALSE,echo=FALSE----------------------------------------------------
# 
# library(BiocParallel)
# library(BatchJobs)
# library(ngsPipeR)
# 
# loadConfig()
# ## register SLURM cluster instructions from the template file
# 
# #setConfig(debug = TRUE, fs.timeout=20)
# 
# #funs <- makeClusterFunctionsSLURM("simple.tmpl")
# funs <- makeClusterFunctionsSLURM("uberSimple.tmpl")
# 
# param <- BatchJobsParam(resources=list(ncpus=1,ntasks=1,walltime=2000),
#                         cluster.functions=funs)
# param$log <- T
# register(param)
# 
# library(GenomicFeatures)
# nonOverlappingExons <- disjointExons(TxDb.Hsapiens.UCSC.hg19.knownGene)
# names(nonOverlappingExons) <- paste(mcols(nonOverlappingExons)$gene_id,
#                                     mcols(nonOverlappingExons)$exonic_part,
#                                     sep="_")
# 
# 
# # myBams <- c("/rugpfs/fs0/ruit/scratch/tcarroll/treg/rnaseq/BAMs/Sorted_T_reg_act_1.bam","/rugpfs/fs0/ruit/scratch/tcarroll/treg/rnaseq/BAMs/Sorted_T_reg_1.bam","/rugpfs/fs0/ruit/scratch/tcarroll/treg/rnaseq/BAMs/Sorted_T_reg_act_2.bam","/rugpfs/fs0/ruit/scratch/tcarroll/treg/rnaseq/BAMs/Sorted_T_reg_2.bam","/rugpfs/fs0/ruit/scratch/tcarroll/treg/rnaseq/BAMs/Sorted_T_reg_act_3.bam")
# myBams <- dir("/rugpfs/fs0/brc/scratch/tcarroll/Test/ptbp1/subreadAlign_salmonDE_txdbAnno/BAMs/",pattern="Sorted.*.bam$",full.names = TRUE)
# senescence_ExonCounts <- summarizeOverlaps(nonOverlappingExons,
#                                     myBams,
#                                     ignore.strand = TRUE,
#                                     inter.feature=FALSE)
# 
# senescence_ExonCounts <- senescence_ExonCounts[,grep("shPTBP1_53|vec_4OHT",colnames(senescence_ExonCounts))]
# 
# colData(senescence_ExonCounts)$condition <- c("shPTBP1_53", "shPTBP1_53", "shPTBP1_53","senescence","senescence","senescence")
# rownames(senescence_ExonCounts) <- NULL
# 
# table(rowSums(assay(senescence_ExonCounts)) > 10)
# senescence_ExonCounts <- senescence_ExonCounts[rowSums(assay(senescence_ExonCounts)) > 10,]
# table(apply(assay(senescence_ExonCounts),1,function(x)sum(x > 10)) > 2)
# senescence_ExonCounts <- senescence_ExonCounts[apply(assay(senescence_ExonCounts),1,function(x)sum(x > 10)) > 2]
# 
# 
# 
# dxd <- DEXSeqDataSetFromSE(senescence_ExonCounts,design= ~ sample + exon + condition:exon)
# dxd <- estimateSizeFactors(dxd)
# dxd <- estimateDispersions(dxd)
# dxd <- testForDEU(dxd, reducedModel = ~ sample + exon)
# 
# dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")
# 
# dxr1 = DEXSeqResults( dxd )
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# myBams <- c("/rugpfs/fs0/ruit/scratch/tcarroll/Tissues/rnaseq/BAMs/Sorted_Heart_1.bam","/rugpfs/fs0/ruit/scratch/tcarroll/Tissues/rnaseq/BAMs/Sorted_Heart_2.bam","/rugpfs/fs0/ruit/scratch/tcarroll/Tissues/rnaseq/BAMs/Sorted_Liver_1.bam","/rugpfs/fs0/ruit/scratch/tcarroll/Tissues/rnaseq/BAMs/Sorted_Liver_2.bam","/rugpfs/fs0/ruit/scratch/tcarroll/Tissues/rnaseq/BAMs/Sorted_Kidney_1.bam","/rugpfs/fs0/ruit/scratch/tcarroll/Tissues/rnaseq/BAMs/Sorted_Kidney_2.bam")
# 
# 
# tissueExonCounts <- summarizeOverlaps(nonOverlappingExons,
#                                     myBams,
#                                     ignore.strand = TRUE,
#                                     inter.feature=FALSE)
# 
# tissueExonCounts <- tissueExonCounts
# colData(tissueExonCounts)$tissue <- c("Heart", "Heart", "Liver","Liver")
# rownames(tissueExonCounts) <- NULL
# 
# table(rowSums(assay(tissueExonCounts)) > 10)
# tissueExonCounts <- tissueExonCounts[rowSums(assay(tissueExonCounts)) > 10,]
# table(apply(assay(tissueExonCounts),1,function(x)sum(x > 10)) > 2)
# tissueExonCounts <- tissueExonCounts[apply(assay(tissueExonCounts),1,function(x)sum(x > 10)) > 2]
# 
# 
# dxd <- DEXSeqDataSetFromSE(tissueExonCounts,design= ~ sample + exon + tissue:exon)
# dxd <- estimateSizeFactors(dxd)
# dxd <- estimateDispersions(dxd)
# dxd <- testForDEU(dxd, reducedModel = ~ sample + exon)
# 
# dxd = estimateExonFoldChanges( dxd, fitExpToVar="tissue")
# 
# dxr1 = DEXSeqResults(dxd)
# 
# 
# 
# 
# 
# load("../../../exonCounts.RData")
# library(limma)
# library(edgeR)
# tregExonCounts <- tregExonCounts[rowSums(assay(tregExonCounts)) > 10,]
# 
# cefcef <- names(table(as.vector(unlist(rowData(tregExonCounts)$gene_id)))[table(as.vector(unlist(rowData(tregExonCounts)$gene_id))) > 1])
# 
# #tregExonCounts <- tregExonCounts[as.vector(unlist(rowData(tregExonCounts)$gene_id)) %in% cefcef,]
# dge <- DGEList(counts=assay(tregExonCounts))
# dge$genes <- data.frame(GeneID=as.vector(unlist(rowData(tregExonCounts)$gene_id)),ExonID=paste0(as.vector(unlist(rowData(tregExonCounts)$gene_id)),"_",as.vector(unlist(rowData(tregExonCounts)$exonic_part))))
# dge <- calcNormFactors(dge)
# 
# 
# f <- factor(c("Act", "Naive", "Act","Naive","Act"))
# design <- model.matrix(~0+f)
# colnames(design) <- c("Act","Naive")
# 
# v <- voom(dge, design, plot=TRUE)
# fit <- lmFit(v, design)
# ex <- diffSplice(fit, geneid="GeneID",exonid="ExonID")
# fdefe <- topSplice(ex,test = "F",number = 50)
# 
# function (fit, geneid, exonid = NULL, robust = FALSE, verbose = TRUE)
# {
#   fit <- fit
#   geneid <- "GeneID"
#    exonid = "ExonID"
#    robust = FALSE
#    verbose = TRUE
#     exon.genes <- fit$genes
#     if (is.null(exon.genes))
#         exon.genes <- data.frame(ExonID = 1:nrow(fit))
#     if (length(geneid) == 1) {
#         genecolname <- as.character(geneid)
#         geneid <- exon.genes[[genecolname]]
#     }else {
#         exon.genes$GeneID <- geneid
#         genecolname <- "GeneID"
#     }
#     if (is.null(exonid)) {
#         exoncolname <- NULL
#     }else {
#         if (length(exonid) == 1) {
#             exoncolname <- as.character(exonid)
#             exonid <- exon.genes[[exoncolname]]
#         }else {
#             exon.genes$ExonID <- exonid
#             exoncolname <- "ExonID"
#         }
#     }
#     if (anyNA(geneid)) {
#         isna <- which(is.na(geneid))
#         geneid[isna] <- paste0("NA", 1:length(isna))
#     }
#     if (is.null(exonid))
#         o <- order(geneid)
#     else o <- order(geneid, exonid)
#     geneid <- geneid[o]
#     exon.genes <- exon.genes[o, , drop = FALSE]
#     exon.coefficients <- fit$coefficients[o, , drop = FALSE]
#     exon.stdev.unscaled <- fit$stdev.unscaled[o, , drop = FALSE]
#     exon.df.residual <- fit$df.residual[o]
#     exon.s2 <- fit$sigma[o]^2
#     exon.stat <- cbind(1, exon.df.residual, exon.s2)
#     gene.sum <- rowsum(exon.stat, geneid, reorder = FALSE)
#     gene.nexons <- gene.sum[, 1]
#     gene.df.residual <- gene.sum[, 2]
#     gene.s2 <- gene.sum[, 3]/gene.sum[, 1]
#     if (verbose) {
#         cat("Total number of exons: ", length(geneid), "\\n")
#         cat("Total number of genes: ", length(gene.nexons), "\\n")
#         cat("Number of genes with 1 exon: ", sum(gene.nexons ==
#             1), "\\n")
#         cat("Mean number of exons in a gene: ", round(mean(gene.nexons),
#             0), "\\n")
#         cat("Max number of exons in a gene: ", max(gene.nexons),
#             "\\n")
#     }
#     squeeze <- squeezeVar(var = gene.s2, df = gene.df.residual,
#         robust = robust)
#     gene.keep <- gene.nexons > 1
#     ngenes <- sum(gene.keep)
#     if (ngenes == 0)
#         stop("No genes with more than one exon")
#     exon.keep <- rep(gene.keep, gene.nexons)
#     geneid <- geneid[exon.keep]
#     exon.genes <- exon.genes[exon.keep, , drop = FALSE]
#     exon.coefficients <- exon.coefficients[exon.keep, , drop = FALSE]
#     exon.stdev.unscaled <- exon.stdev.unscaled[exon.keep, , drop = FALSE]
#     exon.df.residual <- exon.df.residual[exon.keep]
#     gene.nexons <- gene.nexons[gene.keep]
#     gene.df.test <- gene.nexons - 1
#     gene.df.residual <- gene.df.residual[gene.keep]
#     if (robust)
#         squeeze$df.prior <- squeeze$df.prior[gene.keep]
#     gene.df.total <- gene.df.residual + squeeze$df.prior
#     gene.df.total <- pmin(gene.df.total, sum(gene.df.residual))
#     gene.s2.post <- squeeze$var.post[gene.keep]
#     u2 <- 1/exon.stdev.unscaled^2
#     u2.rowsum <- rowsum(u2, geneid, reorder = FALSE)
#     gene.betabar <- rowsum(exon.coefficients * u2, geneid, reorder = FALSE)/u2.rowsum
#     g <- rep(1:ngenes, times = gene.nexons)
#     exon.coefficients <- exon.coefficients - gene.betabar[g,
#         , drop = FALSE]
#     exon.t <- exon.coefficients/exon.stdev.unscaled/sqrt(gene.s2.post[g])
#     gene.F <- rowsum(exon.t^2, geneid, reorder = FALSE)/gene.df.test
#     exon.1mleverage <- 1 - (u2/u2.rowsum[g, , drop = FALSE])
#     exon.coefficients <- exon.coefficients/exon.1mleverage
#     exon.t <- exon.t/sqrt(exon.1mleverage)
#     exon.p.value <- 2 * pt(abs(exon.t), df = gene.df.total[g],
#         lower.tail = FALSE)
#     gene.F.p.value <- pf(gene.F, df1 = gene.df.test, df2 = gene.df.total,
#         lower.tail = FALSE)
#     out <- new("MArrayLM", list())
#     out$genes <- exon.genes
#     out$genecolname <- genecolname
#     out$exoncolname <- exoncolname
#     out$coefficients <- exon.coefficients
#     out$t <- exon.t
#     out$p.value <- exon.p.value
#     out$gene.df.prior <- squeeze$df.prior
#     out$gene.df.residual <- gene.df.residual
#     out$gene.df.total <- gene.df.total
#     out$gene.s2 <- gene.s2[gene.keep]
#     out$gene.s2.post <- gene.s2.post
#     out$gene.F <- gene.F
#     out$gene.F.p.value <- gene.F.p.value
#     gene.lastexon <- cumsum(gene.nexons)
#     gene.firstexon <- gene.lastexon - gene.nexons + 1
#     no <- logical(nrow(exon.genes))
#     isdup <- vapply(exon.genes, duplicated, no)[-gene.firstexon,
#         , drop = FALSE]
#     isgenelevel <- apply(isdup, 2, all)
#     out$gene.genes <- exon.genes[gene.lastexon, isgenelevel,
#         drop = FALSE]
#     out$gene.genes$NExons <- gene.nexons
#     out$gene.firstexon <- gene.firstexon
#     out$gene.lastexon <- gene.lastexon
#     penalty <- rep_len(1L, length(g))
#     penalty[gene.lastexon] <- 1L - gene.nexons
#     penalty <- cumsum(penalty)[-gene.lastexon]
#     penalty <- penalty/rep(gene.nexons - 1L, gene.nexons - 1L)
#     g2 <- g[-gene.lastexon]
#     out$gene.simes.p.value <- gene.F.p.value
#     for (j in 1:ncol(fit)) {
#         o <- order(g, exon.p.value[, j])
#         p.adj <- pmin(exon.p.value[o, j][-gene.lastexon]/penalty,
#             1)
#         o <- order(g2, p.adj)
#         out$gene.simes.p.value[, j] <- p.adj[o][gene.firstexon -
#             0L:(ngenes - 1L)]
#     }
#     out
# }


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
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



## ----gC1,eval=F,echo=FALSE,warning=FALSE,message=FALSE------------------------
# load("data/gC_TissueFull.RData")
# geneCounts <- geneCounts_Tissue
# colData(geneCounts)$Tissue <- c("Heart", "Heart","Kidney","Kidney","Liver","Liver")
# geneCounts <- geneCounts[rowSums(assay(geneCounts)) > quantile(rowSums(assay(geneCounts)),0.4)]
# 
# save(geneCounts, file = "data/gC_TissueFull.RData")


## -----------------------------------------------------------------------------
load("data/gC_TissueFull.RData")
geneCounts



## ----eval=F-------------------------------------------------------------------
# setClassUnion("ExpData", c("matrix", "SummarizedExperiment"))


## ----tissue3,eval=TRUE,echo=TRUE,fig.height=5,fig.width=7,warning=FALSE,message=FALSE----
dds <- DESeqDataSet(geneCounts, design = ~Tissue)
dds <- DESeq(dds)


## ----tissue4,eval=TRUE,echo=TRUE,fig.height=5,fig.width=7---------------------
LiverVskidney <- results(dds,c("Tissue","Liver","Kidney"))
heartVsLiver <- results(dds,c("Tissue","Heart","Liver"))
heartVskidney <- results(dds,c("Tissue","Heart","Kidney"))
heartVskidney


## ----tissue5,eval=TRUE,echo=TRUE,fig.height=5,fig.width=7---------------------
heartVsLiverDF <- as.data.frame(heartVsLiver)
heartVskidneyDF <- as.data.frame(heartVskidney)

heartVsLiverDF <- heartVsLiverDF[!is.na(heartVsLiverDF$padj),]
heartVskidneyDF <- heartVskidneyDF[!is.na(heartVskidneyDF$padj),]

heartVsLiverDF <- heartVsLiverDF[,c("log2FoldChange","padj")]
heartVskidneyDF <- heartVskidneyDF[,c("log2FoldChange","padj")]


## ----tissue6,eval=TRUE,echo=TRUE,fig.height=5,fig.width=7---------------------
colnames(heartVskidneyDF) <- paste0("HeartVsKidney","_",colnames(heartVskidneyDF))
colnames(heartVsLiverDF) <- paste0("HeartVsLiver","_",colnames(heartVsLiverDF))
fullTable <- merge(heartVsLiverDF, heartVskidneyDF, by=0)
fullTable[1:2,]


## ----tissue7,eval=TRUE,echo=TRUE,fig.height=5,fig.width=7---------------------
upInHeart <- fullTable$HeartVsLiver_log2FoldChange > 0 &
             fullTable$HeartVsKidney_log2FoldChange > 0 &
             fullTable$HeartVsLiver_padj < 0.05 &
             fullTable$HeartVsKidney_padj < 0.05
upInHeartTable <- fullTable[upInHeart,]
upInHeartTable[1:2,]


## ----tissue9,eval=TRUE,echo=TRUE,fig.height=5,fig.width=7---------------------
forVenn <- data.frame(UpvsLiver=fullTable$HeartVsLiver_log2FoldChange > 0 &
                       fullTable$HeartVsLiver_padj < 0.05,
                     UpvsKidney=fullTable$HeartVsKidney_log2FoldChange > 0 &
                       fullTable$HeartVsKidney_padj < 0.05)
forVenn[1:3,]


## ----euler,fig.height=3,fig.width=6-------------------------------------------
library(eulerr)
fit <- euler(forVenn)
plot(fit)


## ----eval=T, echo=F-----------------------------------------------------------

# x <- LiverVskidney
# y <- heartVsLiver
# 
# length(which(x$padj<0.05 & 
#         y$padj<0.05 &
#         x$log2FoldChange>0 &
#         y$log2FoldChange<0))
# 
# x <- heartVskidney
# y <- heartVsLiver
# length(which(x$padj<0.05 & 
#         y$padj<0.05 &
#         x$log2FoldChange>0 &
#         y$log2FoldChange>0))
# 
# x <- heartVskidney
# y <- LiverVskidney
# length(which(x$padj<0.05 & 
#         y$padj<0.05 &
#         x$log2FoldChange<0 &
#         y$log2FoldChange<0))





## ----echo=F, warning=FALSE, message=FALSE,fig.height=5,fig.width=7------------
plotPCA(rlog(dds), intgroup="Tissue")



## ----echo=F, warning=FALSE, message=FALSE,fig.height=4,fig.width=5------------
plotPCA(rlog(dds), intgroup="Tissue")



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Visualizing Counts

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Visualizing Counts

---
"    
  )
  
}


## ----gnT,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE----------------------
normLog2Counts <- normTransform(dds)
normLog2Counts


## ----gnTM,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE,fig.height=4,fig.width=7----
matrixOfNorm <- assay(normLog2Counts)
boxplot(matrixOfNorm, las=2, names=c("Heart_1","Heart_2", "Kidney_1","Kidney_2", "Liver_1", "Liver_2"))


## ----gnTMP,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE,fig.height=4,fig.width=7----
library(vsn)
vsn::meanSdPlot(matrixOfNorm)


## ----gRL,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE----------------------
rlogTissue <- rlog(dds)
rlogTissue


## ----gRLM,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE,fig.height=4,fig.width=5----
rlogMatrix <- assay(rlogTissue)
vsn::meanSdPlot(rlogMatrix)


## ----gnTMP2,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE,fig.height=4,fig.width=5----
library(vsn)
vsn::meanSdPlot(matrixOfNorm)


## ----echo=F,fig.height=3,fig.width=7------------------------------------------
ggplot(as.data.frame(LiverVskidney), aes(x=log2FoldChange, y=-log10(padj))) + geom_point() +theme_minimal() + ggtitle("LiverVskidney")



## ----fig.height=3,fig.width=7-------------------------------------------------
ggplot(as.data.frame(LiverVskidney), aes(x=log2FoldChange, y=-log10(padj))) + geom_point() +theme_minimal() + ggtitle("LiverVskidney")



## ----fig.height=3,fig.width=7-------------------------------------------------
library(EnhancedVolcano)

EnhancedVolcano(LiverVskidney,
    lab = rownames(LiverVskidney),
    x = 'log2FoldChange',
    y = 'padj')


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Dimension reduction

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Dimension reduction

---
"    
  )
  
}


## ----gPCA,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE,fig.height=3,fig.width=7----
plotPCA(rlogTissue,
        intgroup="Tissue",
        ntop=nrow(rlogTissue))


## ----gPCA3,eval=TRUE,echo=FALSE,warning=FALSE,message=FALSE,fig.height=4,fig.width=7----
plotPCA(rlogTissue,
        intgroup="Tissue",
        ntop=nrow(rlogTissue))


## ----gPRcomp,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE------------------
pcRes <- prcomp(t(rlogMatrix))
class(pcRes)
pcRes$x[1:2,]


## ----gPRcosmp,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE,fig.height=4,fig.width=7----
plot(pcRes$x,
     col=colData(rlogTissue)$Tissue,
     pch=20,
     cex=2)
legend("top",legend = c("Heart","Kidney","Liver"),
       fill=unique(colData(rlogTissue)$Tissue))


## ----gPRloading,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE---------------
pcRes$rotation[1:5,1:4]


## ----gPRload2,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE-----------------
PC2markers <- sort(pcRes$rotation[,2],decreasing = FALSE)[1:100]
PC2markers[1:10]


## ----gPRcompRot,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE,fig.height=4,fig.width=7----
PC2_hVsl <- heartVsLiver$stat[rownames(heartVsLiver) %in% names(PC2markers)]
PC2_hVsk <- heartVskidney$stat[rownames(heartVskidney) %in% names(PC2markers)]
PC2_LVsk <- LiverVskidney$stat[rownames(LiverVskidney) %in% names(PC2markers)]



## ----gPRcompRot2,eval=FALSE,echo=FALSE,warning=FALSE,message=FALSE------------
# PC1markers <- sort(pcRes$rotation[,1],decreasing = FALSE)[1:100]
# PC1_hVsl <- heartVsLiver$stat[rownames(heartVsLiver) %in% names(PC1markers)]
# PC1_hVsk <- heartVskidney$stat[rownames(heartVskidney) %in% names(PC1markers)]
# PC1_LVsk <- LiverVskidney$stat[rownames(LiverVskidney) %in% names(PC1markers)]
# boxplot(PC1_hVsl,PC1_hVsk,PC1_LVsk,names=c("HeartVsLiver","HeartVsKidney","LiverVsKidney"))


## ----gPRcompRot3,eval=FALSE,echo=FALSE,warning=FALSE,message=FALSE------------
# PC1markers <- sort(pcRes$rotation[,1],decreasing = TRUE)[1:100]
# PC1_hVsl <- heartVsLiver$stat[rownames(heartVsLiver) %in% names(PC1markers)]
# PC1_hVsk <- heartVskidney$stat[rownames(heartVskidney) %in% names(PC1markers)]
# PC1_LVsk <- LiverVskidney$stat[rownames(LiverVskidney) %in% names(PC1markers)]
# boxplot(PC1_hVsl,PC1_hVsk,PC1_LVsk,names=c("HeartVsLiver","HeartVsKidney","LiverVsKidney"))


## ----gPRcompRotf,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE,fig.height=4,fig.width=7----
boxplot(PC2_hVsl,PC2_hVsk,PC2_LVsk,
        names=c("Heart/Liver","Heart/Kidney","Liver/Kidney"),ylab="log2FC")


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Sample-to-Sample correlation

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Sample-to-Sample correlation

---
"    
  )
  
}


## ----gSampleDista,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE-------------
sampleCor <- cor(rlogMatrix)
sampleCor


## ----gSampleDistb,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE-------------

sampleDists <- as.dist(1-cor(rlogMatrix))
sampleDistMatrix <- as.matrix(sampleDists)



## ----gSampleDistc,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE,fig.height=4,fig.width=7----
library(pheatmap)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)


## ----gSampleDistca,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE,fig.height=3,fig.width=7----
library(RColorBrewer)
blueColours <- brewer.pal(9, "Blues")
colors <- colorRampPalette(rev(blueColours))(255)
plot(1:255,rep(1,255),
     col=colors,pch=20,cex=20,ann=FALSE,
     yaxt="n")



## ----gSampleDistd,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE,fig.height=4,fig.width=7----
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         color = colors)


## ----gSampleDiste,eval=T,echo=TRUE,warning=FALSE,message=FALSE,fig.height=3,fig.width=7----
annoCol <- as.data.frame(colData(dds))
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         color = colors,annotation_col = annoCol)


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Batch issues

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Batch issues

---
"    
  )
  
}


## ----eval=FALSE,echo=FALSE, include=FALSE-------------------------------------
# 
# in_mat <- read.delim("~/Downloads/GSE49712_raw_counts_GRCh38.p13_NCBI.tsv", sep ="\t")[,2:11]
# rownames(in_mat) <- read.delim("~/Downloads/GSE49712_raw_counts_GRCh38.p13_NCBI.tsv", sep ="\t")[,1]
# 
# mycoldata <- data.frame(treat = rep(c("A","B"),5),
#                         replicate = c(1,1,2,2,3,3,4,4,5,5))
# rownames(mycoldata) <- colnames(in_mat)
# 
# dds <- DESeqDataSetFromMatrix(in_mat, design=~treat,colData=mycoldata)
# dds <- DESeq(dds)
# 
# myrlog <- rlog(dds)
# 
# plotPCA(myrlog, intgroup="treat")
# 
# 
# plotPCA(myrlog, intgroup="replicate")


## ----eval=FALSE,echo=FALSE, include=FALSE-------------------------------------
# library(magrittr)
# library(tibble)
# library(ggplot2)
# temp <- read.delim("~/Downloads/GSE190411_Yuan2021_PAP_SCC_RNAseq_counts.txt")
# temp$gene <- make.unique(temp$gene)
# temp %<>% column_to_rownames("gene") %>% as.matrix() %>% ceiling()
# CellLine <- grepl("PAP",colnames(temp)) %>% ifelse("PAP","SCC") %>% as.factor()
# Treatment <- grepl("pos",colnames(temp)) %>% ifelse("Pos","Neg") %>% as.factor()
# Batch <- c("B1","B2","B1","B2","B1","B2","B1","B2")
# my_pca <- prcomp(data.matrix(temp))
# my_pca_rot <- as.data.frame(my_pca$rotation)
# my_pca_rot$sample <- rownames(my_pca_rot)
# 
# my_pca_rot$Batch <- Batch
# my_pca_rot$rep <- rep(c(1,2),4)
# my_pca_rot$Treatment <- Treatment
# my_pca_rot$CellLine <- CellLine


## ----eval=FALSE,echo=FALSE, include=FALSE-------------------------------------
# my_pca_rot %>% ggplot(aes(PC1, PC2, color=rep)) + geom_point()
# my_pca_rot %>% ggplot(aes(PC1, PC2, color=Batch)) + geom_point()
# my_pca_rot %>% ggplot(aes(PC1, PC2, color=CellLine)) + geom_point()
# my_pca_rot %>% ggplot(aes(PC1, PC2, color=Treatment)) + geom_point()
# 


## ----eval=FALSE,echo=FALSE, include=FALSE-------------------------------------
# coldata <- data.frame(CellLine=CellLine,
#                       Treatment=Treatment,
#                       Batch=c("Rep1","Rep2","Rep1","Rep2","Rep1","Rep2","Rep1","Rep2"),
#                       row.names = colnames(temp))
# dds <- DESeqDataSetFromMatrix(countData = temp,colData = coldata,design = ~Batch+Treatment)
# myrlog <- rlog(dds)
# a <- plotPCA(myrlog, intgroup="Batch")+ ggtitle("Batch Effects Across Replicates")
# 
# b <- plotPCA(myrlog, intgroup="CellLine")+ ggtitle("Batch Effects Across Cell Line")
# c <- plotPCA(myrlog, intgroup="Treatment")+ ggtitle("Batch Effects Across Treatment")
# 
# a
# b
# c
# # ggsave(a,"batch1.png")
# # ggsave(b,"batch2.png")
# # ggsave(c,"batch3.png")


## ----eval=FALSE,echo=FALSE, include=FALSE-------------------------------------
# 
# test <- read.delim("~/Downloads/GSE190411_Yuan2021_PAP_SCC_RNAseq_counts.txt")
# 
# 
# my_pca <- prcomp(data.matrix(test[,-1]))
# my_pca_rot <- as.data.frame(my_pca$rotation)
# my_pca_rot$sample <- rownames(my_pca_rot)
# 
# my_pca_rot$rep <- rep(c(1,2),4)
# my_pca_rot$neg <- rep(c("neg","neg","pos","pos"),2)
# my_pca_rot$scc <- c(rep("PAP",4),rep("SCC",4))


## ----eval=FALSE,echo=FALSE, include=FALSE-------------------------------------
# my_pca_rot %>% ggplot(aes(PC1, PC2, color=rep)) + geom_point()
# my_pca_rot %>% ggplot(aes(PC1, PC2, color=scc)) + geom_point()
# my_pca_rot %>% ggplot(aes(PC1, PC2, color=neg)) + geom_point()
# 


## ----eval=FALSE,echo=FALSE, include=FALSE-------------------------------------
# my_pca <- prcomp(log2(data.matrix(test[,-1]+1)))
# my_pca_rot <- as.data.frame(my_pca$rotation)
# my_pca_rot$sample <- rownames(my_pca_rot)
# 
# my_pca_rot$rep <- rep(c(1,2),4)
# my_pca_rot$neg <- rep(c("neg","neg","pos","pos"),2)
# my_pca_rot$scc <- c(rep("PAP",4),rep("SCC",4))
# 


## ----eval=FALSE,echo=FALSE, include=FALSE-------------------------------------
# my_pca_rot %>% ggplot(aes(PC1, PC2, color=rep)) + geom_point()
# my_pca_rot %>% ggplot(aes(PC1, PC2, color=scc)) + geom_point()
# my_pca_rot %>% ggplot(aes(PC1, PC2, color=neg)) + geom_point()
# 


## ----eval=T,echo=FALSE, warning=FALSE, message=F------------------------------

library(magrittr)
library(ggplot2)

my_counts <- read.delim("data/GSE190411_Yuan2021_PAP_SCC_RNAseq_counts.txt")
rownames(my_counts) <- make.unique(my_counts$gene)
my_counts <- data.matrix(my_counts[,-1])
my_counts <- ceiling(my_counts)

CellLine <- grepl("PAP",colnames(my_counts)) %>% ifelse("PAP","SCC") %>% as.factor()
Treatment <- grepl("pos",colnames(my_counts)) %>% ifelse("Pos","Neg") %>% as.factor()
Rep <- factor(c("R1","R2","R1","R2","R1","R2","R1","R2"))

coldata <- data.frame(CellLine=CellLine,
                      Treatment=Treatment,
                      Batch=Rep,
                      row.names = colnames(my_counts))
dds <- DESeqDataSetFromMatrix(countData = my_counts,colData = coldata,design =~Treatment)
myrlog <- rlog(dds)

a <- plotPCA(myrlog, intgroup="Batch")+ ggtitle("Batch Effects Across Replicates")
b <- plotPCA(myrlog, intgroup="CellLine")+ ggtitle("Batch Effects Across Cell Line")
c <- plotPCA(myrlog, intgroup="Treatment")+ ggtitle("Batch Effects Across Treatment")


## ----echo=F, eval=T-----------------------------------------------------------
a


## ----echo=F, eval=T-----------------------------------------------------------
b


## ----echo=F, eval=T-----------------------------------------------------------
c


## ----eval=F, echo =T----------------------------------------------------------
# 
# library(ggplot2)
# 
# my_counts <- read.delim("data/GSE190411_Yuan2021_PAP_SCC_RNAseq_counts.txt")
# rownames(my_counts) <- make.unique(my_counts$gene)
# my_counts <- data.matrix(my_counts[,-1])
# my_counts <- ceiling(my_counts)
# head(my_counts)
# 


## ----eval=T, echo=F-----------------------------------------------------------

head(my_counts)


## ----eval=F, echo =T----------------------------------------------------------
# 
# CellLine <- factor(c("PAP","PAP","PAP","PAP","SCC","SCC","SCC","SCC"))
# Treatment <- factor(c("Neg","Neg","Pos","Pos","Neg","Neg","Pos","Pos"))
# Rep <- factor(c("R1","R2","R1","R2","R1","R2","R1","R2"))
# 
# coldata <- data.frame(CellLine=CellLine,
#                       Treatment=Treatment,
#                       Batch=Rep,
#                       row.names = colnames(my_counts))
# head(coldata)


## ----eval=T, echo=F-----------------------------------------------------------

head(coldata)


## ----eval=F, echo =T----------------------------------------------------------
# 
# dds <- DESeqDataSetFromMatrix(countData = my_counts,colData = coldata,design =~Treatment)
# myrlog <- rlog(dds)


## ----fig.height=3,fig.width=5, warning=F, message=F---------------------------
plotPCA(myrlog, intgroup="Batch")+ ggtitle("Batch Effects Across Replicates")



## -----------------------------------------------------------------------------
library(limma)

corrected_matrix <- removeBatchEffect(assay(myrlog), batch=coldata$Batch)




## -----------------------------------------------------------------------------

my_corrected_rlog <- myrlog
assay(my_corrected_rlog) <- corrected_matrix



## ----fig.height=3,fig.width=5, warning=F, message=F---------------------------

plotPCA(my_corrected_rlog, intgroup="Batch")+ ggtitle("Batch Effects Across Replicates - Corrected")


## ----fig.height=3,fig.width=5, warning=F, message=F---------------------------

plotPCA(myrlog, intgroup="Batch")+ ggtitle("Batch Effects Across Replicates - Uncorrected")


## ----fig.height=3,fig.width=5, warning=F, message=F---------------------------

plotPCA(my_corrected_rlog, intgroup="CellLine")+ ggtitle("Batch Effects Across Replicates - Corrected")


## ----fig.height=3,fig.width=5, warning=F, message=F---------------------------

plotPCA(myrlog, intgroup="CellLine")+ ggtitle("Batch Effects Across Replicates - Uncorrected")


## -----------------------------------------------------------------------------

plotPCA(my_corrected_rlog, intgroup="Treatment")+ ggtitle("Batch Effects Across Replicates - Corrected")


## -----------------------------------------------------------------------------

plotPCA(myrlog, intgroup="Treatment")+ ggtitle("Batch Effects Across Replicates - Uncorrected")


## ----warning=F, message=FALSE-------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = my_counts, colData = coldata, design =~Treatment)
dds <- DESeq(dds)


## ----warning=F, message=FALSE-------------------------------------------------
dds_corrected <- DESeqDataSetFromMatrix(countData = my_counts, colData = coldata, design =~Batch+Treatment)
dds_corrected <- DESeq(dds_corrected)



## -----------------------------------------------------------------------------
my_res <- results(dds, contrast = c("Treatment","Pos","Neg" ))
summary(my_res)


## -----------------------------------------------------------------------------
my_res_corrected <- results(dds_corrected, contrast = c("Treatment","Pos","Neg" ))
summary(my_res_corrected)


## ----wanring=F, message=F, eval=F---------------------------------------------
# library(pcaExplorer)
# 
# pcaExplorer(dds = dds, dst = myrlog)


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Clustering Analysis

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Clustering Analysis

---
"    
  )
  
}


## -----------------------------------------------------------------------------
load("data/gC_TissueFull.RData")
dds <- DESeqDataSet(geneCounts, design = ~Tissue)



## ----tissue11,eval=TRUE,echo=TRUE,message=FALSE,warning=FALSE-----------------

dds2 <- DESeq(dds,test="LRT",reduced=~1)
acrossGroups <- results(dds2)
acrossGroups <- acrossGroups[order(acrossGroups$pvalue),]
acrossGroups[1:3,]


## ----tissue12,eval=TRUE,echo=TRUE,fig.width=6,fig.height=5--------------------
plotCounts(dds2, gene="17888", intgroup = "Tissue")


## ----he1,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE----------------------
sigChanges <- rownames(acrossGroups)[acrossGroups$padj < 0.01 & !is.na(acrossGroups$padj)]
sigMat <- rlogMatrix[rownames(rlogMatrix) %in% sigChanges,]
nrow(rlogMatrix)
nrow(sigMat)


## ----he2,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE,fig.height=4,fig.width=7----
library(pheatmap)
pheatmap(sigMat,
         scale="row",
         show_rownames = FALSE)


## ----he2a,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE,fig.height=6,fig.width=7----
pheatmap(sigMat,
         scale="row",
         show_rownames = FALSE)


## ----km1,eval=FALSE,echo=TRUE,warning=FALSE,message=FALSE---------------------
# library(pheatmap)
# set.seed(42)
# k <-   pheatmap(sigMat,
#            scale="row", kmeans_k = 7)


## ----km2,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE,fig.height=4,fig.width=7----
library(pheatmap)
set.seed(42)
k <-   pheatmap(sigMat,
           scale="row",kmeans_k = 7)


## ----km3,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE----------------------
names(k$kmeans)
clusterDF <- as.data.frame(factor(k$kmeans$cluster))
colnames(clusterDF) <- "Cluster"
clusterDF[1:10,,drop=FALSE]


## ----km4,eval=FALSE,echo=TRUE,warning=FALSE,message=FALSE---------------------
# 
# OrderByCluster <- sigMat[order(clusterDF$Cluster),]
# 
# pheatmap(OrderByCluster,
#            scale="row", annotation_row = clusterDF,
#            show_rownames = FALSE, cluster_rows = FALSE)
# 


## ----km4r,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE,fig.height=5,fig.width=7----
OrderByCluster <- sigMat[order(clusterDF$Cluster),]

pheatmap(OrderByCluster,
           scale="row", annotation_row = clusterDF,
           show_rownames = FALSE, cluster_rows = FALSE)



## ----kms5,eval=T,echo=TRUE,warning=FALSE,message=FALSE------------------------
library(NbClust)
rowScaledMat <- t(scale(t(sigMat)))
clusterNum <- NbClust(rowScaledMat,distance = "euclidean",
          min.nc = 2, max.nc = 12, 
          method = "kmeans", index ="silhouette")

clusterNum$Best.nc


## ----kmsh,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE---------------------
clusterNum$Best.partition[1:10]
orderedCluster <- sort(clusterNum$Best.partition)
sigMat <- sigMat[match(names(orderedCluster),rownames(sigMat)),]


## ----kmsha,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE--------------------
pheatmap(sigMat,
           scale="row",annotation_row = clusterDF,
           show_rownames = FALSE,cluster_rows = FALSE)

