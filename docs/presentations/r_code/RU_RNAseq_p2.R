params <-
list(isSlides = "no")

## ----include=FALSE------------------------------------------------------------
suppressPackageStartupMessages(require(knitr))
knitr::opts_chunk$set(echo = TRUE, tidy = T)

condaSalmon <- CondaSysReqs::install_CondaTools("salmon","salmon")
pathToSalmon <- file.path(dirname(dirname(condaSalmon$pathToConda)),"envs",condaSalmon$environment,"bin","salmon")


## ----neQants,include=FALSE,eval=FALSE-----------------------------------------
## 
## salmonExec <- paste0(pathToSalmon," quant")
## fq <- "ENCFF070QMF.fastq.gz"
## outDir <- "TReg_2_Quant"
## salmonQuantCmd <- paste(salmonExec,
##                         "-i",indexName,
##                         "-o",outDir,
##                         "-l A",
##                         "-r",fq)
## salmonQuantCmd
## system(salmonQuantCmd, wait = TRUE)
## 
## salmonExec <- paste0(pathToSalmon," quant")
## fq <- "ENCFF144YYI.fastq.gz"
## outDir <- "TReg_act_1_Quant"
## salmonQuantCmd <- paste(salmonExec,
##                         "-i",indexName,
##                         "-o",outDir,
##                         "-l A",
##                         "-r",fq)
## salmonQuantCmd
## system(salmonQuantCmd, wait = TRUE)
## 
## salmonExec <- paste0(pathToSalmon," quant")
## fq <- "ENCFF042XBW.fastq.gz"
## outDir <- "TReg_act_2_Quant"
## salmonQuantCmd <- paste(salmonExec,
##                         "-i",indexName,
##                         "-o",outDir,
##                         "-l A",
##                         "-r",fq)
## salmonQuantCmd
## system(salmonQuantCmd, wait = TRUE)
## 
## salmonExec <- paste0(pathToSalmon," quant")
## fq <- "ENCFF053CFZ.fastq.gz"
## outDir <- "TReg_act_3_Quant"
## salmonQuantCmd <- paste(salmonExec,
##                         "-i",indexName,
##                         "-o",outDir,
##                         "-l A",
##                         "-r",fq)
## salmonQuantCmd
## system(salmonQuantCmd, wait = TRUE)
## 
## 


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides != "yes"){
  cat("# RNAseq (part 2)


"    
  )
  
}



## ----setup, include=FALSE-----------------------------------------------------
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicAlignments)
library(DESeq2)
library(tximport)
library(org.Mm.eg.db)


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Counting with multiple RNAseq datasets

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 


"    
  )
}else{
  cat("# Counting with multiple RNAseq datasets


"    
  )
  
}



## ----eval=FALSE,echo=TRUE-----------------------------------------------------
## library(Rsamtools)
## bamFilesToCount <- c("Sorted_Treg_1.bam","Sorted_Treg_2.bam",
##                      "Sorted_Treg_act_1.bam","Sorted_Treg_act_2.bam",
##                      "Sorted_Treg_act_3.bam")
## names(bamFilesToCount) <- c("Sorted_Treg_1","Sorted_Treg_2",
##                      "Sorted_Treg_act_1","Sorted_Treg_act_2",
##                      "Sorted_Treg_act_3")
## myBams <- BamFileList(bamFilesToCount,yieldSize = 10000)
## 


## ----eval=FALSE,echo=TRUE-----------------------------------------------------
## library(TxDb.Mmusculus.UCSC.mm10.knownGene)
## library(GenomicAlignments)
## geneExons <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene,by="gene")
## geneCounts <- summarizeOverlaps(geneExons,myBams,
##                                     ignore.strand = TRUE)
## geneCounts


## ----gC1,eval=TRUE,echo=FALSE,cache=TRUE--------------------------------------
load("data/GeneCounts.RData")
geneCounts


## ----bp,eval=FALSE,echo=TRUE,cache=TRUE---------------------------------------
## library(BiocParallel)


## ----bp1,eval=FALSE,echo=TRUE,cache=TRUE,dependson="bp"-----------------------
## paramMulti <- MulticoreParam(workers=2)
## paramSerial <- SerialParam()
## register(paramSerial)


## ----gC2,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gC1"-----------------------
assay(geneCounts)[1:2,]


## ----gC3,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gC1"-----------------------
rowRanges(geneCounts)[1:2,]


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Differential gene expression analysis

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 


"    
  )
}else{
  cat("# Differential gene expression analysis


"    
  )
  
}



## ----de1,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gC1"-----------------------
metaData <- data.frame(Group=c("Naive","Naive","Act","Act","Act"),
                       row.names = colnames(geneCounts))
metaData


## ----de22,eval=TRUE,echo=TRUE,cache=TRUE,dependson="de1"----------------------
countMatrix <- assay(geneCounts)
countGRanges <- rowRanges(geneCounts)
dds <- DESeqDataSetFromMatrix(countMatrix,
                              colData = metaData,
                              design = ~Group,
                              rowRanges=countGRanges)
dds


## ----de2,eval=TRUE,echo=TRUE,cache=TRUE,dependson="de1"-----------------------
colData(geneCounts)$Group <- metaData$Group
geneCounts


## ----de3,eval=TRUE,echo=TRUE,cache=TRUE,dependson="de1"-----------------------
dds <- DESeqDataSet(geneCounts,design = ~Group)
dds


## ----de4,eval=TRUE,echo=TRUE,cache=TRUE,dependson="de3"-----------------------
dds <- DESeq(dds)


## ----de5,eval=TRUE,echo=TRUE,cache=TRUE,dependson="de4"-----------------------
normCounts  <- counts(dds, normalized=TRUE)
normCounts[1:2,]


## ----de6,eval=TRUE,echo=TRUE,cache=TRUE,dependson="de4",fig.width=8,fig.height=4.5----
plotDispEsts(dds)


## ----de7,eval=TRUE,echo=TRUE,cache=TRUE,dependson="de4"-----------------------
myRes <-results(dds,contrast = c("Group","Act","Naive"))
myRes <- myRes[order(myRes$pvalue),]
myRes[1:3,]


## ----dr1,eval=TRUE,echo=TRUE,cache=TRUE,dependson=c("de7")--------------------
summary(myRes)


## ----drccd1a,eval=FALSE,echo=TRUE---------------------------------------------
## plotMA(myRes)


## ----dr1a,eval=TRUE,echo=FALSE,cache=TRUE,dependson=c("de7"),fig.height=4,fig.width=7----
DESeq2::plotMA(myRes)


## ----dr1b,eval=TRUE,echo=TRUE,cache=TRUE,dependson=c("de7"),fig.height=4,fig.width=7,message=FALSE,warning=FALSE----
myRes_lfc <- lfcShrink(dds, coef =  "Group_Naive_vs_Act")
DESeq2::plotMA(myRes_lfc)


## ----dr1c,eval=TRUE,echo=FALSE,cache=TRUE,dependson=c("dr1b")-----------------
myRes <-results(dds,contrast = c("Group","Act","Naive"))


## ----dr2,eval=TRUE,echo=TRUE,cache=TRUE,dependson=c("de7")--------------------
myResAsDF <- as.data.frame(myRes)
myResAsDF[1:2,]


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Significance and Multiple Testing

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 


"    
  )
}else{
  cat("# Significance and Multiple Testing


"    
  )
  
}



## ----dr22a,eval=TRUE,echo=TRUE,cache=TRUE,dependson=c("de7")------------------
table(is.na(myResAsDF$padj))


## ----dr22b,eval=TRUE,echo=TRUE,cache=TRUE,dependson=c("de7")------------------
myResAsDF$newPadj <- p.adjust(myResAsDF$pvalue)
myResAsDF[1:3,]


## ----dr22,eval=TRUE,echo=TRUE,cache=TRUE,dependson=c("de7")-------------------
myResAsDF <- myResAsDF[!is.na(myResAsDF$padj),]
myResAsDF <- myResAsDF[order(myResAsDF$pvalue), ]
myResAsDF[1:3,]


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# DEseq2 and Salmon

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 


"    
  )
}else{
  cat("# DEseq2 and Salmon


"    
  )
  
}



## ----tx1,eval=TRUE,echo=TRUE,cache=TRUE,dependson="de4"-----------------------
library(tximport)


## ----tx2,eval=TRUE,echo=TRUE,cache=TRUE,dependson="tx1"-----------------------
temp <- read.delim("data/Salmon/TReg_2_Quant/quant.sf")

temp[1:3,]


## ----tx3a,eval=TRUE,echo=TRUE,cache=TRUE,dependson="tx2",warning=FALSE,message=FALSE----

Tx2Gene <- select(TxDb.Mmusculus.UCSC.mm10.knownGene,
                  keys = as.vector(temp[,1]),
                  keytype = "TXNAME",
                  columns = c("GENEID","TXNAME"))
Tx2Gene <- Tx2Gene[!is.na(Tx2Gene$GENEID),]
Tx2Gene[1:10,]


## ----tx4,eval=TRUE,echo=TRUE,cache=TRUE,dependson="tx3a"----------------------
salmonQ <- dir("data/Salmon/",recursive = T,
               pattern = "quant.sf",full.names = T)
salmonCounts <- tximport(salmonQ,
                         type="salmon",
                         tx2gene = Tx2Gene)


## ----tx5,eval=TRUE,echo=TRUE,cache=TRUE,dependson="tx4"-----------------------
salmonCounts$abundance[1:2,]
salmonCounts$counts[1:2,]


## ----tx6,eval=TRUE,echo=TRUE,cache=TRUE,dependson="tx5"-----------------------
ddsSalmon <- DESeqDataSetFromTximport(salmonCounts,
                                      colData = metaData,
                                      design = ~Group)


## ----tx7,eval=TRUE,echo=TRUE,cache=TRUE,dependson=c("tx6"),message=FALSE,warning=FALSE----
ddsSalmon <- DESeq(ddsSalmon)
myResS <-results(ddsSalmon,contrast = c("Group","Act","Naive"))
myResS <- myResS[order(myResS$pvalue),]
myResS[1:3,]


## ----dr1ss,eval=TRUE,echo=FALSE,cache=TRUE,dependson=c("tx7","de7")-----------
swde <- merge(as.data.frame(myRes),as.data.frame(myResS),by=0,all=FALSE)

toCompare <- cbind(abs(swde$log2FoldChange.x) > 1 & swde$padj.x < 0.05 & !is.na(swde$padj.x) & !is.na(swde$padj.y), abs(swde$log2FoldChange.y) > 1 & swde$padj.y < 0.05 & !is.na(swde$padj.y) & !is.na(swde$padj.x))
colnames(toCompare) <- c("summarizeOverlaps","Salmon")
limma::vennDiagram(toCompare)


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Adding Annotation

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 


"    
  )
}else{
  cat("# Adding Annotation


"    
  )
  
}



## ----anno,eval=TRUE,echo=TRUE,cache=TRUE,dependson=c("tx7","de7"),message=FALSE,warning=FALSE,tidy=FALSE----
library(org.Mm.eg.db)
eToSym <- select(org.Mm.eg.db,
                 keys = rownames(myResAsDF),
                 keytype = "ENTREZID",
                 columns="SYMBOL")
eToSym[1:10,]


## ----anno2,eval=TRUE,echo=TRUE,cache=TRUE,dependson="anno",tidy=FALSE---------
annotatedRes <- merge(eToSym,myResAsDF,
                      by.x=1,
                      by.y=0,
                      all.x=FALSE,
                      all.y=TRUE)
annotatedRes <- annotatedRes[order(annotatedRes$pvalue),]
annotatedRes[1:3,]

