params <-
list(isSlides = "no")

## ----include=FALSE------------------------------------------------------------
suppressPackageStartupMessages(require(knitr))
knitr::opts_chunk$set(echo = TRUE, tidy = T)

# condaSalmon <- CondaSysReqs::install_CondaTools("salmon","salmon")
# pathToSalmon <- file.path(dirname(dirname(condaSalmon$pathToConda)),"envs",condaSalmon$environment,"bin","salmon")



## ---- echo=F, eval=F----------------------------------------------------------
## if(!file.exists("~/ENCFF332KDA.fastq.gz")){
##   download.file("https://www.encodeproject.org/files/ENCFF332KDA/@@download/ENCFF332KDA.fastq.gz","~/ENCFF332KDA.fastq.gz")
## }
## # if(!file.exists("~/ENCFF070QMF.fastq.gz")){
## #   download.file("https://www.encodeproject.org/files/ENCFF070QMF/@@download/ENCFF070QMF.fastq.gz","~/ENCFF070QMF.fastq.gz")
## # }
## 


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides != "yes"){
  cat("# RNAseq (part 1)

---
"    
  )
  
}



## ---- eval=F------------------------------------------------------------------
## setwd("Path/to/Download/RU_RNAseq-master")
## 
## 


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Working with raw RNAseq data

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Working with raw RNAseq data

---
"    
  )
  
}



## ----shortreada,include=FALSE-------------------------------------------------
library(ShortRead)
library(Rfastp)
library(ggplot2)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)


## ----tregFQ1, eval=F, echo=T--------------------------------------------------
## 
## fq<-"https://www.encodeproject.org/files/ENCFF332KDA/@@download/ENCFF332KDA.fastq.gz"
## download.file(fq, "ENCFF332KDA.fastq.gz")
## 


## -----------------------------------------------------------------------------
json_report <- rfastp(read1 = "data/ENCFF332KDA_sampled.fastq.gz", outputFastq = "ENCFF332KDA_rfastp")


## ----tregShow,eval=T, echo=F--------------------------------------------------
dir(pattern = "ENCFF332KDA_rfastp")


## -----------------------------------------------------------------------------
qcSummary(json_report)


## -----------------------------------------------------------------------------

curvePlot(json_report)


## -----------------------------------------------------------------------------
curvePlot(json_report, curve="content_curves")



## -----------------------------------------------------------------------------
?rfastp



## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Aligning RNAseq data

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Aligning RNAseq data

---
"    
  )
  
}



## ----fa1q, include=FALSE------------------------------------------------------
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)


## ----fa1, echo=TRUE-----------------------------------------------------------
library(BSgenome.Mmusculus.UCSC.mm10)
BSgenome.Mmusculus.UCSC.mm10


## ----fa2,cache=FALSE,echo=TRUE------------------------------------------------
mainChromosomes <- paste0("chr",c(1:19,"X","Y","M"))
mainChrSeq <- lapply(mainChromosomes,
                     function(x) BSgenome.Mmusculus.UCSC.mm10[[x]])
names(mainChrSeq) <- mainChromosomes
mainChrSeqSet <- DNAStringSet(mainChrSeq)
mainChrSeqSet


## ----fa3, echo=TRUE,eval=FALSE------------------------------------------------
## 
## writeXStringSet(mainChrSeqSet,
##                 "BSgenome.Mmusculus.UCSC.mm10.mainChrs.fa")


## ---- echo=TRUE,eval=FALSE----------------------------------------------------
## library(Rsubread)
## buildindex("mm10_mainchrs",
##            "BSgenome.Mmusculus.UCSC.mm10.mainChrs.fa",
##            memory=8000,
##            indexSplit=TRUE)
## 


## ---- echo=TRUE, eval=TRUE----------------------------------------------------
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
myExons <- exons(TxDb.Mmusculus.UCSC.mm10.knownGene,columns=c("tx_id","gene_id"))
myExons <- myExons[lengths(myExons$gene_id) == 1]
myExons


## ---- echo=TRUE, eval=TRUE----------------------------------------------------
dfExons <- as.data.frame(myExons)
SAF <- data.frame(GeneID=dfExons$gene_id,
                  Chr=dfExons$seqnames,
                  Start=dfExons$start,
                  End=dfExons$end,
                  Strand=dfExons$strand)


## ---- echo=TRUE, eval=FALSE---------------------------------------------------
## 
## myMapped <- subjunc("mm10_mainchrs",
##                     "ENCFF332KDA.fastq.gz",
##                     output_format = "BAM",
##                     output_file = "Treg_1.bam",
##                     useAnnotation = TRUE,
##                     annot.ext = SAF,
##                     isGTF=FALSE,
##                     nthreads = 4)
## 


## ----sortindex, echo=TRUE, eval=FALSE-----------------------------------------
## library(Rsamtools)
## sortBam("Treg_1.bam","Sorted_Treg_1")
## indexBam("Sorted_Treg_1.bam")


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Counting with aligned RNAseq data

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Counting with aligned RNAseq data

---
"    
  )
  
}



## ----gm, echo=TRUE,eval=TRUE,cache=TRUE---------------------------------------
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
geneExons <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene,by="gene")
class(geneExons)


## ----gm2,echo=TRUE,eval=TRUE,cache=TRUE,dependson="gm"------------------------
geneExons[1:2]


## ----gm0,echo=TRUE,eval=F-----------------------------------------------------
## library(GenomicAlignments)
## myBam <- BamFile("Sorted_Treg_1.bam",yieldSize = 10000)
## tregGeneCounts <- summarizeOverlaps(geneExons,myBam,
##                                     ignore.strand = TRUE)
## tregGeneCounts


## ----gm3,echo=FALSE,eval=TRUE,cache=TRUE--------------------------------------
load("data/myTregCounts2020.RData")
tregGeneCounts


## ----em2,echo=TRUE,eval=TRUE,cache=TRUE,dependson="em"------------------------
library(GenomicFeatures)
nonOverlappingExons <- disjointExons(TxDb.Mmusculus.UCSC.mm10.knownGene)
names(nonOverlappingExons) <- paste(mcols(nonOverlappingExons)$gene_id,
      mcols(nonOverlappingExons)$exonic_part,
                                    sep="_")
nonOverlappingExons[1:3,]


## ----em3,echo=TRUE,eval=FALSE,cache=TRUE,dependson="em2"----------------------
## tregExonCounts <- summarizeOverlaps(nonOverlappingExons,
##                                     myBam,
##                                     ignore.strand = TRUE,
##                                     inter.feature=FALSE)


## ----em4A,echo=FALSE,eval=TRUE,cache=TRUE-------------------------------------
load("data/myTregExonCounts2020.RData")


## ----em4,echo=TRUE,eval=TRUE,cache=TRUE,dependson=c("em4A","gm3")-------------
geneCounts <- assay(tregGeneCounts)
exonCounts <- assay(tregExonCounts)
head(exonCounts,2)
head(geneCounts,2)


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Transcript quantification with pseudo-alignment

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Transcript quantification with pseudo-alignment

---
"    
  )
  
}



## ---- eval=F------------------------------------------------------------------
## BiocManager::install("Herper")
## library(Herper)
## 


## ---- echo=T, eval=F----------------------------------------------------------
## salmon_paths <- install_CondaTools(tools="salmon", env="RNAseq_analysis")
## salmon_paths


## ---- eval=F, echo=F----------------------------------------------------------
## tempdir2 <- function() {
##     tempDir <- tempdir()
##     if(dir.exists(tempDir)){
##       tempDir <- file.path(tempDir,"rr")
##     }
##     tempDir <- gsub("\\", "/", tempDir, fixed = TRUE)
##     tempDir
## }
## 
## myMiniconda <- file.path(tempdir2(), "Test")
## install_CondaTools(tools="salmon", env="RNAseq_analysis", pathToMiniConda = myMiniconda)
## 


## ----sal1,echo=TRUE,eval=F,cache=TRUE-----------------------------------------
## allTxSeq <- extractTranscriptSeqs(BSgenome.Mmusculus.UCSC.mm10,
##                       TxDb.Mmusculus.UCSC.mm10.knownGene,
##                       use.names=TRUE)
## allTxSeq


## ----sal2,echo=F,eval=F,cache=TRUE,dependson="sal1"---------------------------
## writeXStringSet(allTxSeq,"mm10Trans.fa")


## ----sal0,echo=TRUE,eval=F----------------------------------------------------
## writeXStringSet(allTxSeq,
##                 "mm10Trans.fa")


## ----sal2t,echo=T,eval=F,dependson="sal1",  cache.lazy = FALSE----------------
## mainChromosomes <- paste0("chr",c(1:19,"X","Y","M"))
## mainChrSeq <- lapply(mainChromosomes,
##                      function(x)BSgenome.Mmusculus.UCSC.mm10[[x]])
## names(mainChrSeq) <- mainChromosomes
## mainChrSeqSet <- DNAStringSet(mainChrSeq)
## gentrome <- c(allTxSeq,mainChrSeqSet)


## ----sal0t,echo=TRUE,eval=F---------------------------------------------------
## writeXStringSet(gentrome,
##                 "mm10Gentrome.fa")
## write.table(as.data.frame(mainChromosomes),"decoy.txt",
##             row.names = FALSE,
##             col.names = FALSE,
##             quote=FALSE)


## ----salI_1,echo=TRUE,eval=F,dependson="sal1", warning=F----------------------
## 
## fastaTx <- "mm10Trans.fa"
## indexName <- "mm10Trans"
## 
## with_CondaEnv("RNAseq_analysis",
##                       system2(command="salmon",args = c("index",
##                         "-i",indexName,
##                         "-t",fastaTx),
## 
##                         stdout = TRUE))
## 


## ---- eval=F, echo=F, message=F, warning=F------------------------------------
## 
## with_CondaEnv("RNAseq_analysis",
##                       system2(command="salmon",args = c("quant", "--help-alignment"),
##                         stdout = TRUE),
##               pathToMiniConda = myMiniconda)
## 
## 


## ----salI_1t,echo=TRUE,eval=F-------------------------------------------------
## 
## with_CondaEnv("RNAseq_analysis",
##                       system2(command="salmon",args = c("index",
##                         "-i",indexName,
##                         "-t",fastaTx,
##                         "-d decoy.txt"),
## 
##                         stdout = TRUE))
## 


## ----salQ_1,echo=TRUE,eval=F--------------------------------------------------
## fq <- "ENCFF332KDA.fastq.gz"
## outDir <- "TReg_1_Quant"
## 
## with_CondaEnv("RNAseq_analysis",
##                       system2(command="salmon",args = c("quant",
##                         "-i",indexName,
##                         "-o",outDir,
##                         "-l A",
##                         "-r",fq),
## 
##                         stdout = TRUE))
## 


## ---- eval=F, echo=F----------------------------------------------------------
## 
## fq <- "data/ENCFF332KDA_sampled.fastq.gz"
## outDir <- "TReg_1_Quant"
## 
## with_CondaEnv("RNAseq_analysis",
##                       system2(command="salmon",args = c("quant",
##                         "-i",indexName,
##                         "-o",outDir,
##                         "-l A",
##                         "-r",fq),
## 
##                         stdout = TRUE),
##               pathToMiniConda = myMiniconda)


## ----include=FALSE,eval=FALSE-------------------------------------------------
## salmonExec <- paste0(pathToSalmon," quant")
## fq <- "ENCFF001LDC.fastq.gz"
## outDir <- "Liver1"
## salmonQuantCmd <- paste(salmonExec,
##                         "-i",indexName,
##                         "-o",outDir,
##                         "-l A",
##                         "-r",fq)
## salmonQuantCmd
## system(salmonQuantCmd, wait = TRUE)
## 
## salmonExec <- paste0(pathToSalmon," quant")
## fq <- "ENCFF001LCY.fastq.gz"
## outDir <- "Liver2"
## salmonQuantCmd <- paste(salmonExec,
##                         "-i",indexName,
##                         "-o",outDir,
##                         "-l A",
##                         "-r",fq)
## salmonQuantCmd
## system(salmonQuantCmd, wait = TRUE)
## 
## salmonExec <- paste0(pathToSalmon," quant")
## fq <- "ENCFF001LCD.fastq.gz"
## outDir <- "Heart1"
## salmonQuantCmd <- paste(salmonExec,
##                         "-i",indexName,
##                         "-o",outDir,
##                         "-l A",
##                         "-r",fq)
## salmonQuantCmd
## system(salmonQuantCmd, wait = TRUE)
## 
## salmonExec <- paste0(pathToSalmon," quant")
## fq <- "ENCFF001LCE.fastq.gz"
## outDir <- "Heart2"
## salmonQuantCmd <- paste(salmonExec,
##                         "-i",indexName,
##                         "-o",outDir,
##                         "-l A",
##                         "-r",fq)
## salmonQuantCmd
## system(salmonQuantCmd, wait = TRUE)


## ----salReadIn,echo=TRUE,eval=FALSE,cache=TRUE--------------------------------
## myQuant <- read.delim("TReg_1_Quant/quant.sf")
## myQuant[1:3,]


## ----salReadIn2,echo=FALSE,eval=TRUE,cache=TRUE-------------------------------
myQuant <- read.delim("data/TReg_1_Quant/quant.sf")
myQuant[1:3,]

