params <-
list(isSlides = "no")

## ----include=FALSE------------------------------------------------------------
suppressPackageStartupMessages(require(knitr))
knitr::opts_chunk$set(echo = TRUE, tidy = T)


## ----setup, include=FALSE-----------------------------------------------------
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicAlignments)
library(DESeq2)
library(tximport)
library(org.Mm.eg.db)
library(goseq)


## ---- eval=F------------------------------------------------------------------
## setwd("Path/to/Download/RU_RNAseq-master")
## 
## 


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Gene Sets

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Gene Sets

---
"    
  )
  
}



## ----eval=T,echo=T, eval=F, echo=T, warning=FALSE,tidy=T----------------------
## library(GO.db)
## library(reactome.db)
## library(GSEABase)
## 


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# MSigDB and Gene Set Collections

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# MSigDB and gmt files

---
"    
  )
  
}


## ----gseabase, warning=F, message=F-------------------------------------------
library(GSEABase)
hallMarks <- getGmt(con="data/h.all.v7.1.symbols.gmt")
hallMarks


## ---- GeneSetCollection-------------------------------------------------------
hallMarks[[1]]


## -----------------------------------------------------------------------------
names(hallMarks)


## ---- geneIDs-----------------------------------------------------------------
geneIds(hallMarks)[1:3]


## ---- wehiA-------------------------------------------------------------------
library(msigdbr)

mm_H <- msigdbr(species = "Mus musculus", category = "H")
head(mm_H)


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Testing gene set enrichment

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Testing gene set enrichment

---
"    
  )
  
}


## ---- deres-------------------------------------------------------------------
Activated_minus_Resting <- read.delim(file="data/Group_Activated_minus_Resting.csv",sep=",")
Activated_minus_Resting[1:3,]


## ---- deresFilter-------------------------------------------------------------
Activated_minus_Resting <- Activated_minus_Resting[!is.na(Activated_minus_Resting$padj),]
Activated_minus_Resting[1:3,]



## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Functional enrichment

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Functional enrichment

---
"    
  )
  
}


## ----func,eval=TRUE,echo=TRUE,cache=TRUE,dependson="anno2"--------------------
UpInAct <- Activated_minus_Resting$padj < 0.05 & 
             Activated_minus_Resting$log2FoldChange > 0
UpInAct <- as.integer(UpInAct)
names(UpInAct) <- Activated_minus_Resting$ENTREZID
UpInAct[1:4]
table(UpInAct)


## ----include=FALSE------------------------------------------------------------
library(goseq)


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T----------------------
supGenomes <- supportedGenomes()
supGenomes[1:2,]



## ----func1,eval=TRUE,echo=TRUE,cache=TRUE,dependson="func", warning=F, message=F, fig.height=3.75,fig.width=3.75----
library(goseq)
pwf <- nullp(UpInAct, "mm10", "knownGene", plot.fit = TRUE)


## ----funca,eval=TRUE,echo=FALSE,cache=TRUE,include=FALSE----------------------
load(file = "data/fit.RData")


## ----func2,eval=TRUE,echo=TRUE,cache=TRUE,dependson="func1",warning=FALSE,message=FALSE----
GO_UpInAct <- goseq(pwf,"mm10","knownGene",
                       test.cats=c("GO:BP"))
GO_UpInAct[1:3,]


## ----func3,eval=TRUE,echo=TRUE,cache=TRUE,dependson="funca",warning=FALSE,message=FALSE----
library(org.Mm.eg.db)
ImmuneResponseGenes <- AnnotationDbi::select(org.Mm.eg.db, keytype = "GOALL",
                              keys = "GO:0006955", columns = "ENTREZID")
ImmuneResponseGenes



## ----func4,eval=TRUE,echo=TRUE,cache=TRUE,dependson="func3",warning=FALSE,message=FALSE----
IRG_Entrez <- unique(ImmuneResponseGenes$ENTREZID)
IRG_Res <-  Activated_minus_Resting[Activated_minus_Resting$ENTREZID %in% IRG_Entrez,]
write.table(IRG_Res,
            file="data/ImmuneResponseGeneTable.csv",sep=",",
            row.names = FALSE)
IRG_Res[1:3,]


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# ClusterProfiler

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# ClusterProfiler

---
"    
  )
  
}


## ---- warning=F, message=F----------------------------------------------------
library(clusterProfiler)


## ---- warning=F, message=F----------------------------------------------------

sig_genes <- Activated_minus_Resting[Activated_minus_Resting$padj<0.05,1]
head(sig_genes)


## ---- warning=F, message=F----------------------------------------------------
sig_gene_enr <- enrichGO(sig_genes, OrgDb = org.Mm.eg.db)



## -----------------------------------------------------------------------------
sig_gene_enr



## ---- fig.width=10, fig.height=4, message=FALSE, warning=FALSE----------------
library(ggplot2)
clusterProfiler::dotplot(sig_gene_enr, showCategory = 6) + theme( axis.text.y = element_text(size = 7))



## -----------------------------------------------------------------------------
library(enrichplot)
sig_gene_enr <- pairwise_termsim(sig_gene_enr)
emapplot(sig_gene_enr, showCategory = 15, cex_label_category=0.6) + theme( text = element_text(size = 7))


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# GSEA

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# GSEA

---
"    
  )
  
}


## ---- myRNK-------------------------------------------------------------------

forRNK <- Activated_minus_Resting$stat
names(forRNK) <- Activated_minus_Resting$ENTREZID
  
forRNK <- forRNK[order(forRNK, decreasing = T)]

forRNK[1:6]


## -----------------------------------------------------------------------------

mm_c7 <- msigdbr(species = "Mus musculus", category = "C7")[,c("gs_name","entrez_gene")]
head(mm_c7)


## ---- warning=F, message=F----------------------------------------------------
sig_gene_enr <- GSEA(forRNK, TERM2GENE = mm_c7, eps=1e-100)


## ---- fig.width=10, fig.height=4, message=FALSE, warning=FALSE----------------
clusterProfiler::dotplot(sig_gene_enr, showCategory = 6) + theme( axis.text.y = element_text(size = 7))



## -----------------------------------------------------------------------------
library(enrichplot)
sig_gene_enr <- pairwise_termsim(sig_gene_enr)
emapplot(sig_gene_enr, showCategory = 10, cex_label_category=0.6) + theme( text = element_text(size = 7))


## -----------------------------------------------------------------------------
gseaplot(sig_gene_enr, geneSetID = 1, by = "runningScore", title = "GSE15330_HSC_VS_LYMPHOID_PRIMED_MULTIPOTENT_PROGENITOR_DN")




## -----------------------------------------------------------------------------

write.csv(as.data.frame(sig_gene_enr), "cluster_profiler_GSEA_result.csv")


