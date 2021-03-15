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
## library(KEGG.db)
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


## ---- Readinggmt, eval=TRUE, echo=T-------------------------------------------
hallMarks <- getGmt(con="data/h.all.v7.1.symbols.gmt")
hallMarks


## ---- GeneSetCollection-------------------------------------------------------
hallMarks[1]
hallMarks[[1]]


## ---- names-------------------------------------------------------------------
names(hallMarks)


## ---- geneIDs-----------------------------------------------------------------
geneIds(hallMarks)[1:3]


## ---- gskb--------------------------------------------------------------------
library(gskb)


## ---- data--------------------------------------------------------------------
data(mm_miRNA)
names(mm_miRNA)[1:2]
mm_miRNA[1]


## ---- wehiA-------------------------------------------------------------------
library(msigdbr)

mm_H <- msigdbr(species = "Mus musculus", category = "H")
head(mm_H)


## ---- wehi3-------------------------------------------------------------------
myGeneSetList <- list()
gsets<-unique(mm_H$gs_name)

for(i in 1:length(gsets)){
  genes_in_set <- as.character(mm_H$entrez_gene[which(mm_H$gs_name==gsets[i])])
  myGeneSetList[[i]] <- GeneSet(unique(genes_in_set), 
setName=gsets[i])
}

myGeneSetCollection <- GeneSetCollection(myGeneSetList)
myGeneSetCollection


## ---- wehi4-------------------------------------------------------------------
toGmt(myGeneSetCollection,
      con="mouse_Hallmarks.gmt")


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
pwf = nullp(UpInAct, "mm10", "knownGene", plot.fit = TRUE)


## ----funca,eval=TRUE,echo=FALSE,cache=TRUE,include=FALSE----------------------
load(file="data/fit.RData")


## ----func2,eval=TRUE,echo=TRUE,cache=TRUE,dependson="func1",warning=FALSE,message=FALSE----
GO_UpInAct <- goseq(pwf,"mm10","knownGene",
                       test.cats=c("GO:BP"))
GO_UpInAct[1:3,]


## ----func3,eval=TRUE,echo=TRUE,cache=TRUE,dependson="funca",warning=FALSE,message=FALSE----
library(org.Mm.eg.db)
ImmuneResponseGenes <- AnnotationDbi::select(org.Mm.eg.db,keytype = "GOALL",
                              keys = "GO:0006955",columns = "ENTREZID")
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

forRNK <- data.frame(Activated_minus_Resting$ENTREZID,
                     Activated_minus_Resting$stat)

forRNK[1:3,]


## ---- myRNKwrite--------------------------------------------------------------

write.table(forRNK,
            file="data/Activated_minus_Resting.rnk",
            sep="\t",
            col.names = FALSE,
            row.names = FALSE)


## ---- fgsea-------------------------------------------------------------------
library(fgsea)


## ---- gmtPathways_pub, eval=F, echo=T-----------------------------------------
## mouse_Hallmarks <- gmtPathways("mouse_Hallmarks.gmt")
## class(mouse_Hallmarks)
## names(mouse_Hallmarks)


## ---- gmtPathways, eval=T, echo=F---------------------------------------------
mouse_Hallmarks <- gmtPathways("data/mouse_Hallmarks.gmt")
class(mouse_Hallmarks)
names(mouse_Hallmarks)


## ---- rnkForR-----------------------------------------------------------------
Act_minus_Rest_rnk <- read.delim("data/Activated_minus_Resting.rnk",sep="\t",
                                 h=FALSE,row.names = 1)
Act_minus_Rest_gsea <- Act_minus_Rest_rnk[,1]
names(Act_minus_Rest_gsea) <- rownames(Act_minus_Rest_rnk)
Act_minus_Rest_gsea[1:3]


## ---- fgsefun, warning=F, message=F-------------------------------------------

Act_minus_Rest_gseaRes <- fgsea(mouse_Hallmarks, 
                                Act_minus_Rest_gsea, 
                                minSize=15, maxSize=500, nperm=1000)
Act_minus_Rest_gseaRes <- Act_minus_Rest_gseaRes[order(Act_minus_Rest_gseaRes$NES,
                                                       decreasing = T),]
Act_minus_Rest_gseaRes[1:2,]


## ---- fgsefunews--------------------------------------------------------------
Act_minus_Rest_gseaRes$leadingEdge


## ---- fgsefuneds--------------------------------------------------------------
INTERFERON_Response_LE <- Act_minus_Rest_gseaRes$leadingEdge[[1]]
IR_LE <- Activated_minus_Resting[Activated_minus_Resting$ENTREZID %in% INTERFERON_Response_LE,]
IR_LE[1:2,]


## ---- plotEnrichment, fig.height=4, fig.width=6-------------------------------
plotEnrichment(mouse_Hallmarks[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],
               Act_minus_Rest_gsea)


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# ClusterProfiler for all things GSA

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# ClusterProfiler for all things GSA

---
"    
  )
  
}


## ---- warning=F, message=F----------------------------------------------------
library(clusterProfiler)


## -----------------------------------------------------------------------------

mm_c7 <- msigdbr(species = "Mus musculus", category = "C7")[,c(3,4)]



## ---- wanring=F, message=F----------------------------------------------------

sig_genes <- Activated_minus_Resting[Activated_minus_Resting$padj<0.05,1]

sig_gene_enr <- enricher(na.omit(sig_genes), TERM2GENE=mm_c7)



## -----------------------------------------------------------------------------
sig_gene_enr



## ---- fig.width=10, fig.height=4----------------------------------------------
library(ggplot2)
dotplot(sig_gene_enr, showCategory = 5) + theme( axis.text.y = element_text(size = 7))


