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


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides != "yes"){
  cat("# RNAseq (part 4)

---
"    
  )
  
}



## ----eval=F-------------------------------------------------------------------
# setwd("Path/to/Download/RU_RNAseq-master/r_course")
# 
# 


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
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
# library(GO.db)
# library(reactome.db)
# library(GSEABase)
# 


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
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


## ----msigdb, warning=FALSE, message=FALSE-------------------------------------
library(msigdbr)

mm_H <- msigdbr(species = "Mus musculus", category = "H")
head(mm_H)


## ----wehiA, warning=FALSE, message=FALSE--------------------------------------

ss_C5 <- msigdbr(species = "pig", category = "C5")
head(ss_C5)


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
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


## ----deres--------------------------------------------------------------------
Activated_minus_Resting <- read.delim(file="data/Group_Activated_minus_Resting.csv",sep=",")
Activated_minus_Resting[1:3,]


## ----deresFilter--------------------------------------------------------------
Activated_minus_Resting <- Activated_minus_Resting[!is.na(Activated_minus_Resting$padj),]
Activated_minus_Resting[1:3,]



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
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


## ----echo=F-------------------------------------------------------------------

df <- data.frame( "Genes of Interest" = c(25,75),
                  "Gene Universe" = c(50,1950))
rownames(df) <- c("In Gene Set", "Out of Gene Set")

knitr::kable(df, caption = "Contingency Table for Fisher Test")



## ----warning=F, message=F-----------------------------------------------------
library(clusterProfiler)


## ----warning=F, message=F-----------------------------------------------------

sig_genes <- Activated_minus_Resting[Activated_minus_Resting$padj<0.05,1]
head(sig_genes)


## ----warning=F, message=F-----------------------------------------------------
sig_gene_enr <- enrichGO(sig_genes, OrgDb = org.Mm.eg.db)



## -----------------------------------------------------------------------------
sig_gene_enr



## -----------------------------------------------------------------------------
gsa_out <- as.data.frame(sig_gene_enr)
rio::export(gsa_out, "gsa_result.xlsx")

head(gsa_out)


## ----fig.width=10, fig.height=4, message=FALSE, warning=FALSE-----------------
library(ggplot2)
dotplot(sig_gene_enr, showCategory = 6) + theme( axis.text.y = element_text(size = 7))



## ----fig.width=10, fig.height=4, message=FALSE, warning=FALSE-----------------
library(enrichplot)
sig_gene_enr <- pairwise_termsim(sig_gene_enr)
emapplot(sig_gene_enr, showCategory = 15, cex_label_category=0.6) + theme( text = element_text(size = 7))


## ----warning=FALSE, message=FALSE---------------------------------------------

mm_H <- msigdbr(species = "Mus musculus", category = "H")[,c("gs_name","entrez_gene")]
head(mm_H)


## ----warning=F, message=F, fig.width=10, fig.height=4, message=FALSE, warning=FALSE----
sig_gene_enr <- enricher(sig_genes, TERM2GENE = mm_H)
dotplot(sig_gene_enr)



## ----warning=F, message=F, fig.width=10, fig.height=4, message=FALSE, warning=FALSE----
sig_gene_enr <- enricher(sig_genes, TERM2GENE = mm_H, universe = as.character(Activated_minus_Resting$ENTREZID))
dotplot(sig_gene_enr)





## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
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


## ----myRNK--------------------------------------------------------------------

forRNK <- Activated_minus_Resting$stat
names(forRNK) <- Activated_minus_Resting$ENTREZID
  
forRNK <- forRNK[order(forRNK, decreasing = T)]

forRNK[1:6]


## -----------------------------------------------------------------------------

mm_c7 <- msigdbr(species = "Mus musculus", category = "C7")[,c("gs_name","entrez_gene")]
head(mm_c7)


## ----warning=F, message=F-----------------------------------------------------
sig_gene_enr <- GSEA(forRNK, TERM2GENE = mm_c7)


## ----fig.width=10, fig.height=4, message=FALSE, warning=FALSE-----------------
clusterProfiler::dotplot(sig_gene_enr, showCategory = 6) + theme( axis.text.y = element_text(size = 7))



## ----fig.width=10, fig.height=4, message=FALSE, warning=FALSE-----------------
library(enrichplot)
sig_gene_enr <- pairwise_termsim(sig_gene_enr)
emapplot(sig_gene_enr, showCategory = 10, cex_label_category=0.6) + theme( text = element_text(size = 7))


## ----fig.width=10, fig.height=4, message=FALSE, warning=FALSE-----------------
gseaplot(sig_gene_enr, geneSetID = 1, by = "runningScore", title = "GSE15330_HSC_VS_LYMPHOID_PRIMED_MULTIPOTENT_PROGENITOR_DN")



