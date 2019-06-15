## ----load_annotation_db--------------------------------------------------
library(org.Hs.eg.db)
org.Hs.eg.db


## ----summary_annotation--------------------------------------------------
datasets <- ls("package:org.Hs.eg.db")
head(datasets)


## ----help_annotation, eval=FALSE-----------------------------------------
## help(org.Hs.egCHR)


## ----map_chromosome------------------------------------------------------
cls <- columns(org.Hs.eg.db)


## ----map_chromosome_help-------------------------------------------------
help("SYMBOL")


## ----map_symbol_Keys-----------------------------------------------------
keytypes(org.Hs.eg.db)
k <- head(keys(org.Hs.eg.db, keytype="SYMBOL"))
k


## ----map_symbol_select, tidy=FALSE---------------------------------------
mapGenes <- select(org.Hs.eg.db, keys = "APP",
                   columns=c("MAP","ENTREZID","GO"), 
                   keytype = "SYMBOL")

head(mapGenes)


## ----select_GO-----------------------------------------------------------
genesGO <- AnnotationDbi::select(org.Hs.eg.db,
                 keys="GO:0042594",
                 columns=c("ENTREZID","SYMBOL"), 
                 keytype="GO")
head(genesGO)


## ----library_GO, out.width=35--------------------------------------------
library(GO.db)
def <- AnnotationDbi::select(GO.db, 
                             keys="GO:0042594",
                             columns=c("DEFINITION"),
                             keytype="GOID")
head(def)


## ----compute_hypergeometric----------------------------------------------
N <- 2671
m <- 170
n <- 89
k <- 18
sum(choose(m, k:m) * choose(N-m, n-(k:m)) / choose(N, n)) 


## ----go_using_fisher_test------------------------------------------------
t <- array(c(k,n-k,m-k,N+k-n-m), dim=c(2,2),
           dimnames=list(GS=c("in","out"),
                         DE=c("yes","no")))
t
fisher.test(t, alternative="greater")


## ----load_lmfit_libraries------------------------------------------------
library(GEOquery)
library(sva)
library(limma)


## ----genoqueryRDATAch8, eval=FALSE---------------------------------------
## #load transcriptomic data
## gsm.expr <- getGEO("GSE63061")
## gsm.expr <- gsm.expr[[1]]
## exprGEO <- exprs(gsm.expr)
## 
## #get phenotype data
## pheno <- pData(phenoData(gsm.expr))
## status <- pheno$characteristics_ch1
## status <- gsub("status: ","", as.character(status))
## fstatus <- factor(status)
## levels(fstatus)<-gsub(" ", "", levels(fstatus))
## 
## #select variables
## age <- substr(pheno$characteristics_ch1.2, 6,7)
## age<-as.numeric(age)
## sex <- pheno$characteristics_ch1.3
## phenodat<-data.frame(fstatus, age, sex)
## 
## #build models
## mod0 <- model.matrix( ~ age+sex, data = phenodat)
## mod <- model.matrix( ~ fstatus+age+sex, data = phenodat)
## 
## #compute surrogate variables for batch effects
## svobj <- sva(exprGEO, mod, mod0, n.sv=2)
## design <- model.matrix(~ 0+fstatus+sex+age+svobj$sv[,1]+svobj$sv[,2])
## colnames(design) <- c(levels(fstatus),"age","sex", "sva1","sva2")
## 
## #fit the model for the desired contrast
## fit <- lmFit(exprGEO, design)
## contrast.matrix <- makeContrasts(AD-CTL, MCI-CTL, AD-MCI, levels=design)
## 
## fit2 <- contrasts.fit(fit, contrast.matrix)
## fit2 <- eBayes(fit2)


## ----loadresults, echo=FALSE---------------------------------------------
library(GEOquery)
load("GSE63061.Rdata")
load("fit2.RData")
gsm.expr <- gsm.expr[[1]]


## ----pvalsresultsch6-----------------------------------------------------
results<-fit2$p.value
head(results)

pvalues<-results[,"AD - CTL"]
head(pvalues)


## ----significantresults--------------------------------------------------
pAdj <- p.adjust(pvalues, "fdr")
psig <- pAdj[pAdj < 0.01]
head(psig)


## ----probe2Gene----------------------------------------------------------
genesIDs <- as.character(fData(gsm.expr)$ILMN_Gene)
names(genesIDs)<-rownames(gsm.expr)

psigAnnot<-data.frame(genesIDs=genesIDs[names(psig)], 
                      pvalAdj=psig,stringsAsFactors=FALSE)

head(psigAnnot)


## ----get_gene_universe---------------------------------------------------
genes <- keys(org.Hs.eg.db, keytype="SYMBOL")

geneUniverse <- select(org.Hs.eg.db, keys = genes, 
                     columns=c("ENTREZID"), keytype = "SYMBOL")

head(geneUniverse)


## ----get_DEgenes---------------------------------------------------------
mappedgenes <- psigAnnot$genesIDs
mappedgenes <- intersect(mappedgenes,geneUniverse$SYMBOL)
selmappedgenes <- geneUniverse$SYMBOL%in%mappedgenes
mappedgenesIds <- geneUniverse$ENTREZID[selmappedgenes]
head(mappedgenesIds)


## ----load_GOstats, tidy=FALSE--------------------------------------------
library(GOstats)

params <- new("GOHyperGParams", geneIds=mappedgenesIds,
              universeGeneIds=geneUniverse$ENTREZID,
              annotation="org.Hs.eg.db", ontology="BP",
              pvalueCutoff=0.05, conditional=FALSE,
              testDirection="over")


## ----hyperG_test_over----------------------------------------------------
hgOver <- hyperGTest(params)
hgOver
summary(hgOver)[c(1:5),]


## ----hyperG_test_over_unconditional--------------------------------------
conditional(params) <- TRUE
hgOverCond <- hyperGTest(params)
summary(hgOver)[c(1:3),]


## ----hyperG_test_over_conditional, eval=FALSE----------------------------
## htmlReport(hgOverCond, file="goBPcond.html")


## ----GSEAKEGGmaphuman----------------------------------------------------
#?org.Hs.egPATH
frame <- toTable(org.Hs.egPATH)
keggframeData <- data.frame(frame$path_id, frame$gene_id)
head(keggframeData)


## ----GSEAKEGGgenesetcolection--------------------------------------------
library(GSEABase)
library(KEGG.db)

keggFrame <- KEGGFrame(keggframeData, organism="Homo sapiens")
gsc.KEGG <- GeneSetCollection(keggFrame, setType = KEGGCollection())

gsc.KEGG


## ----GSEAKEGGrun---------------------------------------------------------
KEGG.params.bp <-  GSEAKEGGHyperGParams(name="KEGG",
                geneSetCollection=gsc.KEGG,                  
                geneIds=mappedgenesIds,
                universeGeneIds=geneUniverse$ENTREZID, 
                pvalueCutoff=0.05, 
                testDirection="over")

KEGG.results.bp <- hyperGTest(KEGG.params.bp)
head(summary(KEGG.results.bp))


## ----ClusterProfiler-----------------------------------------------------
library(clusterProfiler)

res.enrich <- enrichKEGG(gene = mappedgenesIds,
                 organism = 'hsa',
                 pvalueCutoff = 0.05)

res.enrich[1:5, 1:7]


## ----pointenrichment, fig.cap="Dotplot for enrichment analysis of top DE genes obtained from transcriptomic data analysis of Alzheimer's disease.", message=FALSE, fig.height=5, fig.width=6----
dotplot(res.enrich)


## ----enricher------------------------------------------------------------
gset <- read.gmt("h.all.v6.2.entrez.gmt")
res.enrich.H <- enricher(gene = mappedgenesIds,
                         TERM2GENE = gset)
res.enrich.H[1:5, 3:5]


## ----load_brca_sig, echo=FALSE-------------------------------------------
load("brcaSigCNV.Rdata")


## ----show_brca_sig-------------------------------------------------------
brca.gr.sig


## ----load_hub, eval=FALSE------------------------------------------------
## library(AnnotationHub)
## ah <- AnnotationHub()
## ahDb <- query(ah, pattern = c("Homo sapiens", "EnsDb"))
## ahDb


## ----get_coord, eval=FALSE-----------------------------------------------
## ahEdb <- ahDb[["AH60977"]]
## hg.genes <- genes(ahEdb)


## ----load_hg_genes, echo=FALSE-------------------------------------------
load("hgGenes.Rdata")


## ----show_seq_names------------------------------------------------------
seqnames(hg.genes)
seqnames(brca.gr.sig)


## ----get_genes, eval=FALSE-----------------------------------------------
## seqlevels(hg.genes) <- paste0("chr", seqlevels(hg.genes))


## ----sel_chr1------------------------------------------------------------
sel.protein <- hg.genes[seqnames(hg.genes)%in%c("chr1")]
sel.protein <- sel.protein[sel.protein$gene_biotype == "protein_coding"]
sel.cnvs <- brca.gr.sig[seqnames(brca.gr.sig)%in%c("chr1")]
length(sel.cnvs)


## ----test_cnvs_enrich----------------------------------------------------
library(regioneR)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
res.prot <- overlapPermTest(A=sel.cnvs, B=sel.protein, 
                            mask=NA, genome="hg19", ntimes=100,
                            per.chromosome=TRUE, count.once=TRUE)
res.prot


## ----enrich_promoters_get, eval=FALSE------------------------------------
## file.prom <- "http://gattaca.imppc.org/regioner/data/UCSC.promoters.hg19.bed"
## promoters <- toGRanges(file.prom)
## sel.prom <- promoters[seqnames(promoters) %in% c("chr1")]


## ----enrich_promoters_get2, echo=FALSE-----------------------------------
load("selProm.Rdata")


## ----enrich_promoters----------------------------------------------------
res.prom <- overlapPermTest(sel.cnvs, sel.prom, 
                            ntimes=100, genome="hg19", 
                            count.once=TRUE)
res.prom


## ----length_genes--------------------------------------------------------
length(mappedgenes)


## ----load_CTDquerier, cache=FALSE----------------------------------------
library( CTDquerier )


## ----query, eval=FALSE---------------------------------------------------
## genesAD <- query_ctd_gene( terms = mappedgenes,
##                            verbose = TRUE )

## ----queryLoad, echo=FALSE, cache=FALSE----------------------------------
data(genesAD, package="brgedata")


## ----CTDgenes, message=FALSE, fig.cap='Genes found in CTD that were obtained from the analysis of GEO data GSE63061. Seven from 514 are not present in CTD.', fig.width=3, fig.height=3----
library( ggplot2 )
plot( genesAD ) + ggtitle( "Number of genes" )


## ----genesAD_lost--------------------------------------------------------
get_terms(genesAD)[[ "lost" ]]


## ----az_show_2-----------------------------------------------------------
genesAD


## ----az_gda_all----------------------------------------------------------
diseases <- get_table(genesAD, index_name = "diseases" )
colnames(diseases)
dim(diseases)


## ----az_diseases_unique--------------------------------------------------
length( unique(diseases$Disease.Name))


## ----table_direct_evidence-----------------------------------------------
table(diseases$Direct.Evidence)


## ----az_diseases_curated-------------------------------------------------
mask <- !is.na(diseases$Direct.Evidence) &
         diseases$Direct.Evidence != ""
diseases_cu <- diseases[mask,]
dim(diseases_cu )
length(unique(diseases_cu$Disease.Name))


## ----genes_diseases_az---------------------------------------------------
adgenes <- diseases[diseases$Disease.Name == "Alzheimer Disease" ,]
o <- order(as.numeric(adgenes$Inference.Score), decreasing = TRUE)
adgenes.o <- adgenes[o, c(1,5:8)]
head(adgenes.o)


## ----air_ctd-------------------------------------------------------------
air <- query_ctd_chem( terms = "Air Pollutants" )
air


## ----load_hugo-----------------------------------------------------------
hgnc_universe <- read.delim(system.file( "extdata", "HGNC_Genes.tsv",
                                         package="CTDquerier" ),
                            sep = "\t", stringsAsFactor = FALSE )


## ----gala_enrich_air-----------------------------------------------------
ans.air <- enrich( genesAD, air,
                   universe = hgnc_universe$Approved.Symbol,
                   use = "all" )
ans.air


## ----clean_ch9, echo=FALSE, results='hide'-------------------------------
#remove variables
rm(list=ls())

#garbage collection
gc(verbose = FALSE)

