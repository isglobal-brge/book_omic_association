## ----GetDatamethylationGEOch7, eval=FALSE--------------------------------
## library(GEOquery)
## gsm.meth <- getGEO("GSE80970", destdir = ".")
## gsm.meth <- gsm.meth[[1]]


## ----load_minfi, echo=FALSE, cache=FALSE---------------------------------
library(minfi)

## ----eSet_GRS, eval=FALSE------------------------------------------------
## library(minfi)
## betas <- exprs(gsm.meth)
## pheno <- pData(gsm.meth)
## gsm.grs <- makeGenomicRatioSetFromMatrix(mat=betas,
##                                          pData=pheno)


## ----load2, eval=FALSE---------------------------------------------------
## gsm.grs <- dropMethylationLoci(gsm.grs)
## gsm.grs <- dropLociWithSnps(gsm.grs)
## rr <- seqnames(rowRanges(gsm.grs))
## gsm.grs.f  <- gsm.grs[!rr%in%c("chrX", "chrY"),]


## ----load_data_filtered, echo=FALSE--------------------------------------
load(file="gsmGRSfil.Rdata")


## ----get_status----------------------------------------------------------
status <- gsm.grs.f$characteristics_ch1.1
status <- gsub("disease status: ","", as.character(status))
gsm.grs.f$status <- as.factor(status)
table(gsm.grs.f$status)


## ----n_na----------------------------------------------------------------
sum(is.na(getBeta(gsm.grs.f)))


## ----noNA, cache=FALSE---------------------------------------------------
gsm.grs.noNA <-  gsm.grs.f[rowSums(is.na(getBeta(gsm.grs.f))) == 0,]
nrow(gsm.grs.f)
nrow(gsm.grs.noNA)


## ----rm, cache=FALSE, echo=FALSE, results='hide'-------------------------
rm(gsm.grs.f, gsm.grs.f1, gsm.grs.f2)
gc()


## ----load_MEAL, cache=FALSE----------------------------------------------
library(MEAL)


## ----dmp_analysis--------------------------------------------------------
methRes <- runPipeline(set = gsm.grs.noNA, 
                       variable_names = "status")
methRes
names(methRes)


## ----get_res_mean--------------------------------------------------------
head(getAssociation(methRes, rid="DiffMean"))


## ----get_res_dmrcate-----------------------------------------------------
head(getAssociation(methRes, rid="dmrcate"))


## ----get_CpGs_by_gene, eval=FALSE----------------------------------------
## getGeneVals(methRes, "ARMS2", genecol = "UCSC_RefGene_Name",
##             fNames = c("chromosome", "start"))


## ----qqplot_dmps---------------------------------------------------------
plot(methRes, rid = "DiffMean", type = "qq")


## ----dmp_analysis_sva, eval = FALSE--------------------------------------
## methRes <- runPipeline(set = gsm.grs.noNA,
##                          variable_names = "status",
##                          sva = TRUE)


## ----ManhattanMean, fig.cap="Manhattan plot corresponding to the comparison of mean CpG values between cases and controls belonging to the Alzheimer's disease dataset (GEO number GSE80970). Top horizontal lines stands for two different significant levels."----
targetRange <- GRanges("chr19:45350000-45550000")
plot(methRes, rid = "DiffMean", type = "manhattan", 
     main = "Differences in Means", 
     highlight = targetRange)


## ----ManhattanVar, fig.cap="Manhattan plot corresponding to the comparison of variance CpG values between cases and controls belonging to the Alzheimer's disease dataset (GEO number GSE80970). Top horizontal lines stands for two different significant levels."----
plot(methRes, rid = "DiffVar", type = "manhattan", 
     main = "Differences in Variances", 
     highlight = targetRange)


## ----RDA_region----------------------------------------------------------
targetRange <- GRanges("chr19:45350000-45550000")
methRange <- runRDA(set = gsm.grs.noNA, 
                    model = ~ status,
                    range = targetRange)


## ----plotRDA, fig.cap="Methylation regional analysis (chr19:45.35-45.55 Mb) of Alzheimer case/control analysis (GEO number GSE80970)"----
plotRDA(object = methRange,
        pheno = as.factor(gsm.grs.noNA$status))


## ----getRDAresults-------------------------------------------------------
getRDAresults(methRange)


## ----RegionalPlot1, fig.cap='DMR analysis (top track after gene annotation) and DMP analysis for means (DiffMean) and variance (DiffVar) of Alzheimer case/control analysis (GEO number GSE80970). The coefficients (lines) and the $P$-values (dots) are also depicted.'----
targetRange <- GRanges("chr19:45350000-45550000")
plotRegion(rset = methRes, range = targetRange)


## ----get_gsmExpr, echo=FALSE---------------------------------------------
load("GSE63061.Rdata")
gsm.expr <- gsm.expr[[1]]


## ----create_proper_data--------------------------------------------------
mask <- fData(gsm.expr)$Chromosome%in%c(1:22)
gsm.expr <- gsm.expr[mask,]

chr <- fData(gsm.expr)$Chromosome 
fData(gsm.expr)$Chromosome <- paste0("chr", chr)

chrPos <- as.character(fData(gsm.expr)$Probe_Coordinates)
ff <- function(x, pos){
  ifelse(length(x)>1, x[pos], 1e10)
}

temp <- sapply(strsplit(chrPos, "-"), FUN=ff, pos=1)
chrStart <- sapply(strsplit(temp, ":"), "[[", 1)

temp <- sapply(strsplit(chrPos, "-"), FUN=ff, pos=2)
chrEnd <- sapply(strsplit(temp, ":"), "[[", 1)              
fData(gsm.expr)$start <- as.numeric(chrStart)
fData(gsm.expr)$end <- as.numeric(chrEnd)


## ----mealExpr------------------------------------------------------------
status <- pData(gsm.expr)$characteristics_ch1
status <- gsub("status: ","", as.character(status))
pData(gsm.expr)$status <- relevel(as.factor(status),3)
table(gsm.expr$status)
exprRes <- runPipeline(gsm.expr, variable_names = "status",
                       sva=TRUE)


## ----plotExprAD----------------------------------------------------------
library(ggplot2)
plot(exprRes, rid = "DiffMean", type = "volcano") + 
  ggtitle("Differences in Means")


## ----plotEpiTrans, fig.cap="Epigenomic (DMRs and DMPs) and transcriptomic results obtained from Alzheimer's disease."----
targetRange <- GRanges("chr19:45350000-45550000")
plotRegion(rset = methRes, rset2 = exprRes, 
           range = targetRange, 
           results = c("DiffMean", "bumphunter"),
           fNames = c("chromosome", "start", "end"),
           fNames2 = c("Chromosome", "start", "end"))


## ----install_meffil, eval=FALSE------------------------------------------
## library(devtools)
## install_github("perishky/meffil")


## ----list_references_cell------------------------------------------------
library(meffil)
meffil.list.cell.type.references()


## ----cell_count_Alzheimer------------------------------------------------
betas <- getBeta(gsm.grs.noNA)
cell.counts <- meffil.estimate.cell.counts.from.betas(betas,
                                  "blood gse35069 complete")
head(cell.counts)


## ----run_Alzheimer_cell--------------------------------------------------
pData(gsm.grs.noNA) <- cbind(pData(gsm.grs.noNA),
                             cell.counts) #add cell count variables
model.cell <- ~ status + Bcell + CD4T + Eos + Mono +
                Neu + NK
methRes.cell <- runDiffMeanAnalysis(gsm.grs.noNA, 
                                    model = model.cell)


## ----res_cell------------------------------------------------------------
head(getAssociation(methRes.cell, rid="DiffMean"))


## ----clean_ch7, echo=FALSE, results='hide'-------------------------------
#remove variables
rm(list=ls())

#garbage collection
gc(verbose = FALSE)

