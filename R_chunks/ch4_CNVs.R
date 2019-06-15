## ----install_gada, eval=FALSE--------------------------------------------
## devtools::install_github("isglobal-brge/R-GADA")


## ----copy_data-----------------------------------------------------------
ss1 <- system.file("extdata/madData", package="brgedata")
dir.create("rawData")
ss2 <- "rawData"
files <- list.files(ss1)
file.copy(file.path(ss1,files), ss2)


## ----show_wd-------------------------------------------------------------
dir("rawData")


## ----load_gada, cache=FALSE----------------------------------------------
library(gada)


## ----setupPar------------------------------------------------------------
cnv.call <- setupParGADA(log2ratioCol = 4,
                        BAFcol = 5) 
cnv.call


## ----SBLcnv--------------------------------------------------------------
parSBL(cnv.call, estim.sigma2=TRUE, aAlpha=0.8)


## ----BackElimCNV---------------------------------------------------------
parBE(cnv.call, T=6, MinSegLen=25)


## ----summaryCNVseg-------------------------------------------------------
summ.cnvs <- summary(cnv.call)
summ.cnvs


## ----finNormalLimits-----------------------------------------------------
findNormalLimits(cnv.call)


## ----summaryCNVseg2------------------------------------------------------
summ2.cnvs <- summary(cnv.call, length.base=c(500,1e6))
summ2.cnvs


## ----GR_cnvs, cache=FALSE------------------------------------------------
library(GenomicRanges)
cnvs <- getCNVs(summ.cnvs)
cnvs


## ----subset_cnvs---------------------------------------------------------
subset(cnvs, sample=="CASE377")


## ----subset_cnvs2--------------------------------------------------------
subset(cnvs, seqnames=="chr6")


## ----load_Gviz, cache=FALSE----------------------------------------------
library(Gviz)
library(Homo.sapiens)


## ----cnvsPlot, fig.cap='CNV gains (dark gray) and losses (light gray) in the region chr8:45-60Mb of samples of \\textit{brgedata} package.', cache=FALSE----
rr <- GRanges("chr8:45.0e6-60.0e6")
plotCNVs(cnvs, range=rr, drawGenes = TRUE, 
         col.cnvs = c("darkgray", "lightgray"))


## ----casecont------------------------------------------------------------
casecont <- rep("control", length(cnvs))
casecont[grep("CASE", cnvs$sample)] <- "case"
cnvs$casecont <- as.factor(casecont)


## ----counts_casecont-----------------------------------------------------
cnvs.cases <- subset(cnvs, casecont=="case")
cnvs.controls <- subset(cnvs, casecont=="control")


## ----test_CNV------------------------------------------------------------
ss <- unique(cnvs.controls$sample)
cnvs.controls.list <- GRangesList(lapply(ss, function(x)
  subset(cnvs.controls, sample==x)))
count.controls <- countOverlaps(cnvs, cnvs.controls.list)

ss <- unique(cnvs.cases$sample)
cnvs.cases.list <- GRangesList(lapply(ss, function(x) 
  subset(cnvs.cases, sample==x)))
count.cases <- countOverlaps(cnvs, cnvs.cases.list)
cnvs$counts <- cbind(count.controls, count.cases)
cnvs[, c("counts")]


## ----overlap-------------------------------------------------------------
subsetByOverlaps(cnvs, cnvs[2,])


## ----numbercasescontrols-------------------------------------------------
agg <- aggregate(rep(1, length(cnvs)), 
                 by=list(casecont=cnvs$casecont, 
                         sample=cnvs$sample), 
                 FUN=sum)
agg
n <- table(agg$casecont)
n


## ----test_cnv------------------------------------------------------------
testCNV <- function(x, n) {
  tt <- matrix(c(x[1], n[1] - x[1], x[2], n[2] - x[2]), ncol=2)
  ans <- try(fisher.test(tt), TRUE)
  if (inherits(ans, "try-error"))
    out <- NA
  else
    out <- ans$p.value
  out
}  
  
cnvs$pvalue <- apply(cnvs$counts, 1, testCNV, n=n)
cnvs$BH <- p.adjust(cnvs$pvalue, method="BH")
cnvs[,c("counts", "pvalue", "BH")]


## ----load_CNV, cache=FALSE-----------------------------------------------
library(RTCGA.CNV)


## ----get_brca_cnv--------------------------------------------------------
head(BRCA.CNV)


## ----dlb_filter----------------------------------------------------------
dim(BRCA.CNV)
BRCA.CNV <- subset(BRCA.CNV, (Segment_Mean > log2(3/2) |
                     Segment_Mean < log2(1/2)) &
                     Num_Probes >= 5 & Num_Probes < 500)
dim(BRCA.CNV)


## ----createCNVobject-----------------------------------------------------
state <- ifelse(BRCA.CNV$Segment_Mean>0, 1, -1)

BRCA.CNV$Chromosome <- paste0("chr", BRCA.CNV$Chromosome)
BRCA.CNV$Chromosome[BRCA.CNV$Chromosome=="chr23"] <- "chrX"
brca.gr <- GRanges(seqnames = BRCA.CNV$Chromosome,
                   ranges = IRanges(start = BRCA.CNV$Start,
                                    end = BRCA.CNV$End),
                   LenProbe = BRCA.CNV$Num_Probes,
                   MeanAmp = BRCA.CNV$Segment_Mean, 
                   State = state,
                   sample = BRCA.CNV$Sample)

brca.gr$tumor <- rep(NA, length(brca.gr))
brca.gr$tumor[grep("10A", brca.gr$sample)] <- "normal"
brca.gr$tumor[grep("01A", brca.gr$sample)] <- "tumor"
brca.gr <- brca.gr[!is.na(brca.gr$tumor), ]

brca.gr$sample <- substr(brca.gr$sample, 1, 12)
brca.gr


## ----samplesBreast-------------------------------------------------------
length(unique(brca.gr$sample))


## ----plotBRCAcnv, fig.cap='CNVs gain (dark gray) and losses (light gray) in the region chr7:5.0-8.0Mb of samples belonging to the BRCA dataset from TCGA.', cache=FALSE----
rr <- GRanges("chr7:5e6-8e6")
plotCNVs(brca.gr, range = rr, drawGenes = TRUE,
         col.cnvs=c("darkgray", "lightgray"))


## ----analysis_brca_cnvs--------------------------------------------------
ans <- getCounts(brca.gr, "tumor")
brca.gr$counts <- ans$counts
n <- ans$n

brca.gr$pvalue <- apply(brca.gr$counts, 1, testCNV, n=n)
brca.gr$BH <- p.adjust(brca.gr$pvalue, method="BH")
brca.gr.sig <- subset(brca.gr[,c("counts", "pvalue", "BH")], 
                      brca.gr$BH < 0.01)
brca.gr.sig


## ----show_region---------------------------------------------------------
rr <- GRanges("chr16:653459-2025678")
subsetByOverlaps(brca.gr, rr)


## ----plotBRCAcnv2, fig.cap='CNVs in the region chr16:0.65-2.02Mb of one normal (dark gray) and tumor samples (light gray) beloning to BRCA dataset from TCGA.'----
rr <- GRanges("chr16:653459-2025678")
plotCNVs(brca.gr, range = rr, group="tumor", drawGenes = TRUE,
         col.group = c("darkgray", "lightgray"))


## ----select_gains_loses--------------------------------------------------
brca.gr.gains <- brca.gr[brca.gr$State==1,]
brca.gr.loses <- brca.gr[brca.gr$State==-1,]
length(brca.gr)
length(brca.gr.gains)
length(brca.gr.loses)


## ----load_dataMLPA-------------------------------------------------------
data(dataMLPA, package="CNVassocData")
head(dataMLPA)


## ----load_CNVassoc2, cache=FALSE-----------------------------------------
library(CNVassoc)


## ----MLPAsignal, fig.cap="Signal distributions for Gene 1 and Gene 2 from MLPA data example."----
par(mfrow=c(2,2),mar=c(3,4,3,1))
hist(dataMLPA$Gene1,main="Gene 1 signal histogram",xlab="",ylab="frequency")
hist(dataMLPA$Gene2,main="Gene 2 signal histogram",xlab="",ylab="frequency")
par(xaxs="i")
plot(density(dataMLPA$Gene1), main="Gene 1 signal density function",xlab="",
     ylab="density")
plot(density(dataMLPA$Gene2), main="Gene 2 signal density function",xlab="",
     ylab="density")


## ----plotSignalcasecon, fig.cap="Signal distribution for cases and controls of Gene 2 from MLPA data example."----
with(dataMLPA, plotSignal(Gene2, case.control=casco))


## ----plotSignalquanti, fig.cap="Signal distribution of Gene 2 from MLPA data example for a quantitative trait."----
with(dataMLPA, plotSignal(Gene2, case.control = quanti))


## ----calling_cnv1--------------------------------------------------------
CNV.1 <- cnv(x = dataMLPA$Gene1, threshold.0 = 0.06, 
             num.class = 3, mix.method = "mixdist")


## ----print_cnv1----------------------------------------------------------
CNV.1


## ----CNVcalling1, fig.cap="CNV calling of Gene 1 from the MLPA example.", fig.height=5, fig.width=6----
plot(CNV.1, case.control = dataMLPA$casco, main = "Gene 1")


## ----call_cnv2-----------------------------------------------------------
CNV.2 <- cnv(x = dataMLPA$Gene2, threshold.0 = 0.01, 
             mix.method = "mixdist")
CNV.2


## ----CNVcalling2, fig.cap="CNV calling of Gene 2 from the MLPA example.", fig.height=5, fig.width=6----
plot(CNV.2, case.control = dataMLPA$casco, main = "Gene 2")


## ----getProbs------------------------------------------------------------
probs.2 <- getProbs(CNV.2)
head(probs.2)


## ----assign_probs--------------------------------------------------------
CNV.2probs <- cnv(probs.2)


## ----quality_cnvtools----------------------------------------------------
CNVassoc::getQualityScore(CNV.1, type = "CNVtools")
CNVassoc::getQualityScore(CNV.2, type = "CNVtools")


## ----quality_class-------------------------------------------------------
CNVassoc::getQualityScore(CNV.1, type = "class")
CNVassoc::getQualityScore(CNV.2, type = "class")


## ----quality_canary------------------------------------------------------
CNVassoc::getQualityScore(CNV.1, type = "CANARY")
CNVassoc::getQualityScore(CNV.2, type = "CANARY")


## ----cnv_assoc_caco------------------------------------------------------
model1mul <- CNVassoc(casco ~ CNV.1, data = dataMLPA, model = "mul")
model2mul <- CNVassoc(casco ~ CNV.2, data = dataMLPA, model = "mul")


## ----print_assoc_mult----------------------------------------------------
model1mul


## ----summary_cnv_assoc---------------------------------------------------
summary(model1mul)
summary(model2mul)


## ----assoc_quanti--------------------------------------------------------
mod <- CNVassoc(quanti ~ CNV.2 + cov, family = "gaussian",
                data = dataMLPA, model = 'add', emsteps = 10)
coefficients(mod)


## ----summary_cnv_quanti--------------------------------------------------
summary(mod)


## ----global_cnv_test-----------------------------------------------------
CNVtest(model1mul, type = "Wald")
CNVtest(model1mul, type = "LRT")


## ----get_coef------------------------------------------------------------
coef(summary(model1mul, ref = 2))


## ----test_cnv_add--------------------------------------------------------
model2add <- CNVassoc(casco ~ CNV.2, data = dataMLPA, model = "add")
model2add


## ----summary_cnv_add-----------------------------------------------------
summary(model2add)


## ----anova_cnv-----------------------------------------------------------
anova(model2mul, model2add)


## ----load_NeveData-------------------------------------------------------
data(NeveData, package="CNVassocData")
intensities <- NeveData$data
pheno <- NeveData$pheno


## ----call_aCGH, eval=FALSE-----------------------------------------------
##   ######################################################
##   ### chunk number 1: Class of aCGH data
##   ######################################################
##   library(CGHcall)
##   Neve <- make_cghRaw(intensities)
## 
##   ######################################################
##   ### chunk number 2: Preprocessing
##   ######################################################
##   cghdata <- preprocess(Neve, maxmiss = 30, nchrom = 22)
## 
##   ######################################################
##   ### chunk number 3: Normalization
##   ######################################################
##   norm.cghdata <- normalize(cghdata, method = "median",
##                             smoothOutliers = TRUE)
## 
##   ######################################################
##   ### chunk number 4: Segmentation
##   ######################################################
##   seg.cghdata <- segmentData(norm.cghdata, method = "DNAcopy")
## 
##   ######################################################
##   ### chunk number 5: Calling
##   ######################################################
##   NeveCalled <- CGHcall(seg.cghdata, nclass = 3)
##   NeveCalled <- ExpandCGHcall(NeveCalled, seg.cghdata)


## ----called_load, echo=FALSE---------------------------------------------
data(NeveCalled, package="CNVassocData")


## ----getProbs_Neve-------------------------------------------------------
probs <- getProbs(NeveCalled)
probs[1:5, 1:7]


## ----cghregions, eval=FALSE----------------------------------------------
## library(CGHregions)
## NeveRegions <- CGHregions(NeveCalled)


## ----load_cghregions, echo=FALSE-----------------------------------------
data(NeveRegions, package="CNVassocData")


## ----get_prob_reg--------------------------------------------------------
probsRegions <- getProbsRegions(probs, NeveRegions, intensities)


## ----wg_assoc_cnv--------------------------------------------------------
pvalsCNV <- multiCNVassoc(probsRegions, formula = "pheno ~ CNV", 
                         model = "mult", num.copies = 0:2, 
                         cnv.tol = 0.01)


## ----get_sig_bh----------------------------------------------------------
pvalsBH <- getPvalBH(pvalsCNV)
head(pvalsBH)


## ----table6_gonz---------------------------------------------------------
cumsum(table(cut(pvalsBH[, 2], c(-Inf, 1e-5, 1e-4, 1e-3, 1e-2, 0.05))))


## ----install_mad, eval=FALSE---------------------------------------------
## devtools::install_github("isglobal-brge/MAD")


## ----load_mad, cache=FALSE-----------------------------------------------
library(mad)


## ----copy_data_no, eval=FALSE--------------------------------------------
## ss1 <- system.file("extdata/madData", package="brgedata")
## dir.create("rawData")
## ss2 <- "rawData"
## files <- list.files(ss1)
## file.copy(file.path(ss1,files), ss2)


## ----setup_mad-----------------------------------------------------------
mosaic <- setupParGADA.B.deviation(NumCols=6, GenoCol=6, 
                                    BAFcol=5, log2ratioCol=4)


## ----show_setup_mad------------------------------------------------------
mosaic


## ----parSBL_BE-----------------------------------------------------------
parSBL(mosaic, estim.sigma2=TRUE, aAlpha=0.8)
parBE.B.deviation(mosaic, T=9, MinSegLen=100)


## ----getMosaics----------------------------------------------------------
mosaic.gr <- getMosaics(mosaic)


## ----mosaics-------------------------------------------------------------
mosaic.gr


## ----plotMosaic, fig.cap="LRR and BAF corresponding to sample CONTROL191 from \\textit{brgedata} which was detected as mosaic loss.", cache=FALSE----
# argument 'col' is not necessary
mosaic <- addAnnot(mosaic) # add annotation
plotQMosaic(mosaic, sample="CONTROL191", chr=20, 
            regions=mosaic.gr,
            col.dots=c("black", "gray50"))


## ----plotMosaicZoom, fig.height=7, fig.width=4, fig.cap='LRR and BAF corresponding to sample CONTROL191 from \\textit{brgedata} detected as mosaic loss having the altered region delimited by vertical lines.', cache=FALSE----
plotZoomQMosaic(mosaic, sample="CONTROL191", chr=20,
                regions = mosaic.gr,
                col.dots=c("black", "gray50"))


## ----plotSegments, fig.dim=c(6, 4), fig.cap='Copy number alterations at genome level for different samples from \\textit{brgedata}. The legend corresponds to: UPD (1), deletion (2), duplication (3), trisomy (4).', cache=FALSE----
library(ggplot2)
plotSegments(mosaic.gr)


## ----plotMosaics, fig.cap='Mosaic alterations UPD (light gray) and LOH (dark gray) in the region chr8:63.5-67Mb of samples belonging to data available at the \\textit{brgedata} package.', cache=FALSE----
rr <- GRanges("chr16:63.5e6-67.0e6")
plotCNVs(mosaic.gr, range=rr, drawGenes = TRUE, mosaic=TRUE)


## ----install_madloy, eval=FALSE------------------------------------------
## devtools::install_github("isglobal-brge/MADloy")


## ----load_madloy, cache=FALSE--------------------------------------------
library(MADloy)


## ----copy_data2, results='hide'------------------------------------------
ss1 <- system.file("extdata/rawData", package="brgedata")
dir.create("LOY/rawData")
ss2 <- "LOY/rawData"
files <- list.files(ss1)
file.copy(file.path(ss1,files), ss2)


## ----get_data_path-------------------------------------------------------
rawDataPath <- "LOY/rawData"
files <- dir(rawDataPath)
length(files)
files[1:5]


## ----check_sex-----------------------------------------------------------
sex <- checkSex(rawDataPath, LRRCol=4)


## ----plotSex, fig.cap='Cluster of samples using LRR of the chromosome X and Y on individuals belonging to the \\textit{brgedata} data example.'----
plot(sex)


## ----print_sex-----------------------------------------------------------
sex


## ----sample_females------------------------------------------------------
sex$par$files[sex$class=="FEMALE"]


## ----remove_females------------------------------------------------------
files.males <- sex$par$files[sex$class!="FEMALE"]


## ----madloy--------------------------------------------------------------
mLRR <- madloy(files.males)
mLRR


## ----plotMLRRY, fig.cap = "Plot of MADloy object of males samples from \\textit{brgedata}."----
plot(mLRR, print.labels=TRUE, threshold=-0.3)


## ----loySample1, fig.cap="LRR (dots) in the mLRR-Y region (shaded) of SAMPLE\\_1 (normal sample) from \\textit{brgedata} data example. The solid horizontal line represents the median LRR values in the mLRR-Y region while the horizontal bottom dashed line is the expected value in a normal sample.", cache=FALSE----
plotIndLRR(mLRR, sample="SAMPLE_1")


## ----plotLOYLRR,  fig.cap="LRR in the mLRR-Y region (shaded) of SAMPLE\\_29 from the \\textit{brgedata} data example.", cache=FALSE----
plotIndLRR(mLRR, sample="SAMPLE_29")


## ----getLOY--------------------------------------------------------------
mLRR.call <- getLOY(mLRR)
mLRR.call

head(mLRR.call$res)


## ----LOY_Fosberg---------------------------------------------------------
table(mLRR.call$res$Fosberg)


## ----loy-----------------------------------------------------------------
loy <- checkBdev(mLRR.call)
loy


## ----get_discordant------------------------------------------------------
loyDat <- loy$data
sel <- which(loyDat$class == "other")
rownames(loyDat)[sel]


## ----discordant, fig.cap="LRR and BAF of a sample being tagged as 'other' by MADloy calling"----
plotIndSNPX(mLRR, sample="SAMPLE_98")


## ----plotLOY, fig.cap='LOY calling on the \\textit{brgedata} data example.'----
plot(mLRR.call, ylim=c(-1.5, 0.5), print.labels=TRUE)


## ----clean_cnvs, echo=FALSE, results='hide'------------------------------
#remove variables
rm(list=ls())

#garbage collection
gc(verbose = FALSE)

