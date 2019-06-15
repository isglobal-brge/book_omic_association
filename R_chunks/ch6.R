## ----Unormalizedch6------------------------------------------------------
exprMat<- read.table("GSE63061_non-normalized.txt", 
                     header=TRUE, check.names=FALSE, 
                     as.is=TRUE)

probesID<-exprMat[,1]
exprMat<-exprMat[,-1]
subsID<-colnames(exprMat)
exprMat<-as.matrix(exprMat)


## ----Nomalizationch6-----------------------------------------------------
library(limma)

exprNormdat<-normalizeBetweenArrays(exprMat, method="quantile")

exprNormdat[1:5,1:5]


## ----contNomalizationch6-------------------------------------------------
L2exprMat<-log2(exprMat)
L2exprNorm<-log2(exprNormdat)

comp<-cbind(L2exprMat,L2exprNorm)
nsubs<-ncol(L2exprMat)

gr<-factor(rep(c("Not Normalized","Normalized"),each=nsubs))


## ----plotnormch6,fig.cap="Distribution of probe intensities for each subject in the expression data setGSE63061. In grey, the not normalized distribution and in black the normalized ones."----
plotDensities(comp, group = gr, col=c("black","grey"), legend="topright")



## ----NormalizedExpSetch6-------------------------------------------------
library(Biobase)

phenoDF <- data.frame(subsID, row.names=subsID)
phenodat <- AnnotatedDataFrame(data=phenoDF)

featureDF <- data.frame(probesID, row.names=probesID)
featureData <- AnnotatedDataFrame(data=featureDF)

exprNorm <- ExpressionSet(assayData=L2exprNorm, 
                          phenoData=phenodat, 
                          featureData=featureData)

show(exprNorm)


## ----genoquerych6, eval=FALSE--------------------------------------------
## library(GEOquery)
## gsm.expr <- getGEO("GSE63061", destdir = ".")
## gsm.expr <- gsm.expr[[1]]
## show(gsm.expr)


## ----genoqueryRDATAch6, echo=FALSE---------------------------------------
load("GSE63061.Rdata")
gsm.expr <- gsm.expr[[1]]
show(gsm.expr)


## ----probesIDsch6--------------------------------------------------------
probesIDsGEO<-as.character(fData(gsm.expr)$ID)
probesID<-as.character(fData(exprNorm)$probesID)


## ----selprobesch6--------------------------------------------------------
selprobes<-probesID%in%probesIDsGEO
nselprobes<-sum(selprobes)

exprNorm.sel<-exprNorm[selprobes,]


## ----comprenomalizations,fig.cap="Distribution of probe intensities for each subject in the expression data set GSE63061. In grey, the normalized distributions as reported in GEO and in black the normalized distributions as obtained here."----
exprGEO<-exprs(gsm.expr)
expr<-exprs(exprNorm.sel)

comp<-cbind(exprGEO,expr)

nsubs<-ncol(exprGEO)
gr<-factor(rep(c("Normalized GEO","Normalized"),each=nsubs))

plotDensities(comp, group = gr, col=c("black", "grey"), legend="topright")


## ----comprenomalizationssclaed,fig.cap="Distribution of probe intensities for each subject in the expression data set GSE63061 displaced to zero mean. In grey, the normalized distributions as reported in GEO and in black the normalized distributions as obtained here."----
exprScaledGEO <- scale(exprGEO, scale=FALSE)
exprScaled <- scale(expr, scale=FALSE)

comp<-cbind(exprScaledGEO,exprScaled)
plotDensities(comp, group = gr, col=c("black", "grey"), legend="topright")


## ----filterplot,fig.cap="Histogram for the difference of hybridization means between filtered-in probes and filtered-out probes across subjects."----
exprNorm.notsel<-exprNorm[!selprobes,]
expr.notsel<-exprs(exprNorm.notsel)

hist(colMeans(expr.notsel)-colMeans(expr), 
     main="Differences between means \n (filtered - not
     filtered probes)", xlab="Differece between means")


## ----filterprobes--------------------------------------------------------
expr.all<-exprs(exprNorm)

qnthreshold<-quantile(expr.all,0.5)
checkquants<-sapply(as.data.frame(t(expr.all)), quantile,0.8)

tb<-table(selprobes, checkquants>qnthreshold)
tb
sum(diag(tb))/sum(tb)


## ----associationch6------------------------------------------------------
#expression data
exprGEO<-exprs(gsm.expr)

#get phenotype data
pheno <- pData(phenoData(gsm.expr))


## ----defineStatuschr6----------------------------------------------------
status <- pheno$characteristics_ch1
status <- gsub("status: ","", as.character(status))
fstatus <- factor(status)
levels(fstatus)<-gsub(" ", "", levels(fstatus))

table(fstatus)


## ----sexandage-----------------------------------------------------------
age <- substr(pheno$characteristics_ch1.2, 6,7)
age<-as.numeric(age)
sex <- pheno$characteristics_ch1.3


## ----sva4assocch6--------------------------------------------------------
library(sva)
phenodat<-data.frame(fstatus, age, sex)

mod0 <- model.matrix( ~ age+sex, data = phenodat)
mod <- model.matrix( ~ fstatus+age+sex, data = phenodat)

svobj <- sva(exprGEO, mod, mod0, n.sv=2)


## ----desingmatassocch6---------------------------------------------------
sv1 <- svobj$sv[,1]
sv2 <- svobj$sv[,1]
design <- model.matrix(~ 0 + fstatus + sex + age + sv1 + sv2)
colnames(design) <- c(levels(fstatus),"age","sex", "sva1","sva2")


## ----fitassocch6---------------------------------------------------------
fit <- lmFit(exprGEO, design)

contrast.matrix <- makeContrasts(AD-CTL, MCI-CTL, AD-MCI, levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)


## ----vennDiag, fig.cap="Venn diagram indicating the number of genes that are differentially expressed for each comparison in the Alzheimer's disease example.", fig.dim=c(4,4)----
results <- decideTests(fit2)
vennDiagram(results, cex=c(1, 0.7, 0.5))


## ----getCoef-------------------------------------------------------------
colnames(fit2$coef)


## ----contrast1-----------------------------------------------------------
tt<-topTable(fit2, coef=1, adjust="BH")

genesIDs <- as.character(fData(gsm.expr)$ILMN_Gene)
names(genesIDs)<-rownames(gsm.expr)

data.frame(genesIDsgenes=genesIDs[rownames(tt)], 
           logFC=tt$logFC, pvalAdj=tt$adj.P.Val)


## ----volcanoch6, fig.cap="Volcano plot showing the P-values of association against the log2 fold expression change."----
volcanoplot(fit2, coef = 1, highlight=5, names=genesIDs, cex=0.2)


## ----microVSrnaseq, echo=FALSE, fig.cap="Gene expression distribution obtained from RNA-seq and microarray data from a hypothetical gene."----
  
dd <- rnbinom(10000, mu=2, size=100)
dd2 <- rnorm(10000, mean=15, sd=8)

par(mfrow=c(1,2))
plot(prop.table(table(dd)), type="h", lwd=4,  
     main="RNAseq data", xlab="counts",  
     ylab="Probability")

plot(density(dd2), main="Microarray data", xlab="intensity", 
     lwd=3, col="gray90")



## ----get_hapmap----------------------------------------------------------
library(tweeDEseqCountData)
data(pickrell)
pickrell.eset


## ----read_and_normalize_breast-------------------------------------------
counts <- exprs(pickrell.eset)
lib.size <- colSums(counts)
NormByTotalNrReads <- sweep(counts, 2, FUN="/", lib.size)


## ----get_annotation------------------------------------------------------
data(annotEnsembl63)
head(annotEnsembl63)
genes.ok <- intersect(as.character(rownames(counts)), 
                      as.character(rownames(annotEnsembl63)))
geneAnnot <- annotEnsembl63[genes.ok,]
counts.ok <- counts[genes.ok,]
identical(rownames(geneAnnot), rownames(counts.ok))


## ----norm_rpkm-----------------------------------------------------------
width <- geneAnnot$Length
NormByRPKM <- t(t(counts.ok / width *1000)
                    /colSums(counts.ok)*1e6)


## ----tmm-----------------------------------------------------------------
library(tweeDEseq)
NormByTMM <- normalizeCounts(counts.ok, method="TMM")


## ----cqn-----------------------------------------------------------------
library(cqn)       
annotation <- geneAnnot[,c("Length", "GCcontent")]
NormByCQN <- normalizeCounts(counts.ok, method="cqn",
                             annot=annotation)


## ----checkNorm, fig.cap = "Comparison of log-ratio count intensity of samples 1 and 2 from the Pickrell dataset."----
MbyT <- log2(NormByTotalNrReads[, 1] / NormByTotalNrReads[, 2])
MbyRPKM <- log2(NormByRPKM[, 1] / NormByRPKM[, 2])

par(mfrow=c(1,2))
hist(MbyT, xlab="log2-ratio", main="Total reads")
abline(v=0, col="red")
hist(MbyRPKM, xlab="log2-ratio", main="RPKM") 
abline(v=0, col="red")


## ----MAraw, fig.cap='MA-plot of raw (e.g. non-normalized) Pickrell data on samples 1 and 3.', cache=FALSE----
library(edgeR)
maPlot(counts[,1], counts[,2], pch=19, cex=.5, ylim=c(-8,8),
       allCol="darkgray", lowess=TRUE,  
       xlab=expression( A == log[2] (sqrt(Sample1 %.% Sample3))  ),
       ylab=expression(M == log[2](Sample1/Sample3)))
       grid(col="black")


## ----load_edgeR, echo=FALSE, cache=FALSE---------------------------------
library(edgeR)


## ----MAall, fig.cap="MA-plot of Pickrell data on samples 1 and 3 for raw data and normalized data using RPKM, TMM and CQN methods.", echo=FALSE, cache=FALSE----
par(mfrow=c(2,2))
maPlot(counts[,1], counts[,2], pch=19, cex=.5, ylim=c(-8,8),
       allCol="darkgray", lowess=TRUE,  
       xlab=expression( A == log[2] (sqrt(Sample1 %.% Sample3))  ),
       ylab=expression(M == log[2](Sample1/Sample3)))
       grid(col="black")
title("Raw data") 
              
NormByRPKM <- NormByRPKM[complete.cases(NormByRPKM),]       
maPlot(NormByRPKM[,1], NormByRPKM[,2], pch=19, cex=.5, ylim=c(-8,8),
       allCol="darkgray", lowess=TRUE,  
       xlab=expression( A == log[2] (sqrt(Sample1 %.% Sample3))  ),
       ylab=expression(M == log[2](Sample1/Sample3)))
       grid(col="black")       
title("RPKM normalization")       

maPlot(NormByTMM[,1], NormByTMM[,2], pch=19, cex=.5, ylim=c(-8,8),
       allCol="darkgray", lowess=TRUE,  
       xlab=expression( A == log[2] (sqrt(Sample1 %.% Sample3))  ),
       ylab=expression(M == log[2](Sample1/Sample3)))
       grid(col="black")  
title("TMM normalization")              

maPlot(NormByTMM[,1], NormByTMM[,2], pch=19, cex=.5, ylim=c(-8,8),
       allCol="darkgray", lowess=TRUE,  
       xlab=expression( A == log[2] (sqrt(Sample1 %.% Sample3))  ),
       ylab=expression(M == log[2](Sample1/Sample3)))
       grid(col="black")        
title("CQN normalization")       


## ----filtering-----------------------------------------------------------
NormByCQN.f <- filterCounts(NormByCQN)
dim(NormByCQN)
dim(NormByCQN.f)


## ----edegR_common_dispersion---------------------------------------------
library(edgeR)
sex <- pData(pickrell.eset)$gender
d <- DGEList(counts = NormByCQN, group = as.factor(sex))
d <- estimateCommonDisp(d)
d$common.dispersion


## ----show_dispersion-----------------------------------------------------
sqrt(d$common.dispersion)


## ----edgeR_tagwise-------------------------------------------------------
d <- estimateTagwiseDisp(d)
names(d)


## ----dseq_object---------------------------------------------------------
library(DESeq)
cds <- newCountDataSet(NormByCQN, as.factor(sex))
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds


## ----edgeR_common--------------------------------------------------------
resEdgeR.common <- exactTest(d, pair=c("female", "male"), 
                              dispersion="common")
resEdgeR.tagwise <- exactTest(d, pair=c("female", "male"), 
                              dispersion="tagwise")


## ----resEdgeR_common-----------------------------------------------------
topTags(resEdgeR.common)


## ----resEdgeR_tag--------------------------------------------------------
topTags(resEdgeR.tagwise)


## ----DESeq---------------------------------------------------------------
resDESeq <- nbinomTest(cds, "female", "male")


## ----sig_resEdgeR--------------------------------------------------------
out <- resDESeq[resDESeq$padj < 0.05,]
out.o <- out[order(out$padj),]
head(out.o)


## ----DESeq2--------------------------------------------------------------
library(DESeq2)
counts <- exprs(pickrell.eset)
pheno <- pData(pickrell.eset)
dds <- DESeqDataSetFromMatrix(countData = counts,
           colData = pheno, 
           design = ~ gender)


## ----ref-----------------------------------------------------------------
table(pheno$gender)


## ----filtering2----------------------------------------------------------
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]


## ----DESeq2_results------------------------------------------------------
dds <- DESeq(dds)


## ----visualize_results---------------------------------------------------
resDESeq2  <- results(dds, pAdjustMethod = "fdr")
head(resDESeq2)


## ----visualize_results2--------------------------------------------------
resDESeq2 <- resDESeq2[order(resDESeq2$padj), ]
resDESeq2


## ----maplotDESEq2, fig.cap='MA-plot comparing both conditions of Pickrell data.'----
plotMA(resDESeq2)


## ----plotGeneTop, fig.cap='Gene expression levels of females and males corresponding to the gene most associated with gender status in the Pickrell dataset.'----
plotCounts(dds, gene = "ENSG00000129824",
           intgroup = "gender")


## ----voom, fig.cap='Mean-variance relationsip corresponding to the Pickrell dataset.'----
library(limma)
design <- model.matrix( ~ gender, data=pheno)
v <- voom(counts, design=design, plot=TRUE)


## ----voomLimma, fig.cap='Mean-variance relationship of voom transformed genes from Pickrell dataset.'----
fit <- lmFit(v, design)
fit <- eBayes(fit)
topTable(fit, coef=ncol(design))


## ----clean_ch6, echo=FALSE, results='hide'-------------------------------
#remove variables
rm(list=ls())

#garbage collection
gc(verbose = FALSE)

