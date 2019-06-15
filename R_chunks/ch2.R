## ----load_LOAD, tidy=FALSE-----------------------------------------------
data <- read.delim(
file = "phs000168.v2.pht000707.v2.p2.c1.LOAD610K_Subject_Phenotypes.GRU.txt",
comment.char = "#")

names(data)[1:18]


## ----extract_raceandPCA--------------------------------------------------
race <- as.factor(data$Race)
table(race)


## ----PCA, fig.cap="Genome-wide PCA of LOAD-NIA dataset.", cache=FALSE----
mycols <- c("gray90", "black", "gray70", "gray50", "gray20", "white")
cols <- as.character(factor(race, labels=mycols))
plot(data$AllEthnicity_PC1, data$AllEthnicity_PC2,
     type="n", main="PCA LOAD-NIA", 
     xlab="PC1", ylab="PC2")
points(data$AllEthnicity_PC1, data$AllEthnicity_PC2,
       col = cols, pch=0:5)

legend("bottomright", c("white", "black", "american indian",
                        "asian", "other", "missing"),
       col=mycols, pch=0:5)


## ----load_ResultsLOAD----------------------------------------------------
results <- read.table(
           file="Run_precompute_1877_unrelated_samples.assoc", 
           header=TRUE)
head(results)


## ----QQ, fig.cap="QQ plot LOAD-NIA.", message=FALSE----------------------
library(snpStats)
qq.chisq(results$CHISQ)


## ----pval----------------------------------------------------------------
pval <- results[, c("SNP", "CHR", "BP", "OR", "P" )]
selSig <- p.adjust(pval$P, method="bonf")
subset(pval, selSig<0.05)


## ----get_R_functions_book, echo=FALSE, cache=FALSE-----------------------
dd <- "../R"
ff <- dir(dd)
for (i in ff)
  source(file.path(dd,i))


## ----manhattanLOAD, fig.cap="GWAS results on unrelated individuals of the LOAD-NIA study. Each point is a SNP in the genomic location described in the x axis. On the left y axis, the -$log_{10}(P)$ is plotted. The dashed line indicates significance threshold that corrects from multiple testing (i.e. genome-wide significance level).", fig.height=7, fig.width=9----
library(tidyverse)
library(ggplot2)
library(ggrepel)

# remove missing P-values
pval.nona <- pval[!is.na(pval), ]
# select SNPs with p-value lower than 0.01
# to speed up the plot without loosing information
pval.nonasig <- pval.nona[pval.nona$P<1e-2, ]

manhattanPlot(pval.nonasig, color=c("gray90", "gray40"))


## ----select_region, echo=FALSE, eval=FALSE-------------------------------
## sel <- pval[,"CHR"]==19 & pval[,"BP"]>=50021054-100000 &
##        pval[,"BP"]<=50106291+100000 & !is.na(pval[,"P"])
## 
## write.table(pval[which(sel),],
##             file="sigP.txt", quote=FALSE,
##             row=FALSE, col=TRUE)


## ----genoquery-----------------------------------------------------------
library(GEOquery)
gsm.expr <- getGEO("GSE63061", destdir = ".")
gsm.expr <- gsm.expr[[1]]


## ----genoquery2----------------------------------------------------------
show(gsm.expr)


## ----getdata-------------------------------------------------------------

#get transcriptomic data and apply log2
expr <- log2(exprs(gsm.expr))
dim(expr)

#get phenotype data
pheno <- pData(phenoData(gsm.expr))
status <- pheno$characteristics_ch1
status <- gsub("status: ","", as.character(status))
table(status)

#create case control variable from AD and CLT labels
selcaco <- status%in%c("AD","CTL")
caco <- rep(NA, length(status))
caco[status=="AD"] <- 1
caco[status=="CTL"] <- 0


## ----analysisGEO, message=FALSE------------------------------------------
#locate rows in expression data for transcripts in APP, TOMM40 and APOC1
genesIDs <- fData(gsm.expr)[13]
genesIDs <- as.character(unlist(genesIDs))

selAPP <- which(genesIDs%in%c("APP"))[2:3]
selTOMM40 <- which(genesIDs%in%c("TOMM40"))
selAPOC1 <- which(genesIDs%in%c("APOC1"))

selTranscripts <- c(selAPP, selTOMM40, selAPOC1)
labTranscripts <- c("APP", "APP", "TOMM40", "APOC")


## ----ExpressionGEO, fig.cap="Transcript differences in APP, TOMM40 and APOC between Alzheimer's cases and controls", message=FALSE----
library(vioplot)
par(mfrow=c(2,2), mar=c(2,4,2,1))

for(trans in 1:4){
  x <- selTranscripts[trans]
  exprsel <- expr[x,]
  
  #association
  mod <- glm(caco ~ exprsel, family="binomial")
  pval <- summary(mod)$coeff[2,4]
  
  #plot
  lab <- paste0(labTranscripts[trans], "\n P = ", 
              as.character(round(pval,3)))
  vioplot(exprsel[which(caco==0)], exprsel[which(caco==1)], 
          col="gray80", names=c("controls", "cases"))
  title(lab, cex.main=0.8)
  }


## ----GetDatamethylationGEO-----------------------------------------------
gsm.meth <- getGEO("GSE80970", destdir = ".")
gsm.meth <- gsm.meth[[1]]


## ----AnalDatamethylationGEO----------------------------------------------
met <- log2(exprs(gsm.meth))

#get phenotype data
pheno <- pData(phenoData(gsm.meth))
status <- pheno$characteristics_ch1.1
status <- gsub("disease status: ","", as.character(status))
table(status)

#create case control variable from Alzheimer's disease and control labels
caco <- as.numeric(status=="Alzheimer's disease")

#locate CpGs in APP
genesIDs <- fData(gsm.meth)$UCSC_RefGene_Name
genesIDs <- as.character(unlist(genesIDs))

selAPP <- grep("APP;APP",genesIDs)
length(selAPP)


## ----PvalmethylationGEO--------------------------------------------------
pval <- sapply(selAPP, function(probes){
  mod <- glm(caco ~ met[probes,], family="binomial")
  summary(mod)$coeff[2,4]
})
names(pval) <- rownames(met[selAPP,])
sigPval <- which(pval < 0.05/length(pval))
names(sigPval)


## ----methylationGEO, fig.cap="Methylation differences near APP's promoter between Alzheimer's cases and controls.", message=FALSE, fig.dim = c(4, 4)----

#significant CpG info
cg <- names(sigPval)
P <- as.character(round(min(pval),3))

#plot
metapp <- met[cg,]
lab <- paste0(cg, "\n 2.5Kb from APP's promoter\n P = ", P)
vioplot(metapp[which(caco==0)], metapp[which(caco==1)], 
          col="gray90", names=c("controls", "cases"))
title(lab)


## ----PvalmethylationGEOAdjust--------------------------------------------
#retrieve covariates
tissue <- pheno$characteristics_ch1
sex <- pheno$characteristics_ch1.3
age <- pheno$characteristics_ch1.4
age <- as.numeric(substr(as.character(age),11,12))

#model of signficant finding adjusted by covariates
probesig <- selAPP[sigPval]
mod <- glm(caco ~ met[probesig,]+tissue+sex+age, family="binomial")
round(summary(mod)$coefficients,4)


## ----loadRecount---------------------------------------------------------
library("recount")


## ----GTEXdownload, eval=FALSE--------------------------------------------
## #download and load data
## download_study('SRP012682', outdir = ".")


## ----GTEXdata, eval=FALSE------------------------------------------------
## load("rse_gene.Rdata")
## #scale RNA-seq data
## rse_gene <- scale_counts(rse_gene)
## rse_gene


## ----GTEXdata1, message=FALSE, eval=FALSE--------------------------------
## #obtain expression data and annotation
## recountCounts <- assays(rse_gene)$counts
## recountMap <- rowRanges(rse_gene)
## 
## #locate data for APP
## indGene <- which(unlist(recountMap$symbol)=="APP")
## selAPP <- grep(names(indGene),rownames(recountCounts))
## recountCounts2APP <- log2(recountCounts[selAPP,]+1)


## ----GTEXdata2, eval=FALSE-----------------------------------------------
## #select data for four different brain tissues
## gtexPd <- colData(rse_gene)
## 
## mask1 <- which(gtexPd$smtsd=="Brain - Cerebellum")
## recountCrb <- recountCounts2APP[mask1]
## 
## mask2 <- which(gtexPd$smtsd=="Brain - Cortex")
## recountCtx <- recountCounts2APP[mask2]
## 
## mask3 <- which(gtexPd$smtsd=="Brain - Frontal Cortex (BA9)")
## recountFctx <- recountCounts2APP[mask3]
## 
## mask4 <- which(gtexPd$smtsd=="Brain - Hippocampus")
## recountHipp <- recountCounts2APP[mask4]


## ----GTEXload4plot, echo=FALSE-------------------------------------------
#save(recountCrb,recountCtx,recountFctx, recountHipp, file="recount.RData"
load("recount.RData")


## ----APPexpressionPLot, fig.cap="Expression of APP from RNA-seq count data across four brain tissues (CRB: Cerebellum, CRT: Cortex, FRT: Frontal, HIP: Hippocampus) as measured in the GTEx project.", message=FALSE, tidy=FALSE, fig.dim = c(5, 5)----
library(vioplot)
vioplot(recountCrb, recountCtx, recountFctx, recountHipp,
        col="gray90", 
        names=c("CRB", "CRT", "FRT", "HIP"))
title("Scaled count data (RNA-seq) for APP")


## ----RTCGA.rnaseq--------------------------------------------------------
library(RTCGA.rnaseq)
library(Biobase)

BRCA.rnaseq_ExpressionSet <- convertTCGA(BRCA.rnaseq)  
genenms <- rownames(BRCA.rnaseq_ExpressionSet)

#locate the full gene names
genenms[grep("XBP1",genenms)]
genenms[grep("ESR1",genenms)]

#get expression data
expr <- expressionsTCGA(BRCA.rnaseq, 
                        extract.cols = c("XBP1|7494","ESR1|2099")) 

#get ID data
barcode <- expr$bcr_patient_barcode

#locate tumor data encoded with “01” string
cancertissue <- substr(barcode, 14, 15)=="01"

#locate healthy tissue data encoded with “11” string

controltissue <- substr(barcode, 14, 15)=="11"

#extract ID string relative to subjects
Subidscancer<-substr(barcode[cancertissue],1,12)
Subidscontrol<-substr(barcode[controltissue],1,12)

#find which subjects have both cancer and healthy tissue data
commonids<-intersect(Subidscancer, Subidscontrol)
selectids<-unlist(sapply(commonids, grep, barcode))

#get expression data for subjects and genes
expr<-expr[selectids,]

#transfrom
exprLogXBP1<-log2(expr$`XBP1|7494`+1)
exprLogESR1<-log2(expr$`ESR1|2099`+1)


## ----TCGAscater, fig.cap="Correlation between RNA-seq expression levels of XBP1 and ESR1 in tumorous breast tissue.", message=FALSE----
plot(exprLogXBP1, exprLogESR1, xlab="log-expression XBP1",
     ylab="log-expression XBP1", pch=21, bg="gray90")


## ----TCGassocDataforAssoc------------------------------------------------
#outcome variable
selbarcode <- barcode[selectids]
status <- rep(NA, length(selectids))

#locate tumor 01 and healty 11 data in selected subjects
selcancertissue <- substr(selbarcode, 14, 15)=="01"
selcontroltissue <- substr(selbarcode, 14, 15)=="11"

status[selcancertissue] <- 1
status[selcontroltissue] <- 0


## ----TCGassocAssocXBP1---------------------------------------------------
#associations
mod1 <- glm(status ~ exprLogXBP1, family="binomial")
res1 <- summary(mod1)$coefficients["exprLogXBP1",]
res1


## ----TCGassocAssocESR1---------------------------------------------------
mod2 <- glm(status ~ exprLogESR1, family="binomial")
res2 <- summary(mod2)$coefficients["exprLogESR1",]
res2


## ----TCGAviolin, fig.cap="RNA-seq expression in XBP1 and ESR1 in tumorous breast tissue (tum) and healthy tissue (cont).", message=FALSE----
library(vioplot)

#violin plot
##labels
lab1 <- paste0("XBP1 expression \n P = ", 
              as.character(format(res1["Pr(>|z|)"], sci=TRUE,dig=3)))

lab2 <- paste0("ESR1 expression \n P = ", 
             as.character(round(res2["Pr(>|z|)"],3)))

#remove missing values
notna<-!is.na(exprLogXBP1) & !is.na(exprLogESR1) & !is.na(status)

par(mfrow=c(1,2), mar=c(5,3,4,1))
vioplot(exprLogXBP1[which(status==0 & notna)],
        exprLogXBP1[which(status==1 & notna)],
        col="gray90", names=c("cont", "tum"))
title(lab1)

vioplot(exprLogESR1[which(status==0 & notna)],
        exprLogESR1[which(status==1 & notna)],
        col="gray90", names=c("cont", "tum"))
title(lab2)



## ----NHANESlibrary, eval=FALSE-------------------------------------------
## install.packages("RNHANES")


## ----NHANESfiles, warning= FALSE, message = FALSE------------------------
library(RNHANES)
files <- nhanes_data_files()
#variables <- nhanes_variables()
nhanes_search(files, "environmental phenols")


## ----NHANESdataPAH_G-----------------------------------------------------
nhanes_dat <- nhanes_load_data("PAH_G", "2011-2012", demographics = TRUE)
rownames(nhanes_dat)<-as.character(nhanes_dat$SEQN)

dim(nhanes_dat)

nhanes_dat[1:10, 1:10]


## ----NHANEScorrelation---------------------------------------------------
library(survey)
library(jtools)

des <- nhanes_survey_design(nhanes_dat)

svycor(~ log(URXP01) + log(URXP02) + log(URXP03) + log(URXP04) + 
         log(URXP05) , design = des, na.rm = TRUE)


## ----NHANESdataCANCER----------------------------------------------------
medical <- nhanes_load_data("MCQ", "2011-2012", demographics = TRUE)
rownames(medical) <- as.character(medical$SEQN)


## ----NHANESdataCANCERMerge-----------------------------------------------
#common names
nms <- intersect(rownames(medical), rownames(nhanes_dat))

#locate exposures that start with URX
wh <- grep("URX", colnames(nhanes_dat))

#define PHA matrix exposure
PAH <- log(nhanes_dat[nms,wh])
vars <- colnames(PAH)

#select medical variables 
cancer<-medical[nms,c("MCQ220","RIAGENDR","RIDAGEYR")]

#select complete cases
dat <- data.frame(PAH, cancer)
selrow <- complete.cases(dat)
dat<-dat[selrow,]

#Recode cancer=1, control=0
dat$MCQ220 <- 2-dat$MCQ220


## ----NHANESdataCANCERAdjustedage-----------------------------------------
ExWAS<-sapply(vars,  function(pha){
  fit <- summary(glm(MCQ220 ~ dat[,pha] + RIAGENDR + RIDAGEYR, 
                     data = dat, family = "binomial"))
  
  fit$coef[2,c(1,2,4)]
})
   
t(ExWAS)  


## ----clean_ch2, echo=FALSE, results='hide'-------------------------------
#remove variables
rm(list=ls())

#garbage collection
gc(verbose = FALSE)


