## ----genoquerych5, eval=FALSE--------------------------------------------
## library(GEOquery)
## gsm.expr <- getGEO("GSE63061", destdir = ".")
## gsm.expr <- gsm.expr[[1]]
## show(gsm.expr)


## ----genoqueryRDATAch5, echo=FALSE---------------------------------------
load("GSE63061.Rdata")
gsm.expr <- gsm.expr[[1]]
show(gsm.expr)


## ----Trancriptomedatach5-------------------------------------------------
expr <- exprs(gsm.expr)
dim(expr)


## ----phenotypedatach5----------------------------------------------------
#get phenotype data
pheno <- pData(phenoData(gsm.expr))
status <- pheno$characteristics_ch1
status <- gsub("status: ","", as.character(status))
table(status)


## ----phenofactorch5------------------------------------------------------
fstatus <- factor(status)
levels(fstatus) <- c("case","case", "cont", "cont", "case", "case","case")
table(fstatus)


## ----modelsofSVAch5------------------------------------------------------
library(sva)
mod0 <- model.matrix( ~ 1, data = fstatus)
mod <- model.matrix( ~ fstatus)


## ----nsvSVAch5-----------------------------------------------------------
n.sv <- num.sv(expr, mod, method="leek")
n.sv


## ----svach5SVAch5--------------------------------------------------------
svobj <- sva(expr, mod, mod0, n.sv=2)
names(svobj)


## ----SVA, fig.cap="SVA components of GSE63061. Cases are black dots and controls are grey.", tidy=FALSE----
col <- fstatus
levels(col) <- c("black","grey")
plot(svobj$sv[,1:2],col=as.character(col), xlab="sva 1", 
     ylab="sva 2", pch=16)


## ----svach5notNormalized-------------------------------------------------
exprNN <- read.table("GSE63061_non-normalized.txt",
                     header=TRUE, check.names=FALSE, as.is=TRUE)
probesID <- exprNN[,1]
selprobes <- probesID%in%rownames(expr)

exprNN <- exprNN[selprobes,as.character(pheno$title)]
exprNN <- log2(exprNN)


## ----nsvNNormalizedch5---------------------------------------------------
n.sv <- num.sv(exprNN, mod, method="leek")
n.sv


## ----svFitNNormalizedch5-------------------------------------------------
exprNNmat <- as.matrix(exprNN)
svobjN <- sva(exprNNmat, mod, mod0, n.sv=n.sv)


## ----SVAUnorm, fig.cap="SVA components of unnormalized data for GSE63061. Cases are black dots and controls are grey.", tidy=FALSE----
plot(svobjN$sv[,1:2],col=as.character(col), xlab="sva 1", 
     ylab="sva 2", pch=16)


## ----assocNotNomrSVA-----------------------------------------------------
ADstatus <- factor(status)
levels(ADstatus)
levels(ADstatus) <- c("1","0", "0", "0", "0", "0","1")

summary(glm(ADstatus ~ svobjN$sv, family="binomial" ))


## ----assocSVAstatus------------------------------------------------------
summary(glm(svobj$sv[,1]  ~ status))


## ----expbrainch5---------------------------------------------------------
expr <- read.table(file="expression_data_Brain.txt", header=TRUE, 
                   as.is=TRUE )

probesIDS<-expr[,1]
expr <- as.matrix(log2(expr[,-1]))

expr[1:5,1:5]


## ----expbrainch5missings-------------------------------------------------
selna<-!is.na(rowSums(expr))
expr<-expr[selna,]


## ----covbrainch5---------------------------------------------------------
cov <- read.table(file="covariate_data_Brain.txt", header=TRUE, 
                  as.is=TRUE )
names(cov)
table(cov$institute_sample_source)


## ----covmodelbrainch5----------------------------------------------------
modcombat <- model.matrix(~1, data=cov)

expr_corrected <- ComBat(dat=expr, batch=cov$institute_sample_source, 
                         mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

expr_corrected[1:5,1:5]


## ----nsvabrainfbrainch5--------------------------------------------------
fbrain <- as.factor(cov$brain_region)
table(fbrain)

mod0 <- model.matrix( ~ 1, data = fbrain)
mod <- model.matrix( ~ fbrain)


## ----exprsvabrainch5-----------------------------------------------------
n.sv <- num.sv(expr, mod, method="leek")
n.sv
svobj <- sva(expr, mod, mod0, n.sv=n.sv)


## ----svacorrectedbrainch5------------------------------------------------
n.sv_corrected <- num.sv(expr_corrected, mod, method="leek")
n.sv
svobj_corrected <- sva(expr_corrected, mod, mod0, n.sv=2)


## ----assocsvabrainbatchch5-----------------------------------------------
batch <-factor(cov$institute_sample_source)
summary(aov(lm(svobj$sv[,2] ~ batch)))


## ----assocsvabraincorrectedch5-------------------------------------------
summary(aov(lm(svobj_corrected$sv[,2] ~ batch)))


## ----plotassocsvabrainch5, fig.cap="SVA components of corrected and uncorrected for batch effects the expression data of GSE8919. SVAs of corrected data are more uniform across batches despite differences between them being still present.", tidy=FALSE----
par(mfrow=c(2,2))
boxplot(svobj$sv[,1] ~ batch, ylab="sva 1", ylim=c(-0.2,0.2), 
        xlab="batch", main="Uncorrected")
lines(c(0,20),c(0,0), lty=2)

boxplot(svobj_corrected$sv[,1] ~ batch, ylab="sva 1", ylim=c(-0.2,0.2), 
        xlab="batch", main="Corrected")
lines(c(0,20),c(0,0), lty=2)

boxplot(svobj$sv[,2] ~ batch, ylab="sva 2", ylim=c(-0.2,0.2), 
        xlab="batch",main="Uncorrected")
lines(c(0,20),c(0,0), lty=2)

boxplot(svobj_corrected$sv[,2] ~ batch, ylab="sva 2", ylim=c(-0.2,0.2), 
        xlab="batch", main="Corrected")
lines(c(0,20),c(0,0), lty=2)


## ----clean_ch5, echo=FALSE, results='hide'-------------------------------
#remove variables
rm(list=ls())

#garbage collection
gc(verbose = FALSE)


