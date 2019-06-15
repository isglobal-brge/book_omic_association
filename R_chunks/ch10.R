## ----load_data-----------------------------------------------------------
data(breastMulti, package="brgedata")


## ----targetAZ------------------------------------------------------------
library(Biobase)
library(GenomicRanges)
targetRange <- GRanges("chr19:350000-45550000")
breastMulti.target <- breastMulti[, , targetRange]


## ----dims2---------------------------------------------------------------
dims(breastMulti.target)


## ----corr_meth_exp-------------------------------------------------------
library(MEAL)
methExprs <- correlationMethExprs(breastMulti.target)
head(methExprs)


## ----corrMethTrans-------------------------------------------------------
methy <- breastMulti[["methylation"]]
resER <- runPipeline(methy, variable_names = "ER.Status",
                     sva=TRUE)
dmps <- getProbeResults(resER)
head(dmps)

cpgs <- rownames(dmps)[1:10] #select top-10 CpGs
methExprs <- correlationMethExprs(breastMulti,
                                  sel_cpgs=cpgs)
methExprs


## ----massive, fig.cap="Significant correlation between the top CpG, associated with ER status, and gene transcription from breast TCGA data."----
library(vioplot)
par(mfrow=c(1,2))
cpg <- assay(methy)["cg19996355",]
boxplot(cpg[methy$ER.Status=="Positive"],
        cpg[methy$ER.Status=="Negative"],
        col="gray80", names=c("Positive", "Negative"),
        ylab="cg19996355", las=2)

geneExpr <- breastMulti[["expression"]]
gene <- assay(geneExpr)["PBX4",]
plot(cpg, gene, xlab="cg19996355",
     ylab="PBX4", type="n")
points(cpg, gene, pch=21, bg="gray90")
abline(lm(gene ~ cpg), lwd=2)


## ----load_lusc-----------------------------------------------------------
data(lusc, package="brgedata")


## ----cit-----------------------------------------------------------------
library(cit)
L <- lusc[,"LOY"]               # instrumental variable
G <- lusc[, "DDX3Y"]            # mediator
T <- lusc[, "Cancer"]           # outcome
C <- lusc[, "age"]              # confounders  
cit.bp(L, G, T, C, rseed=1234)


## ----mod_step1-----------------------------------------------------------
mod <- glm(Cancer ~ LOY + age , data=lusc, family="binomial")
summary(mod)


## ----mod_step2-----------------------------------------------------------
mod.M <- glm(DDX3Y ~ LOY + age, data=lusc)
summary(mod.M)


## ----mod_step3-----------------------------------------------------------
mod.Y <- glm(Cancer ~ LOY + DDX3Y + age, data=lusc, family="binomial")
summary(mod.Y)


## ----mediation-----------------------------------------------------------
library(mediation)
set.seed(12345)
res <- mediate(mod.M, mod.Y, treat='LOY', mediator='DDX3Y')
summary(res)


## ----load_haploR, cache=FALSE--------------------------------------------
library(haploR)


## ----haploR--------------------------------------------------------------
x <- queryHaploreg(query=c("rs6931936", "rs12201995"))
x


## ----haploR2-------------------------------------------------------------
sel <- as.numeric(x$r2) > 0.9
x2 <- x[sel, c("pos_hg38", "rsID", "Motifs", "gwas",
               "GENCODE_name", "Enhancer_histone_marks")]
x2


## ----regulome------------------------------------------------------------
xx <- queryRegulome(c("rs6931936", "rs12201995"))
xx$res.table


## ----haploR_3------------------------------------------------------------
x3 <- x[, c("eQTL", "rsID")]
x3


## ----PCAex, fig.cap='Hypothetical representation of gene expression of two genes on 45 indivuals.', echo=FALSE----
set.seed(12345)
y <- rnorm(45, 4, 1)
x <- 2 * y + rnorm(45)
plot(x,y, cex.lab=1.2, cex.axis = 1.15, xlab="GATA3", ylab="XPBD1")
points(x,y, pch=16, cex=1.15)
abline(lm(y~x), col="gray50", lwd=2)
arrows(6.5, 2.5, 6.5, 3.2, length=0.15, col="gray50")
text(6.5, 2.3, "First PCA", col="gray50")
segments(10, 3.4, 7.6, 5.0, col="gray50")
arrows(11.7, 3.6, 10, 3.6, length=0.15, col="gray50")
text(12.0, 3.6, "Second \n PCA", adj=0, col="gray50")


## ----load_ade4-----------------------------------------------------------
library(ade4)
geneExpr <- breastMulti[["expression"]]
dim(geneExpr)


## ----pcaRNAseq-----------------------------------------------------------
breastPCA <- dudi.pca(assay(geneExpr),
                      scannf=FALSE, nf=2)
summary(breastPCA)


## ----pcaBreast, fig.cap='Gene expression principal component analysis of the BRCA dataset from the TCGA project.'----
library(made4)
group <- geneExpr$ER.Status
out <- ord(assay(geneExpr), type="pca", classvec=group)

par(mfrow=c(2,1))
cols <- c("black", "gray50")
mycols <- ifelse(group=="Positive", cols[1], cols[2])
plotarrays(out$ord$co, classvec=group, arraycol = cols)
plotgenes(out, col="gray50")


## ----list_genes----------------------------------------------------------
ax1 <- topgenes(out, axis=1, n=5, ends="pos")
ax2 <- topgenes(out, axis=2, n=5, ends="neg")
cbind(pos=ax1, neg=ax2)


## ----sPCA, results='hide'------------------------------------------------
library(PMA)
X <- t(assay(geneExpr))
spca <- SPC(X, sumabsv=3, K=2)


## ----sPCAcv--------------------------------------------------------------
subabsvs.grid <- seq(10, sqrt(ncol(X)),len=20)
spca.cv <- SPC.cv(X, sumabsvs=subabsvs.grid)
spca.cv


## ----sPCA2, results='hide'-----------------------------------------------
spca <- SPC(X, sumabsv=spca.cv$bestsumabsv,
            K=2)


## ----sPCA_out------------------------------------------------------------
rownames(spca$u) <- rownames(X)
rownames(spca$v) <- colnames(X)
head(spca$u)
head(spca$v)


## ----SPCAplot, fig.cap='Gene expression sparse principal component analysis of the BRCA dataset from the TCGA project.'----
plot(spca$u, type="n", xlab="First sPCA",
     ylab="Second sPCA")
cols <- c("black", "gray50")
mycols <- ifelse(group=="Positive", cols[1], cols[2])
points(spca$u, pch=16, cex=1.3, col=mycols)
legend("bottomright", c("Positive", "Negative"),
       col=cols, pch=16)
abline(h=0, lty=2, lwd=2)
abline(v=0, lty=2, lwd=2)
grid()


## ----sPCA_genes----------------------------------------------------------
ss <- spca$v[,1]
ss.nonzero <- ss[ss!=0]
ss.pos <- ss[order(ss, decreasing = TRUE)][1:5]
ss.pos


## ----scca, results='hide'------------------------------------------------
library(PMA)
df1 <- t(assay(breastMulti[["expression"]]))
df2 <- t(exprs(breastMulti[["RPPA"]]))

ddlist <- list(df1, df2)
perm.out <- MultiCCA.permute(ddlist,
                             type=c("standard", "standard"),
                             trace=FALSE)
resMultiCCA <- MultiCCA(ddlist, ncomponents=1,
                        penalty=perm.out$bestpenalties,
                        ws=perm.out$ws.init,
                        type=c("standard", "standard"),
                        trace=FALSE)


## ----multiCCA_out--------------------------------------------------------
rownames(resMultiCCA$ws[[1]]) <- colnames(df1)
rownames(resMultiCCA$ws[[2]]) <- colnames(df2)
head(resMultiCCA$ws[[1]])
head(resMultiCCA$ws[[2]])


## ----nonzeroSCCA1--------------------------------------------------------
genes <- resMultiCCA$ws[[1]]
genes.sig <- genes[genes[,1]!=0,]
head(genes.sig)
length(genes.sig)


## ----nonzeroSCCA2--------------------------------------------------------
proteins <- resMultiCCA$ws[[2]]
proteins.sig <- proteins[proteins[,1]!=0,]
head(proteins.sig)
length(proteins.sig)


## ----getcia--------------------------------------------------------------
library(made4)
library(omicade4)
resCIA <- cia(assay(breastMulti[["expression"]]), 
              exprs(breastMulti[["RPPA"]]))


## ----ciaBreast, fig.cap='Gene expression and protein coinertia analysis of the BRCA dataset from the TCGA project.'----
mycols <- c("black", "gray50")
cols <- ifelse(group=="Positive", 
               mycols[1], mycols[2])
plot(resCIA, classvec=group, nlab=3, clab=0, 
     cpoint=3, col=cols)


## ----top_features_pos----------------------------------------------------
topVar(resCIA, axis=1, topN=5, end="positive")


## ----top_features_neg----------------------------------------------------
topVar(resCIA, axis=1, topN=5, end="negative")


## ----mcia----------------------------------------------------------------
ll <- as.list(breastMulti)
names(ll)
resMCIA <- mcia(ll[ c("methylation", "expression", 
                      "miRNA", "RPPA") ] )


## ----mciaBreast, fig.cap='Multiple coinertia analysis of the BRCA dataset from the TCGA project.'----
mycols <- c("black", "gray50")
cols <- ifelse(group=="Positive", mycols[1], mycols[2])
plot(resMCIA, axes=1:2, sample.lab=FALSE, 
     sample.legend=FALSE, phenovec=group, 
     gene.nlab=2, sample.col=cols, df.pch=2:5,
     df.color=c("black", "gray30", "gray50", "gray90"))


## ----top_features_m------------------------------------------------------
topVar(resMCIA, end="positive", axis=1, topN=5)


## ----rgcca_pca, fig.show='hide', out.extra=''----------------------------
library(RGCCA)

X <- t(assay(breastMulti[["expression"]]))

# one omic X 
# Design matrix C
# Shrinkage parameters tau = c(tau1, tau2)
pca.with.rgcca <- rgcca(A = list(X, X),
                        C = matrix(c(0, 1, 1, 0), 2, 2),
                        tau = c(1, 1),
                        ncomp = c(2, 2))


## ----get_R_functions_book_2, echo=FALSE----------------------------------
dd <- "../R"
ff <- dir(dd)
for (i in ff)
  source(file.path(dd,i))


## ----plotPCArgcca, fig.cap='Principal component analysis using RGCCA on the BRCA dataset from the TCGA project.'----
plotInd(pca.with.rgcca, group=group,
        col.list=c("gray50", "black"))


## ----rgcca_cca, fig.show='hide'------------------------------------------
X <- t(assay(breastMulti[["expression"]]))
Y <- t(exprs(breastMulti[["miRNA"]]))
cca.with.rgcca <- rgcca(A= list(X, Y), 
                            C = matrix(c(0, 1, 1, 0), 2, 2),
                            tau = c(0, 0),
                            ncomp=c(2,2))


## ----gcca_rgcca, eval = FALSE--------------------------------------------
## # X1 = omic1, ..., XJ = omicJ, X_{J+1} = [X1, ..., XJ]
## # (J+1)*(J+1) Design matrix C
## C = matrix(c(0, 0, 0, ..., 0, 1,
##              0, 0, 0, ..., 0, 1,
##              ...,
##              1, 1, 1, ..., 1, 0), J+1, J+1)
## #Shrinkage parameters tau = c(tau1, ...,  tauJ, tau_{J+1})
## gcca.with.rgcca = rgcca(A= list(X1, ..., XJ, cbind(X1, ..., XJ)),
##                         C = C, tau = rep(0, J+1),
##                         scheme = "factorial")


## ----mcia_rgcca, eval = FALSE--------------------------------------------
## # X1 = omic1, ..., XJ = omicJ, X_{J+1} = [X1, ..., XJ]
## # (J+1)*(J+1) Design matrix C
## C = matrix(c(0, 0, 0, ..., 0, 1,
##              0, 0, 0, ..., 0, 1,
##              ...,
##              1, 1, 1, ..., 1, 0), J+1, J+1)
## # Shrinkage parameters tau = c(tau1, ...,  tauJ, tau_{J+1})
## mcia.with.rgcca = rgcca(A= list(X1, ..., XJ, cbind(X1, ..., XJ)),
##                         C = C, tau = c(rep(1, J), 0),
##                         scheme = "factorial")


## ----rgcca_1-------------------------------------------------------------
X <- t(assay(breastMulti[["expression"]]))
Y <- t(exprs(breastMulti[["miRNA"]]))
Z <- t(exprs(breastMulti[["RPPA"]]))
A <- list(X,Y,Z)
A <- lapply(A, function(x) scale2(x, bias = TRUE))


## ----C-------------------------------------------------------------------
C <- matrix(c(1,0,1,0,1,1,0,0,1), nrow=3, byrow = TRUE)
C


## ----rgcca.factorial, fig.show='hide'------------------------------------
rgcca.factorial <- rgcca(A, C=C, tau = rep(0, 3),
                         scheme ="factorial", 
                         ncomp=c(2,2,2), scale = FALSE,
                         verbose = FALSE)


## ----sol-----------------------------------------------------------------
head(rgcca.factorial$a[[1]]) #for RNAseq
head(rgcca.factorial$a[[2]]) #for miRNA
head(rgcca.factorial$a[[3]]) #for proteins


## ----plotRGCCAfac, fig.cap='Multiomic data analysis using RGCCA on the BRCA dataset from the TCGA project.'----
plotInd(rgcca.factorial, group=group,
        col.list=c("gray50", "black"))


## ----topVarsRight--------------------------------------------------------
topVars(rgcca.factorial, axis=1, end="pos", topN=5)


## ----topVarsLeft---------------------------------------------------------
topVars(rgcca.factorial, axis=1, end="neg", topN=5)


## ----ave-----------------------------------------------------------------
rgcca.factorial$AVE


## ----clean_ch10, echo=FALSE, results='hide'------------------------------
#remove variables
rm(list=ls())

#garbage collection
gc(verbose = FALSE)

