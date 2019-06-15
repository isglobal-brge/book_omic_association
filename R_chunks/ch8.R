## ----load_brge_2, cache=FALSE--------------------------------------------
library(brgedata)
ls("package:brgedata")


## ----get_expo------------------------------------------------------------
brge_expo


## ----load_rexposome, cache=FALSE-----------------------------------------
library(rexposome)


## ----table_missing_phenotype---------------------------------------------
tableMissings(brge_expo, set = "phenotypes")


## ----plotMissingExpo, fig.cap='Missing data pattern of the \\textit{brgedata} exposome dataset.'----
library(ggplot2)
plotMissings(brge_expo, set = "exposures") +
    ggtitle("Missing Values in INMA-Sabadell Exposome")


## ----imputation----------------------------------------------------------
brge_expo_complete <- imputation(brge_expo)
tableMissings(brge_expo_complete, set="exposures")


## ----plotAll, fig.cap='Distribution of exposures accross families on brge exposome dataset.'----
plotFamily(brge_expo, family = "all")


## ----plotAllSex, fig.cap='Distribution of PCBs stratified by sex the \\textit{brgedata} exposome dataset.'----
plotFamily(brge_expo, family = "PCBs", group = "Sex")


## ----standar-------------------------------------------------------------
brge_expo_standard <- standardize(brge_expo)
brge_expo_standard


## ----plotStandar, fig.cap='Distribution of PCBs estratified by sex on the \\textit{brgedata} exposome dataset.', echo=FALSE----
plotFamily(brge_expo_standard, family = "PCBs", group = "Sex")


## ----plotPCAexpo, fig.cap='Principal component analysis on the \\textit{brgedata} exposome standardized data', fig.width=3.5, fig.height=3.5----
brge_expo_pca <- pca(brge_expo_standard)
rexposome::plotPCA(brge_expo_pca, set = "all")


## ----plotPCAsex, fig.cap='Principal component analysis on the \\textit{brgedata} exposome standardized data using gender as an illustrative variable.'----
rexposome::plotPCA(brge_expo_pca, 
                   set = "samples", 
                   phenotype = "Sex")


## ----plotCor, fig.cap='Correlation matrix between the \\textit{brgedata} exposure datasets.', cache=FALSE----
exp_cr <- correlation(brge_expo_standard, 
                      use = "pairwise.complete.obs", 
                      method.cor = "pearson")
plotCorrelation(exp_cr, type = "matrix")


## ----exwas---------------------------------------------------------------
library(splines)
resExwas <- exwas(brge_expo_standard, 
                  formula = Asthma ~ Sex + Age, 
                  family = "binomial")
resExwas


## ----manhattanExwas, fig.cap='Manhattan plot of exposures obtained from ExWAS.'----
clr <- rainbow(length(familyNames(brge_expo_standard)))
names(clr) <- familyNames(brge_expo_standard)
plotExwas(resExwas, color = clr) +  
  ggtitle("ExWAS - Univariate Approach")


## ----effectExwas, fig.cap='Effects (beta values) and confidence intervals of exposures corresponding to the asthma EXWAS example.'----
plotEffect(resExwas)


## ----expo_trans, cache=FALSE---------------------------------------------
library(omicRexposome)
library(MultiDataSet)
mds <- createMultiDataSet()
mds <- add_exp(mds, brge_expo_standard)
mds <- add_genexp(mds, brge_gexp)
resGene <- association(mds, formula = ~ Sex + Age, 
                       sva = "fast",
                       expset = "exposures", 
                       omicset = "expression")
resGene


## ----hits----------------------------------------------------------------
tableHits(resGene, th = 1e-06)


## ----lambda--------------------------------------------------------------
tableLambda(resGene, trim = 0.95)


## ----plotVolcanoNO2, fig.cap='Volcano of NO2 and transcriptomic data from \\textit{brgedata} exposome data example.'----
plotAssociation(resGene, rid = "NO2_p", type = "volcano",
                tPV = -log10(1e-04), tFC = 0.1, 
                show.effect=FALSE) +
  ggplot2::ggtitle("NO2_p")


## ----plotQQNO2, fig.cap='QQ plot of NO2 and transcriptomic data from \\textit{brgedata} exposome data example.'----
plotAssociation(resGene, rid = "NO2_p", type = "qq", 
                show.lambda=FALSE) +  ggplot2::ggtitle("NO2_p")


## ----cluster, cache=FALSE------------------------------------------------
brge_expo_complete <- imputation(brge_expo)
brge_expo_standard2 <- standardize(brge_expo_complete)

brge_expo_cluster <- clustering(brge_expo_standard2,
                                method = Mclust, 
                                G=2)
table(rexposome::classification(brge_expo_cluster))


## ----plotCluster, fig.cap='Exposome clustering using the mclust method including exposures from the \\textit{brgedata} example.', fig.dim=c(5,7), cache=FALSE----
plotClassification(brge_expo_cluster)


## ----group_vs_Gene, cache=FALSE------------------------------------------
mds2 <- createMultiDataSet()
mds2 <- add_cls(mds2, brge_expo_cluster)
mds2 <- add_genexp(mds2, brge_gexp)
resGene2 <- association(mds2, formula = ~ Sex + Age, 
                        sva = "fast",
                        expset = "cluster", 
                        omicset = "expression")
tableHits(resGene2, th = 1e-06)


## ----plotQQVolcano2, fig.cap='Volcano and QQ plots of transcriptomic data analysis considering exposome clusters as the factor variable.'----
par(mfrow=c(2,1))
plotAssociation(resGene2, type = "volcano",
                tPV = -log10(1e-06), tFC = 1, 
                show.effect=FALSE) +
  ggplot2::ggtitle("Transcriptome - Cluster group 1 vs 2")

plotAssociation(resGene2, type = "qq", show.lambda=FALSE) +
  ggplot2::ggtitle("NO2_p")


## ----remove_partial, echo=FALSE, results='hide', cache=FALSE-------------
#remove variables
rm(resGene, resGene2, brge_expo_cluster, resExwas, 
   exp_cr, brge_expo_complete)

#garbage collection
gc(verbose = FALSE)


## ----expoMethy, eval=FALSE, cache=FALSE----------------------------------
## mds <- add_methy(mds, brge_methy)
## resMethy <- association(mds, formula = ~ Sex + Age,
##                         sva = "fast",
##                         expset = "exposures",
##                         omicset = "methylation")
## mds2 <- add_methy(mds2, brge_methy)
## resMethy2 <- association(mds2, formula = ~ Sex + Age,
##                          sva = "fast",
##                          expset = "cluster",
##                          omicset = "methylation")


## ----save_methy, echo=FALSE, eval=FALSE----------------------------------
## save(resMethy, resMethy2, file="resMethy.RData")


## ----load_methy, echo=FALSE, cache=FALSE---------------------------------
load("resMethy.Rdata")


## ----plotMethy, fig.cap='Volcano plot of NO2 analysis on epigenomic data.'----
plotAssociation(resMethy, rid = "NO2_p", type = "volcano",
                tPV = -log10(1e-04), tFC = 0.1, 
                show.effect=FALSE) +
  ggplot2::ggtitle("NO2_p")


## ----plotMethyClust, fig.cap='Volcano plot of epigenomic data comparing exposome clusters.'----
plotAssociation(resMethy2, type = "volcano",
                tPV = -log10(1e-08), tFC = 0.3, 
                show.effect=FALSE) +
  ggplot2::ggtitle("Epigenome - Exposome cluster group 1 vs 2")


## ----clean_ch8, echo=FALSE, results='hide', cache=FALSE------------------
#remove variables
rm(list=ls())

#garbage collection
gc(verbose = FALSE)

