## ----load_inveRsion, eval=FALSE, tidy=FALSE------------------------------
## library(VariantAnnotation)
## fl <-
##   "ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
## 
## #generate tabix file if needed
## #idx <- indexTabix(fl, "vcf")
## #tab <- TabixFile(fl, idx)
## 
## rng <- GRanges(seqnames="17", ranges=IRanges(start=43e6, end=47e6))
## vcf <- readVcf(fl, "hg19",param=rng)
## 


## ----snp_genotypes_inversion, eval=FALSE, tidy=FALSE---------------------
## #extract genotypes
## genos <- geno(vcf)$GT
## 
## popinfo <- read.delim("20130606_g1k.ped", as.is=TRUE)
## selCEU <- popinfo$Individual.ID[popinfo$Population=="CEU"]
## 
## #select CEUs
## genosCEU <-genos[, colnames(genos)%in%selCEU]
## selout <- rowMeans(genosCEU == "0|0") != 1
## genosCEU <- genosCEU[selout,]
## 
## #construct haplotype matrix
## hap1 <- matrix(0, ncol=ncol(genosCEU), nrow=nrow(genosCEU))
## hap1[genosCEU == "1|0"]<-1
## hap1[genosCEU == "1|1"]<-1
## 
## hap2<-matrix(0, ncol=ncol(genosCEU), nrow=nrow(genosCEU))
## hap2[genosCEU == "0|1"]<-1
## hap2[genosCEU == "1|1"]<-1
## 
## hap<-lapply(1:ncol(hap1), function(x) rbind(hap1[,x],hap2[,x]))
## hap<-do.call(rbind,hap)
## 
## colnames(hap)<-start(rowRanges(vcf))[selout]
## 
## #write matrix formatted for inveRsion
## write.table(hap, file="hap.txt",
##             col=TRUE, row=FALSE, quote=FALSE, sep="\t")
## 


## ----codeHap_inveRsion, message=FALSE, results='hide', tidy=FALSE--------
library(inveRsion)

hapCode<-codeHaplo(file="hap.txt", 
                   blockSize=5, minAllele=0.1,
                   saveRes=FALSE)


## ----scan_inveRsion, message=FALSE, results='hide', tidy=FALSE-----------
window <- 0.5
scanRes <- scanInv(hapCode,window=window,
                 saveRes=FALSE, geno=FALSE)

invList <- listInv(scanRes, hapCode=hapCode,
                   geno=FALSE,all=FALSE,thBic=0)

invList


## ----inversionScan, fig.cap="Scanning of the 17q21 region for signal of an inversion in the CEU population of the 1000 genomes.", message=FALSE----
plot(scanRes)


## ----load_hapmap, message=FALSE, tidy=FALSE------------------------------
library(snpStats)
fl <- "hapmap3_r2_b36_fwd.consensus.qc.poly"
genosHapMap <- read.plink(fl)

#select CEU
popinfo <- read.delim(file="relationships_w_pops_121708.txt", 
                      as.is=TRUE)
selpops <- popinfo$population%in%"CEU"
ceusIDs <- rownames(genosHapMap$genotypes)%in%popinfo$IID[selpops]

geno <- genosHapMap$genotypes[ceusIDs,]
annot <- genosHapMap$map[c("chromosome", "snp.name", "position")]



## ----invClust, message=FALSE,results='hide'------------------------------
library(invClust)

roi <- data.frame(chr=8, LBP=7897515, RBP=11787032, reg="8p23")

invCall <- invClust(roi=roi, wh = 1, geno=geno, annot=annot, dim=2)



## ----call_invClust, message=FALSE----------------------------------------
invCall
head(invCall["genotypes"])


## ----invCallCEU, fig.cap="Call of inv-8p23 in the CEU population of HapMap using invClust."----
plot(invCall)


## ----load_obesity_scoreInvHap--------------------------------------------
#load genotypes
library(snpStats)
genos.ob <- read.plink("obesity")

#recover PCA analysis
load(file="pca.RData")

eig <- data.frame(c1=pca$eigenvect[,1], c2=pca$eigenvect[,1])
rownames(eig)<-as.character(pca$sample.id)

#select individuals
selids<-rownames(genos.ob$genotypes)%in%as.character(pca$sample.id)

genos.ob$genotypes<-genos.ob$genotypes[selids,]

#confirm subject selection
identical(rownames(genos.ob$genotypes),as.character(pca$sample.id))


## ----inversionCall_scoreInvHap, tidy= FALSE------------------------------
library(scoreInvHap)

callinv <- scoreInvHap(genos.ob, 
                       SNPsR2 =  SNPsR2$inv8_001, 
                       hetRefs = hetRefs$inv8_001, 
                       Refs = Refs$inv8_001)

inv8p23 <- classification(callinv)

head(inv8p23)


## ----association_scoreInvHap---------------------------------------------
library(SNPassoc)

#load obesity phenotypes and merge data with 8p23 inversion call
ob.pheno <- read.delim(file="obesity.txt")
rownames(ob.pheno) <- ob.pheno$id

idspheno <- intersect(names(inv8p23), rownames(ob.pheno))

datainv <- data.frame(inv=inv8p23[idspheno], eig[idspheno,], 
                      ob.pheno[idspheno,])

#set up data for SNPassoc
ob.inv <- setupSNP(data=datainv, colSNPs=1, sep="")

#association test
association(obese ~ inv + c1+ c2+ country + smoke, ob.inv)



## ----clean_inversions, echo=FALSE, results='hide'------------------------
#remove variables
rm(list=ls())

#garbage collection
gc(verbose = FALSE)

