## ----get_R_functions_book2, echo=FALSE, cache=FALSE----------------------
library(ggplot2)
dd <- "../R"
ff <- dir(dd)
for (i in ff)
  source(file.path(dd,i))


## ----load_asthma---------------------------------------------------------
data(asthma, package = "brgedata")
str(asthma, list.len=9)
asthma[1:5, 1:10]


## ----install_git, eval=FALSE---------------------------------------------
## require(devtools)
## devtools::install_github("isglobal-brge/SNPassoc")


## ----prepare data--------------------------------------------------------
library(SNPassoc)
asthma.s <- setupSNP(data=asthma, colSNPs=7:ncol(asthma), sep="")


## ----class SNP-----------------------------------------------------------
head(asthma.s$rs1422993)
class(asthma.s$rs1422993)


## ----summary SNP---------------------------------------------------------
summary(asthma.s$rs1422993)


## ----descriptive---------------------------------------------------------
ss <- summary(asthma.s, print=FALSE)
head(ss)


## ----plotMissing, fig.cap="Missing genotypes (black) of asthma data example.", fig.width=5.5, fig.height=5----
plotMissing(asthma.s, print.labels.SNPs = FALSE)


## ----HWE-----------------------------------------------------------------
hwe <- tableHWE(asthma.s)
head(hwe)


## ----HWE controls--------------------------------------------------------
hwe2 <- tableHWE(asthma.s, casecontrol)

#SNPs is HWE in the whole sample but not controls
snpNHWE <- hwe2[,1]>0.05 & hwe2[,2]<0.05
rownames(hwe2)[snpNHWE]
hwe2[snpNHWE,]



## ----remove_SNPs---------------------------------------------------------
snps.ok <- rownames(hwe2)[hwe2[,2]>=0.001]
pos <- which(colnames(asthma)%in%snps.ok, useNames = FALSE)
asthma.s <- setupSNP(asthma, pos, sep="")


## ----association---------------------------------------------------------
association(casecontrol ~ rs1422993, data = asthma.s)


## ----maxrs1422993--------------------------------------------------------
maxstat(asthma.s$casecontrol, asthma.s$rs1422993)


## ----mode inheritance----------------------------------------------------
association(casecontrol ~ rs1422993, asthma.s, model="dominant")


## ----adjusted------------------------------------------------------------
association(casecontrol ~ rs1422993 + country + smoke, asthma.s)


## ----stratified----------------------------------------------------------
association(casecontrol ~ rs1422993 + survival::strata(gender), asthma.s)


## ----subset--------------------------------------------------------------
association(casecontrol ~ rs1422993, asthma.s, 
                subset=country=="Spain")


## ----quantitative--------------------------------------------------------
association(bmi ~ rs1422993, asthma.s) 


## ----more than 1 SNP-----------------------------------------------------
ans <- WGassociation(casecontrol, data=asthma.s)
head(ans)


## ----more than 1 SNP adjusted, eval=FALSE--------------------------------
## ans.adj <- WGassociation(casecontrol ~ country + smoke, asthma.s)
## head(ans.adj)


## ----more than 1 SNP scan, eval=FALSE------------------------------------
## ans.fast <- scanWGassociation(casecontrol, asthma.s)


## ----plotGWAS, fig.cap="Manhattan-type plots for different genetic models to assess the association between case-control status and SNPs in the asthma example.", fig.width=6, fig.height=7.1----
plot(ans)


## ----max-statistic-------------------------------------------------------
ans.max <- maxstat(asthma.s, casecontrol)
head(ans.max)


## ----max-rs1422993-------------------------------------------------------
#Maximum-statistics for rs1422993
ans.max[,"rs1422993"]


## ----max-Bonferroni------------------------------------------------------
#minimum P-value across SNPs
min(ans.max["Pr(>z)",])


## ----OR_several_SNPs, results='hide'-------------------------------------
infoTable <- WGstats(ans)
infoTable$rs1422993


## ----getNiceTable, results='hide'----------------------------------------
library(xtable)
out <- getNiceTable(ans[c("rs1422993", "rs184448")])

nlines <- attr(out, "nlines")
hlines <- c(-1, -1, 0, cumsum(nlines+1), nrow(out), nrow(out))

print(xtable(out, caption='Genetic association using
                different genetic models from asthma 
                data example of rs1422993 and rs184448 
                SNPs obtained with SNPassoc.',
             label = 'tab-2SNPs'),
      tabular.enviroment="longtable", file="tableSNPs",
      floating=FALSE,  include.rownames = FALSE, 
      hline.after= hlines, sanitize.text.function=identity)


## ----inter_nat2, echo=FALSE, fig.cap="Interaction plot illustrating the different effect of NAT2 genotypes in never and ever smokers.", fig.width=4, fig.height=4----
ever <- c(1, 1.04, 1.24)
never <- c(1, 1.01, 0.96)
nat2 <- c("GG", "AG", "AA")
plot(1:3, ever, type="n", xlab="NAT2 genotypes", 
     ylab="Risk of bladder cancer (OR)", ylim=c(0.7, 1.4),
     axes=FALSE, cex.lab=1.2)
axis(2, cex.axis=1.15)
axis(1, at=1:3, labels=nat2, cex.axis=1.15) 
points(1:3, ever, pch=16, col="black", cex=1.5)
lines(1:3, ever, col="black", lty=2, lwd=2)
points(1:3, never, pch=16, col="gray90", cex=1.5)
lines(1:3, never, col="gray90", lty=2, lwd=2)
abline(h=1, lty=3)
segments(2, 0.96, 2, 1.12, col="black")
segments(2, 0.94, 2, 1.07, col="gray90")
segments(3, 1.16, 3, 1.32, col="black")
segments(3, 0.86, 3, 1.08, col="gray90")
legend("bottomleft", c("Ever smokers", "Never smokers"), pch=16, col=c("black", "gray90"), cex=1.1)


## ----snpxsmoke-----------------------------------------------------------
association(casecontrol ~ dominant(rs1422993)*factor(smoke), 
            data=asthma.s)


## ----snpxsnp-------------------------------------------------------------
association(casecontrol ~ rs1422993*factor(rs184448), 
            data=asthma.s, model.interaction = "dominant" )


## ----load_libraries_hap--------------------------------------------------
library(LDheatmap)
library(genetics)


## ----laodBiomart---------------------------------------------------------
library(biomaRt)
mart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
nrow(listFilters(mart))
head(listFilters(mart))
listFilters(mart)[11,]


## ----get_annot-----------------------------------------------------------
snps <- labels(asthma.s)
snpInfo <- getBM(c("refsnp_id", "chr_name", "chrom_start", "allele"),  
                  filters = c("snp_filter"), 
                  values = snps, mart = mart)
head(snpInfo)


## ----sel_chr7------------------------------------------------------------
mask <- with(snpInfo, chr_name=="7" & chrom_start>34.5e6 & 
               chrom_start<35.0e6)
snps.sel <- snpInfo[mask, "refsnp_id"]
sel <- which(names(asthma)%in%snps.sel)
asthma.hap <- setupSNP(asthma, sel, sep="")


## ----info_hap------------------------------------------------------------
snp.pos <- snpInfo[mask, "chrom_start"]
snp.geno <- data.frame(lapply(asthma.hap[, snps.sel], genotype))


## ----LDplot, fig.cap="Linkage disequilibrium heatmap of selected SNPs located in chr7:34500000-35000000 region belonging to the asthma dataset.", fig.width=5, fig.height=5----
LDplot <- LDheatmap(snp.geno, LDmeasure = "r",
 title = "Pairwise LD in r^2", add.map = TRUE,
 color = grey.colors(20), name = "myLDgrob", add.key = TRUE,
 flip=TRUE, SNP.name=snps.sel)


## ----haplo_em------------------------------------------------------------
library(haplo.stats)
snpsH <- c("rs714588", "rs1023555",  "rs898070")
genoH <- make.geno(asthma.hap, snpsH)
em <- haplo.em(genoH, locus.label = snpsH, miss.val = c(0, NA))
em


## ----hap_assoc-----------------------------------------------------------
trait <- asthma.hap$casecontrol
mod <- haplo.glm(trait ~ genoH,           
                 family="binomial", 
                 locus.label=snpsH,
                 allele.lev=attributes(genoH)$unique.alleles,
                 control = haplo.glm.control(haplo.freq.min=0.05))   
intervals(mod)


## ----sliding_window_hap--------------------------------------------------
snpsH2 <- colnames(snp.geno)[6:15]
genoH2 <- make.geno(asthma.hap, snpsH2)
haplo.score <- list()
for (i in 4:7) {
 trait <- asthma.hap$casecontrol
 haplo.score[[i-3]] <- haplo.score.slide(trait, genoH2, 
                          trait.type="binomial",
                          n.slide=i,
                          simulate=TRUE,
                          sim.control=score.sim.control(min.sim=100,
                                       max.sim=200)) 
 }


## ----plotSliding, fig.cap="Sliding window approach varying haplotype size from 4 up to 7 SNPs located in chr7:34500000-35000000 region from asthma data example"----
par(mfrow=c(2,2))
for (i in 4:7) {
    plot(haplo.score[[i-3]])
    title(paste("Sliding Window=", i, sep=""))
 }


## ----hap_assoc_bestH-----------------------------------------------------
snpsH3 <- snpsH2[4:7]
genoH3 <- make.geno(asthma.hap, snpsH3)
mod <- haplo.glm(trait~genoH3,           
                 family="binomial", 
                 locus.label=snpsH3,
                 allele.lev=attributes(genoH3)$unique.alleles,
                 control = haplo.glm.control(haplo.freq.min=0.05))      
intervals(mod)


## ----lrt_hap-------------------------------------------------------------
lrt <- mod$lrt
pchisq(lrt$lrt, lrt$df, lower=FALSE)


## ----lrt_hap_Adj---------------------------------------------------------
smoke <- asthma.hap$smoke
mod.adj.ref <- glm(trait ~ smoke, family="binomial")
mod.adj <- haplo.glm(trait ~ genoH3 + smoke ,           
                 family="binomial", 
                 locus.label=snpsH3,
                 allele.lev=attributes(genoH3)$unique.alleles,
                 control = haplo.glm.control(haplo.freq.min=0.05))

lrt.adj <- mod.adj.ref$deviance - mod.adj$deviance
pchisq(lrt.adj, mod.adj$lrt$df, lower=FALSE)


## ----score_load----------------------------------------------------------
library(PredictABEL)
data(asthma, package = "brgedata")
dd.s <- setupSNP(asthma, 7:ncol(asthma), sep="")
ans <- WGassociation(casecontrol, dd.s, model="log-add")


## ----sel_SNPs------------------------------------------------------------
sel <- labels(dd.s)[additive(ans)<0.1]
dd.sel <- dd.s[,sel]
head(dd.sel)


## ----create_df-----------------------------------------------------------
dd.sel <- data.frame(lapply(dd.sel, additive))
dd.end <- data.frame(casecontrol=dd.s$casecontrol, dd.sel)
head(dd.end)


## ----stepAIC-------------------------------------------------------------
dd.end.complete <- dd.end[complete.cases(dd.end),]
mod <- stepAIC(glm(casecontrol ~ ., dd.end.complete,
                   family="binomial"),
               method="forward", trace=0)
snps.score <- names(coef(mod))[-1]
snps.score


## ----print_gs------------------------------------------------------------
summary(mod)


## ----create_gs-----------------------------------------------------------
pos <- which(names(dd.end.complete)%in%snps.score)
names(dd.end.complete)
pos
score <- riskScore(mod, data=dd.end.complete, 
                      cGenPreds=pos,
                      Type="unweighted")
table(score)


## ----histGS, fig.cap='Distribution of the genetic score used to predict case/control status in the asthma example.'----
hist(score, col="gray90")


## ----effect_gs-----------------------------------------------------------
mod.lin <- glm(casecontrol~score, dd.end.complete,
               family="binomial")
exp(coef(mod.lin)[2])


## ----rocGS, fig.cap='ROC curve of the genetic score used to predict case/control status in the asthma example.'----
predrisk <- predRisk(mod.lin, dd.end.complete)
plotROC(data=dd.end.complete, cOutcome=1,
        predrisk = predrisk)


## ----load_ob-------------------------------------------------------------
library(snpStats)
path <- system.file("extdata", package="brgedata")
ob.plink <- read.plink(file.path(path, "obesity"))


## ----info_plink----------------------------------------------------------
names(ob.plink)


## ----get_geno_annot------------------------------------------------------
ob.geno <- ob.plink$genotypes
ob.geno

annotation <- ob.plink$map
head(annotation)

family <- ob.plink$fam
head(family)


## ----read_pheno----------------------------------------------------------
ob.pheno <- read.delim(file.path(path, "obesity.txt"))
head(ob.pheno)


## ----ids-----------------------------------------------------------------
rownames(ob.pheno) <- ob.pheno$id


## ----check_ids-----------------------------------------------------------
identical(rownames(ob.pheno), rownames(ob.geno))


## ----select_equals-------------------------------------------------------
ids <- intersect(rownames(ob.pheno), rownames(ob.geno))
geno <- ob.geno[ids, ]
ob <- ob.pheno[ids, ]
identical(rownames(ob), rownames(geno))

family <- family[ids, ]


## ----qc_snps-------------------------------------------------------------
info.snps <- col.summary(geno)
head(info.snps)


## ----qc_snp--------------------------------------------------------------
controls <- ob$obese ==0 & !is.na(ob$obese)
geno.controls <- geno[controls,]
info.controls <- col.summary(geno.controls)

use <- info.snps$Call.rate > 0.95 &
       info.snps$MAF > 0.05 &
       abs(info.controls$z.HWE < 3.3)    
mask.snps <- use & !is.na(use)

geno.qc.snps <- geno[ , mask.snps]
geno.qc.snps

annotation <- annotation[mask.snps, ]


## ----describe_removed_SNPs-----------------------------------------------
# number of SNPs removed for bad call rate 
sum(info.snps$Call.rate < 0.95)
# number of SNPs removed for low MAF
sum(info.snps$MAF < 0.05, na.rm=TRUE)
#number of SNPs that do not pass HWE test
sum(abs(info.controls$z.HWE > 3.3), na.rm=TRUE)    
# The total number of SNPs do not pass QC
sum(!mask.snps)


## ----qc_samples----------------------------------------------------------
info.indv <- row.summary(geno.qc.snps)
head(info.indv)


## ----HetX, fig.cap='Heterozygosity in chromosome X by gender provided in the phenotypic data.'----
geno.X <- geno.qc.snps[,annotation$chromosome=="23" & 
                         !is.na(annotation$chromosome)]
info.X <- row.summary(geno.X)
mycol <- ifelse(ob$gender=="Male", "gray40", "gray80")
plot(info.X$Heterozygosity, col=mycol, 
     pch=16, xlab="Individuals", 
     ylab="Heterozygosity in chromosome X")
legend("topright", c("Males", "Females"), col=mycol,
       pch=16)


## ----sex_disc------------------------------------------------------------
sex.discrep <- (ob$gender=="Male" & info.X$Heterozygosity > 0.2) |  
               (ob$gender=="Female" & info.X$Heterozygosity < 0.2)   


## ----compute_FHet--------------------------------------------------------
MAF <- col.summary(geno.qc.snps)$MAF
callmatrix <- !is.na(geno.qc.snps)
hetExp <- callmatrix %*% (2*MAF*(1-MAF))
hetObs <- with(info.indv, Heterozygosity*(ncol(geno.qc.snps))*Call.rate)
info.indv$hetF <- 1-(hetObs/hetExp)

head(info.indv)


## ----Het2, echo=FALSE, fig.cap='Heterozygosity computed using F statistic (left panel) and using row.summary function (right panel). The horizontal dashed line shows a suggestive value to detect individuals with outlier heterozygosity values.'----
par(mfrow=c(1,2))
plot(info.indv$hetF, ylim=c(-0.15,0.15), ylab="F Heterozygosity", 
     pch=16, cex=0.6, col="gray80")
abline(h=0.1, col="red", lty=2)
abline(h=-0.1, col="red", lty=2)
o <- info.indv$hetF>0.1
wordcloud::textplot(c(1:nrow(info.indv))[o], info.indv$hetF[o], rownames(info.indv)[o], new=FALSE, cex=0.8)
plot(info.indv$Heterozygosity, ylab="Heterozygosity", 
     pch=16, cex=0.6, col="gray80")
o <- info.indv$Heterozygosity<0.32
wordcloud::textplot(c(1:nrow(info.indv))[o], info.indv$Heterozygosity[o], rownames(info.indv)[o], new=FALSE, cex=0.8)


## ----ibd-----------------------------------------------------------------
library(SNPRelate)

# Transform PLINK data into GDS format
snpgdsBED2GDS("obesity.bed", 
              "obesity.fam",
              "obesity.bim", 
              out="obGDS")
genofile <- snpgdsOpen("obGDS")

#Prune SNPs for IBD analysis
set.seed(12345)
snps.qc <- colnames(geno.qc.snps)
snp.prune <- snpgdsLDpruning(genofile, ld.threshold = 0.2,
                          snp.id = snps.qc)
snps.ibd <- unlist(snp.prune, use.names=FALSE)


## ----ibd2----------------------------------------------------------------
ibd <- snpgdsIBDMoM(genofile, kinship=TRUE,
                    snp.id = snps.ibd,
                    num.thread = 1)
ibd.kin <- snpgdsIBDSelection(ibd) 
head(ibd.kin)


## ----related-------------------------------------------------------------
ibd.kin.thres <- subset(ibd.kin, kinship > 0.1)
head(ibd.kin.thres)


## ----qc_related----------------------------------------------------------
ids.rel <-  related(ibd.kin.thres) 
ids.rel


## ----qc_indiv------------------------------------------------------------
use <- info.indv$Call.rate > 0.95 &
       abs(info.indv$hetF) < 0.1 &
       !sex.discrep &
       !rownames(info.indv)%in%ids.rel
mask.indiv <- use & !is.na(use)
geno.qc <- geno.qc.snps[mask.indiv, ]

ob.qc <- ob.pheno[mask.indiv, ]
identical(rownames(ob.qc), rownames(geno.qc))


## ----summary_qc----------------------------------------------------------
# number of individuals removed to bad call rate
sum(info.indv$Call.rate < 0.95)
# number of individuals removed for heterozygosity problems
sum(abs(info.indv$hetF) > 0.1)
# number of individuals removed for sex discrepancies
sum(sex.discrep)
# number of individuals removed to be related with others
length(ids.rel)
# The total number of individuals that do not pass QC
sum(!mask.indiv)


## ----pca-----------------------------------------------------------------
pca <- snpgdsPCA(genofile, sample.id = rownames(geno.qc),  
                           snp.id = snps.ibd, 
                           num.thread=1)


## ----pca_plot, fig.cap="1st and 2nd principal components of obesity GWAS data example."----
with(pca, plot(eigenvect[,1], eigenvect[,2], 
               xlab="1st Principal Component", 
               ylab="2nd Principal Component", 
               main = "Ancestry Plot",
               pch=21, bg="gray90", cex=0.8))


## ----ad_pc---------------------------------------------------------------
ob.qc <- data.frame(ob.qc, pca$eigenvect[, 1:5])


## ----close_gds-----------------------------------------------------------
closefn.gds(genofile)


## ----gwas----------------------------------------------------------------
res <- single.snp.tests(obese, data=ob.qc, 
                               snp.data=geno.qc)
res[1:5,]


## ----gwas_quant----------------------------------------------------------
res.quant <- snp.rhs.tests(age ~ 1,  data=ob.qc, 
                           snp.data=geno.qc,
                           family="Gaussian")
head(res.quant)


## ----qqPlot, fig.cap="QQ-plot corresponding to obesity GWAS data example.", fig.height=3.5, fig.width=3.5----
chi2 <- chi.squared(res, df=1)
qq.chisq(chi2)


## ----gwas_adh------------------------------------------------------------
res.adj <- snp.rhs.tests(obese ~ X1 + X2 + X3 + X4 + X5, 
                         data=ob.qc, snp.data=geno.qc)
head(res.adj)


## ----p-value-------------------------------------------------------------
pval.log10 <- -log10(p.value(res.adj))


## ----manhattan, fig.cap="Manhattan plot of obesity GWAS data example.", fig.height=6.5, fig.width=9----
library(tidyverse)
library(ggplot2)
library(ggrepel)
# Create the required data frame
pvals <- data.frame(SNP=annotation$snp.name, 
                    CHR=annotation$chromosome,
                    BP=annotation$position,
                    P=p.value(res.adj))
# missing data is not allowed
pvals <- subset(pvals, !is.na(CHR) & !is.na(P)) 

manhattanPlot(pvals, color=c("gray90", "gray40"))


## ----subset_pvals--------------------------------------------------------
topPvals <- subset(pvals, P<10e-5)
topSNPs <- as.character(topPvals$SNP)


## ----export_SNPs---------------------------------------------------------
# subset top SNPs
geno.topSNPs <- geno.qc[, topSNPs]
geno.topSNPs
# export top SNPs
write.SnpMatrix(geno.topSNPs, file="topSNPs.txt")
# import top SNPs
ob.top <- read.delim("topSNPs.txt", sep="")
# add phenotypic information(ids are in the same order)
ob.top <- cbind(ob.top, ob.qc)
# prepare data for SNPassoc (SNPs are coded as 0,1,2)
ii <- grep("^rs", names(ob.top))
ob.top.s <- setupSNP(ob.top, colSNPs = ii, 
                     name.genotypes=c(0,1,2))
# run association (all)
WGassociation(obese, ob.top.s)
# run association (one)
association(obese ~ rs2769689, ob.top.s)


## ----postgwas, fig.cap="Genomic annotation of SNPs from asthma data example in chr7:34.5-35.0Mb region joint with SNPs described in the GWAS catalog, eQTLs annotated in Pickrell et al. (2010) and dsQTL from Degner et al. 2012."----
library(gwascat)
library(rtracklayer)
library(Gviz)

data(gwascatalog, package = "brgedata")

gff <- "chr7_33000000_35000000.gff3"

chr <- "chr7"
ir <- IRanges(start=33e6, end=35e6)
region <- GRanges(seqnames=chr, ir)

pos <- snpInfo$chrom_start[snpInfo$chr_name=="7"]
snps.loc <- GRanges(seqnames=chr,
                     IRanges(start=pos,
                     width=1))
 
c7 <- import(gff)

# annotate GWAS catalog
annot <- gwcex2gviz(basegr = gwascatalog, 
                    contextGR=region, plot.it=FALSE)

# annotate dsQTL
c7ds <- c7[ which(mcols(c7)$type == "Degner_dsQTL") ]
ds <- DataTrack( c7ds, chrom=chr, genome="hg19", 
                  data="score", name="dsQTL")

# annotate eQTL
c7e <- c7[ which(mcols(c7)$type == "Pickrell_eqtl") ]
eq <- AnnotationTrack(c7e,  chrom=chr, genome="hg19", 
                          name="eQTL")

# annotate SNPs from association study
snps2 <- AnnotationTrack(snps.loc,  chrom=chr, 
                         genome="hg19", name="asthma SNPs")

# Create plot
displayPars(eq)$col <- "black"
displayPars(ds)$col <- "red"
displayPars(snps2)$col <- "darkgreen"
integ <- list(annot[[1]], eq, ds, snps2, annot[[2]], annot[[3]])
plotTracks(integ)


## ----diseases------------------------------------------------------------
traits <- subsetByOverlaps(gwascatalog, 
                           region)$DISEASE.TRAIT
length(traits)
head(traits, n=12)


## ----load_CNVassoc, cache=FALSE------------------------------------------
library(CNVassoc)


## ----load_snptest--------------------------------------------------------
data(SNPTEST, package="CNVassocData")
dim(cases)
dim(controls)


## ----show_cases----------------------------------------------------------
cases[1:10,1:11]


## ----show_controls-------------------------------------------------------
controls[1:10,1:8]


## ----put_data_SNPTEST----------------------------------------------------
nSNP <- nrow(cases)
probs <- lapply(1:nSNP, function(i) {
  snpi.cases <- matrix(as.double(cases[i, 6:ncol(cases)]),
                       ncol = 3, byrow = TRUE)
snpi.controls <- matrix(as.double(controls[i, 6:ncol(controls)]),
                        
                        ncol = 3, byrow = TRUE)
return(rbind(snpi.cases, snpi.controls))
})


## ----set_casecon---------------------------------------------------------
casecon <- rep(1:0, c(500, 500))


## ----fit_models----------------------------------------------------------
pvals <- multiCNVassoc(probs, formula = "casecon ~ CNV", 
                       model = "add", num.copies = 0:2, 
                       cnv.tol = 0.001)


## ----correct_SNPTEST-----------------------------------------------------
pvalsBH <- getPvalBH(pvals)
head(pvalsBH)


## ----table_pvalues_SNPTEST-----------------------------------------------
table(cut(pvalsBH[, 2], c(-Inf, 1e-3, 1e-2, 0.05, 0.1, Inf)))


## ----remove_dll, echo=FALSE----------------------------------------------
detach("package:snpStats", unload=TRUE)
detach("package:gwascat", unload=TRUE)


## ----clean_snps, echo=FALSE, results='hide'------------------------------
#remove variables
rm(list=ls())

#garbage collection
gc(verbose = FALSE)

