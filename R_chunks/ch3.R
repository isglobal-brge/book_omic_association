## ----load_plink, eval=FALSE----------------------------------------------
## library(snpStats)
## ob.plink <- read.plink(bed = "obesity.bed",
##                        bim = "obesity.bim",
##                        fam = "obesity.fam")


## ----read_plink2, eval=FALSE---------------------------------------------
## snps <- read.plink("obesity")


## ----read_plink3---------------------------------------------------------
path <- system.file("extdata", package="brgedata")
snps <- read.plink(file.path(path, "obesity"))
names(snps)


## ----genotypes-----------------------------------------------------------
geno <- snps$genotypes
geno


## ----individuals---------------------------------------------------------
individuals <- snps$fam
head(individuals)


## ----annotation----------------------------------------------------------
annotation <- snps$map
head(annotation)


## ----subset_snps---------------------------------------------------------
annotationChr1 <- annotation[annotation$chromosome == "1" & 
                               !is.na(annotation$chromosome), ]
genoChr1 <- geno[, rownames(annotationChr1)]
genoChr1


## ----select_controls-----------------------------------------------------
individualsCtrl <- individuals[individuals$affected == 1, ]
genoCtrl <- geno[rownames(individualsCtrl), ]
genoCtrl


## ----genoquerych3--------------------------------------------------------
library(GEOquery)
gsm.expr <- getGEO("GSE63061", destdir = ".")[[1]]
gsm.expr <- gsm.expr


## ----exprgetdatach3------------------------------------------------------
expr <- exprs(gsm.expr)
dim(expr)
expr[1:5,1:5]


## ----phenogetdatach3-----------------------------------------------------
#get phenotype data
pheno <- phenoData(gsm.expr)
pheno
phenoDataFrame <- pData(gsm.expr)
phenoDataFrame[1:5,1:4]

#Alzheimer's case control variable
summary(phenoDataFrame$characteristics_ch1)


## ----fDatagetdatach3-----------------------------------------------------
probes <- fData(gsm.expr)
probes[1:5, 1:5]


## ----summexp_example-----------------------------------------------------
library(brgedata)
brge_methy
extends("GenomicRatioSet")


## ----get_assays----------------------------------------------------------
betas <- assay(brge_methy)
betas[1:5, 1:4]


## ----get_features--------------------------------------------------------
rowData(brge_methy)[,2:5]


## ----createGR------------------------------------------------------------
library(GenomicRanges)
gr <- GRanges(seqnames=c(rep("chr1", 4), rep("chr2", 4)),
              ranges = IRanges(start = c(1000, 1800, 5300, 7900,
                                         1300, 2100, 3400, 6700),
                               end =c(2200, 3900, 5400, 8100,
                                      2600, 3300, 4460, 6850)),
              strand = rep(c("+", "-"), 4),
              disease = c(rep("Asthma",4), rep("Obesity",4)))
gr


## ----gr1-----------------------------------------------------------------
gr[1]


## ----gr2-----------------------------------------------------------------
seqnames(gr)
seqnames(gr)[1] <- "chr2"
gr


## ----gr3-----------------------------------------------------------------
gr$gene_id <- paste0("Gene", 1:8)
gr


## ----gr4-----------------------------------------------------------------
#shift: move all intervals 10 base pair towards the end
shift(gr, 10)

#shift: move each intervals individually
shift(gr, seq(10,100, length=8))

#flank:  recover regions next to the input set. 
#        For a 50 base stretch upstream (negative value for
#        downstream)
flank(gr, 50)


## ----gr5-----------------------------------------------------------------
disjoin(gr, ignore.strand=TRUE)


## ----gr6-----------------------------------------------------------------
reduce(gr, ignore.strand=TRUE)


## ----gr7-----------------------------------------------------------------
coverage(gr)


## ----gr8-----------------------------------------------------------------
target <- GRanges(seqnames="chr1", 
                  range=IRanges(start=1200, 4000))
target
gr.ov <- findOverlaps(target, gr)
gr.ov


## ----gr9-----------------------------------------------------------------
gr[subjectHits(gr.ov)]


## ----names_rowranges-----------------------------------------------------
names(rowData(brge_methy))


## ----rowranges-----------------------------------------------------------
rowRanges(brge_methy)[, "genes"]


## ----sample_metadata-----------------------------------------------------
colData(brge_methy)


## ----sample_metadata_variable--------------------------------------------
brge_methy$sex


## ----get_males-----------------------------------------------------------
brge_methy[, brge_methy$sex == "male"]


## ----get_metadata--------------------------------------------------------
metadata(brge_methy)


## ----data_files----------------------------------------------------------
library(rexposome)
path <- system.file("extdata", package="rexposome")
description <- file.path(path, "description.csv")
exposures <- file.path(path, "exposures.csv")
phenotype <- file.path(path, "phenotypes.csv")


## ----read_exposome-------------------------------------------------------
expo <- readExposome(exposures = exposures,
                      description = description,
                      phenotype = phenotype,
                      exposures.samCol = "idnum",
                      description.expCol = "Exposure",
                      description.famCol = "Family",
                      phenotype.samCol = "idnum"
)


## ----show_expos----------------------------------------------------------
expo


## ----show_names----------------------------------------------------------
head(sampleNames(expo))


## ----show_expos_names----------------------------------------------------
head(exposureNames(expo))


## ----show_family---------------------------------------------------------
familyNames(expo)


## ----show_pheno----------------------------------------------------------
phenotypeNames(expo)


## ----get_fdata-----------------------------------------------------------
head(fData(expo))


## ----get_pheno_data------------------------------------------------------
head(pData(expo))


## ----get_expos-----------------------------------------------------------
expos(expo)[1:10, c("Cotinine", "PM10V", "PM25", "X5cxMEPP")]


## ----mae-----------------------------------------------------------------
library(MultiAssayExperiment)
objlist <- list(expression = brge_gexp, methylation = brge_methy)
mae <- MultiAssayExperiment(objlist)
mae


## ----show_mae------------------------------------------------------------
experiments(mae)
colData(mae)
sampleMap(mae)
metadata(mae)


## ----get_first_list_mae--------------------------------------------------
matrices <- assays(mae)
names(matrices)
matrices[[1]][1:4, 1:5]
matrices[[2]][1:4, 1:5]


## ----mae_samples---------------------------------------------------------
mae[ , c("x0001", "x0002")]


## ----mds-----------------------------------------------------------------
data(brge_methy, package="brgedata")
data(brge_gexp, pckage="brgedata")

library(MultiDataSet)
mds <- createMultiDataSet()
mds <- add_methy(mds, brge_methy)
mds <- add_genexp(mds, brge_gexp)
mds


## ----subset_gr-----------------------------------------------------------
range <- GRanges("chr1:1-999999999")
mds[, , range]


## ----subset_sample-------------------------------------------------------
mds[c("x0001", "x0002"), ]


## ----commonsamples-------------------------------------------------------
commonSamples(mds)


## ----load_breast_list----------------------------------------------------
data(breastMulti_list, package="brgedata")


## ----content_list--------------------------------------------------------
lapply(breastMulti_list, dim)


## ----create_methy--------------------------------------------------------
library(minfi)
methyBreast <- breastMulti_list$Methy
methyBreast <- makeGenomicRatioSetFromMatrix(data.matrix(methyBreast)) 
class(methyBreast)


## ----get_annotation_genes------------------------------------------------
library(biomaRt)
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
symbol.genes <- rownames(breastMulti_list$RNAseq)
annot <- getBM(attributes=c('hgnc_symbol',
                          'chromosome_name', 
                          'start_position', 
                          'end_position'),
             filters=c('hgnc_symbol'),
             values=symbol.genes,
             mart=ensembl)


## ----annot_format--------------------------------------------------------
names(annot) <- c("symbol", "chromosome", "start", "end")
annot$chromosome <- paste0("chr", annot$chromosome)
annot <- annot[!duplicated(annot[,1]),]
rownames(annot) <- annot$symbol
annot <- GenomicRanges::makeGRangesFromDataFrame(annot)


## ----summExp_breast------------------------------------------------------
genes <- intersect(names(annot), symbol.genes)
exprBreast <- breastMulti_list$RNAseq[genes, ]
exprBreast <- SummarizedExperiment(assays=list(genexpr=exprBreast),
                                   rowRanges=annot[genes,])
exprBreast


## ----omics_eset----------------------------------------------------------
miRNABreast <- ExpressionSet(data.matrix(breastMulti_list$miRNA))
miRNApreBreast <- ExpressionSet(data.matrix(breastMulti_list$miRNAprecursor))
proteBreast <- ExpressionSet(data.matrix(breastMulti_list$RPPA))
lohBreast <- ExpressionSet(breastMulti_list$LOH)
cnaBreast <- ExpressionSet(breastMulti_list$CNA)


## ----add_clin------------------------------------------------------------
breastMulti_list$clin <- droplevels(breastMulti_list$clin)

colData(methyBreast) <- colData(exprBreast) <-
  DataFrame(breastMulti_list$clin)
pData(miRNABreast) <- pData(miRNApreBreast) <-  pData(proteBreast) <-
  pData(lohBreast) <- pData(cnaBreast) <- breastMulti_list$clin


## ----MDSbreast-----------------------------------------------------------
library(MultiDataSet)
breastMulti <- createMultiDataSet()
breastMulti <- add_methy(breastMulti, methyBreast)
breastMulti <- add_rse(breastMulti, exprBreast, "expression")
breastMulti <- add_eset(breastMulti, miRNABreast, "miRNA", 
                        GRanges = NA)
breastMulti <- add_eset(breastMulti, miRNApreBreast, "miRNAprecursor", 
                        GRanges = NA)
breastMulti <- add_eset(breastMulti, proteBreast, "RPPA", 
                        GRanges = NA)
breastMulti <- add_eset(breastMulti, lohBreast, "LOH", 
                        GRanges = NA)
breastMulti <- add_eset(breastMulti, cnaBreast, "CNA", 
                        GRanges = NA)
breastMulti


## ----clean_ch3, echo=FALSE, results='hide'-------------------------------
rm(list=ls())
gc(verbose = FALSE)

