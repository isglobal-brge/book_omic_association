## ----CRAN, eval=FALSE----------------------------------------------------
## install.package("nameOfPackage")


## ----CRANlib, eval=FALSE-------------------------------------------------
## library("nameOfPackage")


## ----geoinstall, eval=FALSE----------------------------------------------
## source("https://bioconductor.org/biocLite.R")
## biocLite("nameOfPackage")


## ----geoinstall2, eval=FALSE---------------------------------------------
## library(BiocManager)
## install("GEOquery")


## ----geoload, eval=FALSE-------------------------------------------------
## library("GEOquery")


## ----GitHub, eval=FALSE--------------------------------------------------
## library(devtools)
## install_github("nameOfRepository/nameOfPackage")


## ----GitHub2, eval=FALSE-------------------------------------------------
## library(devtools)
## install_github("isglobal-brge/nameOfPackage")


## ----install_brgech1, cache=FALSE, eval=FALSE----------------------------
## library(devtools)
## install_github("isglobal-brge/brgedata")


## ----load_brge, cache=FALSE, eval=FALSE----------------------------------
## library(brgedata)


## ----path_brge, cache=FALSE, eval=FALSE----------------------------------
## path <- system.file("extdata", package="brgedata")


## ----get_brgedata--------------------------------------------------------
data(asthma, package = "brgedata")
asthma[1:5, 1:10]


## ----ls_brgedata---------------------------------------------------------
library(brgedata)
ls("package:brgedata")

