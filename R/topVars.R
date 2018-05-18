
topVars <- function(x, axis = 1, end = "pos", topN = 5){
  
  if (!inherits(x, "rgcca"))
    stop("x must be an object of class 'rgcca'")
  mm <- match(end, c("pos", "neg"), nomatch = NA)
  if (mm==2)
    side <- FALSE
  else
    side <- TRUE
  if (is.na(mm))
    stop("'end' must be 'pos' or 'neg'")
  a <- x$a
  ntables <- length(a)
  tops <- NULL
  for (i in 1: ntables){
    a.i <- a[[i]]
    top.i <- head(rownames(a.i)[order(a.i[,axis], decreasing = side)], n=topN)
    tops <- cbind(tops, top.i)
  }
  colnames(tops) <- paste0("table_", 1:ntables)
  rownames(tops) <- paste0("top_", 1:topN)
  tops
}  
