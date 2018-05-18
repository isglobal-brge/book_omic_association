selectVars <- function(x, table=1, axis = 1, end = "pos"){
  if (!inherits(x, "sgcca"))
    stop("x must be an object of class 'sgcca'")
  mm <- match(end, c("pos", "neg"), nomatch = NA)
  if (mm==2)
    side <- FALSE
  else
    side <- TRUE
  if (is.na(mm))
    stop("'end' must be 'pos' or 'neg'")
  a.i <- x$a[[table]]
  a <- a.i[,axis]
  a <- a[(order(a, decreasing=side))]
  if (side)
   tops <- names(a)[a>0]
  else
    tops <- names(a)[a<0]
  tops
}  
  
  
  