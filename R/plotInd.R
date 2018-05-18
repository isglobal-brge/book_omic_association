
plotInd <- function(x, group, ax1=1, ax2=2, col.list, tab=1, ...){
  if (!inherits(x, "rgcca") & !inherits(x, "sgcca"))
    stop("x must be an object of class 'rgcca'")
  ll <- levels(group)
  levs <- length(ll)
  comp1 <- x$Y[[tab]][, ax1]
  comp2 <- x$Y[[tab]][, ax2]
  if (missing(col.list)){
   mycols <-  c("red", "blue", "green", "orange", "violet", sample(colors()))  
   col.list <- mycols[1:levs]
  } 
  if (length(col.list)!= levs)
    stop("'col.list' length should be equal to the levels of grouping variable")
  cols <- as.character(factor(group, labels=col.list))
  plot(comp1, comp2, type="n", ...)
  points(comp1, comp2, pch=16, col=cols)
  abline(h=0, lty=2)
  abline(v=0, lty=2)
  grid(lty=3, col="gray80")
  legend("bottomright", legend=ll, pch=16, col=col.list)
  
}
