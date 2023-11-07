get.varcomps <- function(model.object){
  stopifnot(inherits(model.object, "lmerMod"))
  m <- as.data.frame(lem4::VarCorr(model.object))
  t <- data.frame(m$grp, pct.var = sapply(m$vcov, function(x) x / sum(m$vcov)))
  return(t)
}
