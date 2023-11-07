#' get variance components
#'
#' @param model.object a model object from a call to lmer()
#'
#' @return a data.frame of percent variance attributable to each model term
#' @export
#'
#' @examples
#' fm1 <- lme4::lmer(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
#' get.varcomps(fm1)
get.varcomps <- function(model.object){
  stopifnot(inherits(model.object, "lmerMod"))
  m <- as.data.frame(lme4::VarCorr(model.object))
  t <- data.frame(m$grp, pct.var = sapply(m$vcov, function(x) x / sum(m$vcov)))
  return(t)
}

