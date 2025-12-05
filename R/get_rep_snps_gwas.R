#' Get snps that represent LD blocks using pvalues
#'
#' @param bed.file path to bedfile, with .bed extension, need .bim and .fam in directory too
#' @param gwas.res.table table that contains pvalues, will be subset
#' @param pval.column name of pvalue column
#' @param r2thresh LD threshold for clumping
#' @param window.size size of window for clumping
#'
#' @return gwas results table subsetted to include only representative snps.
#' @export
#'
#' @examples
#' # no example
get_rep_snps_gwas <- function(bed.file, gwas.res.table, pval.column, r2thresh = .2, window.size = 1000){
  bed.in <- bigsnpr::bed(bed.file)
  x <- bigsnpr::bed_clumping(obj.bed = bed.in,
                             S = -log10(gwas.res.table[,pval.column]),
                             thr.r2 = r2thresh,
                             size = window.size)
  out <- gwas.res.table[x,]
  return(out)
}
