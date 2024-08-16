#' Title get a gene's distance from a given snp
#'
#' @param s Snp physical position
#' @param start gene start
#' @param end gene end
#'
#' @return minimum distance the start or end of a gene is to the snp
#' @export
#'
#' @examples
#' # no example
get_gene_dist_from_snp <- function(s, start, end){
  a <- abs(s - start)
  b <- abs(s - end)
  return(min(a, b))
}
