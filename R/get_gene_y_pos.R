#' Title Generate equally spaced y values
#'
#' @param from start value
#' @param to end value
#' @param length.out number of things
#'
#' @return vector of values
#' @export
#'
#' @examples
#' # no example
get_gene_y_pos <- function(from, to, length.out) {
  length.out <- length.out + 2
  result <- seq(from, to, length.out = length.out)
  result <- result[-c(1,length(result))]
  return(result)
}
