#' write.csv with row.names = F as default
#'
#' @param ... extra options to base r write.csv() function
#' @param row.names save rownames or not, F is default
#'
#' @return saves csv with no rownames by default
#' @export
#'
#' @examples
ritecsv <- function(...,row.names=FALSE){
  utils::write.csv(..., row.names = row.names)
}
