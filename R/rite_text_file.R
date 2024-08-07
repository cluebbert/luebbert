
#' make text file with correct encoding for use on linux system
#'
#' @param filename name and path of output file (without trailing .txt)
#' @param data vector to be output to text file
#'
#' @return vector as text file, each element on separate line
#' @export
#'
#' @examples # no example
#'
rite_text_file <- function(filename, data){
  output.file <- file(paste0(filename,".txt"), "wb")
  f <- data
  cat(f, file = output.file, sep = "\n")
  close(output.file)
}
