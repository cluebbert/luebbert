#' Make a table of LD around a certain snp
#'
#' @param snp.name Name of snp around which to calculate LD
#' @param window Window around snp in KB
#' @param bedfile bedfile path, no .bed extension
#' @param out.filepath path of outfile
#'
#' @return outputs table of LD, for use within plotting function mostly
#' @export
#'
#' @examples
#' # No example
make_ld <- function(snp.name,
                    window = 1000,
                    bedfile,
                    out.filepath = "./ld_out_temp"){
  system(paste0("~/bin/plink --bfile ",
                bedfile,
                " --r2 --ld-snp ", snp.name,
                " --ld-window-kb ", window,
                " --ld-window 99999 --ld-window-r2 0 --out ", out.filepath))
}
