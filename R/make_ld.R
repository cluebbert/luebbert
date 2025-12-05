#' Make a table of LD around a certain snp
#'
#' @param plink.path path to plink 1.9 executable
#' @param snp.name Name of snp around which to calculate LD
#' @param window Window around snp in KB
#' @param bedfile bedfile path, no .bed extension
#' @param in.dir directory of bedfile, include trailing "/"
#' @param out.filepath name and path of outfile
#'
#' @return outputs table of LD, for use within plotting function mostly
#' @export
#'
#' @examples
#' # No example
make_ld <- function(plink.path = "~/bin/plink",
                    snp.name,
                    window,
                    bedfile,
                    in.dir,
                    out.filepath = "ld_out_temp") {
  system(
    paste0(
      plink.path,
      " --silent --bfile ",
      in.dir,
      bedfile,
      " --r2 --ld-snp ",
      snp.name,
      " --ld-window-kb ",
      window,
      " --ld-window 99999 --ld-window-r2 0 --out ", out.filepath
    )
  )
}
