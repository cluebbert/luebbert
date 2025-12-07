# ld one trait gwas result

# ------------------------------------------------------------------------\
# unexported --------
# ------------------------------------------------------------------------\

# make_ld <- function(plink.path, snp.name, window, bedfile, in.dir){
#   system(paste0(plink.path, " --silent --bfile ", in.dir, bedfile, " --r2 --ld-snp ", snp.name, " --ld-window-kb ", window, " --ld-window 99999 --ld-window-r2 0 --out ld_out_temp"))
# }


# ------------------------------------------------------------------------\
# ld function --------
# ------------------------------------------------------------------------\


#' Make a table of maximum R2 of a list of snps to all snps in a window. Takes gwas output from a single gwas model.
#'
#' @param qtl.df table that includes only significant hits in a qtl with columns (CHR, POS, PVAL), corresponding to (chromosome, physical position, and pvalue). QTL are typically defined as hits grouped by LD by something like `plink --clump`
#' @param window kilobases on either side of top QTL snp to plot
#' @param plink.path path to plink 1.9 executable
#' @param geno.bed prefix of genotype files in plink (bed/bim/fam) format
#' @param in.dir directory where genotype files are located
#' @param pvals.in.log boolean, are pvalues in input data.frames in -log10(p)?
#'
#' @returns
#' Named list with 2 items
#'  - table: table with marker.IDs (CHR-POS) and maximum LD in R2 for each snp to the snps in the qtl.df
#'  - key.snp: marker.ID corresponding to the middle of the window
#' @export
#'
#' @examples
#' # Work in progress

get_ld_in_window <- function(qtl.df,
                             window,
                             plink.path,
                             geno.bed,
                             in.dir,
                             pvals.in.log,
                             verbose = T){

  # clean up top hits table
  this.clump.df <- qtl.df %>%
    mutate(marker.ID = paste(.data$CHR, .data$POS, sep = "-"))

  # ensure pval in correct format
  if(!pvals.in.log){
    this.clump.df$PVAL <- -log10(this.clump.df$PVAL)
  }

  # Might be a tie here sometimes if top snp has two traits with same pvalue, I think the [1] solves it
  key.snp <- this.clump.df$marker.ID[which.max(this.clump.df$PVAL)[1]]

  this.chrom <- unique(this.clump.df$CHR)

  # If multiple traits hit same key snp, need to use unique to get just 1 value
  this.pos <- unique(this.clump.df$POS[which(this.clump.df$marker.ID == key.snp)])
  this.snp.name <- key.snp

  # # get only the top snp for each grouping variable
  # this.top.df <- this.clump.df %>%
  #   group_by(group) %>%
  #   mutate(is.top.pval = PVAL == max(PVAL)) %>%
  #   filter(is.top.pval) %>%
  #   # there are some snps that have exactly the same pvalue
  #   mutate(snp.number = 1:n()) %>%
  #   filter(snp.number == 1) %>%
  #   select(-snp.number)

  ld.table_all <- data.frame()

  do.progress <- nrow(this.clump.df) > 1 & verbose

  if(verbose){
    message(paste("Calculating LD to", nrow(this.clump.df), "snps."))
  }
  if(do.progress){
    pb <- txtProgressBar(min = 1, max = nrow(this.clump.df), style = 3)
  }
  for(i in 1:nrow(this.clump.df)){

    this.snp.name <- this.clump.df$marker.ID[i]

    make_ld(plink.path, this.snp.name, window, geno.bed, in.dir)


    ld.table <- read.table("ld_out_temp.ld", header = T)
    ld.table_sub <- ld.table %>%
      select("marker.ID" = "SNP_B", "R2")
    ld.table_all <- bind_rows(ld.table_all, ld.table_sub)
    if(do.progress){
      setTxtProgressBar(pb, i)
    }
  }
  close(pb)

  ld.table_all_topR <- ld.table_all %>%
    group_by(.data$marker.ID) %>%
    arrange(desc(.data$R2)) %>%
    filter(1:n() == 1)

  return(list(table = ld.table_all_topR,
              key.snp = key.snp))

}
