#' Assign snps to clumps based on LD after gwas
#'
#' @param geno.bed character, bedfile path, no .bed extension
#' @param gwas.res data.frame, table of gwas results with columns (marker.ID, CHR, POS, PVAL)
#' @param pvals.in.log boolean, are pvalues in `gwas.res` in -log10(p) or not
#' @param window integer, window in kilobases either side of snp to look for snps in LD,
#' @param ld.thresh numeric, R2 threshold above which snps will be grouped
#'
#' @returns
#' table with columns marker.ID and clump_num. Clump_num indicates groupings, numberings start from the larges pvalue to smallest. May want to reassign afterwards to be along the genome.
#' @export
#'
#' @examples
#' # work in progress
make_clumps <- function(geno.bed,
                        gwas.res,
                        pvals.in.log = T,
                        window = 500,
                        ld.thresh = .5){

  # convert pvals to log if we need
  if(!pvals.in.log){
    gwas.res <- gwas.res %>%
      mutate(PVAL = -log10(PVAL))
  }

  # get window in basepairs from kilobases
  window.bp <- window * 1000

  # get snps sorted by logpvalue, highest to lowest
  snps.to.test <- gwas.res %>%
    arrange(desc(PVAL)) %>%
    pull(marker.ID) %>%
    unique()

  # make an output df and start a counter
  out <- data.frame()
  i <- 1
  # start a progress bar
  message("Creating clumps...")
  pb.end.value <- length(snps.to.test)
  pb <- txtProgressBar(min = 1, max = pb.end.value, style = 3)
  while(length(snps.to.test) > 0){

    this.snp.name <- snps.to.test[1]
    # check if any snps are in window
    this.snp.info <- filter(gwas.res, marker.ID == this.snp.name) %>%
      select(SNP = marker.ID, POS, CHR) %>%
      distinct()

    pos.range <- c(this.snp.info$POS - window.bp, this.snp.info$POS + window.bp)

    snps.in.range <- gwas.res %>%
      filter(CHR == this.snp.info$CHR) %>%
      filter(between(POS, pos.range[1], pos.range[2]))

    if(nrow(snps.in.range) < 2){
      this.out <- data.frame(snp = this.snp.info$SNP,
                             clump_num = i)
      out <- bind_rows(out, this.out)
      snps.to.test <- snps.to.test[-which(snps.to.test %in% this.out$snp)]
      i <- i+1
      print(paste(length(snps.to.test), "snps remaining."))

    } else {

      make.ld.table(this.snp.name, window, geno.bed)

      ld.table <- read.table("./ld_out_temp.ld", header = T)
      ld.table_sub <- ld.table %>%
        select("SNP" = "SNP_B", "R2") %>%
        filter(R2 >= ld.thresh) %>%
        filter(SNP %in% snps.to.test)

      this.out <- data.frame(marker.ID = ld.table_sub$SNP,
                             clump_num = i)
      out <- bind_rows(out, this.out)
      snps.to.test <- snps.to.test[-which(snps.to.test %in% this.out$marker.ID)]
      i <- i+1
      # make this into a progress bar
      # message(paste(length(snps.to.test), "snps remaining."))
      setTxtProgressBar(pb, pb.end.value - length(snps.to.test))
    }
  }

  return(out)
}

