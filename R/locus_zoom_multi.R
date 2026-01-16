
# ------------------------------------------------------------------------\
# main function --------
# ------------------------------------------------------------------------\

#' Make a locus zoom style plot with points colored by R2 to a qtl and nearby genes annotated. Takes gwas output from multiple gwas models.
#'
#' @param gwas.res table of all gwas results, should contain columns (`CHR`, `POS`, `PVAL`), corresponding to (chromosome, physical position, and pvalue).
#' @param qtl.df table with same columns that includes only significant hits in a qtl. QTL are typically defined as hits grouped by LD by something like `plink --clump`.
#' @param annotation.table table with annotations with columns (geneID, CHR, start, end, annotation). start and end correspond to base-pair coordinates of start and end of gene. CHR is chromosome of gene.
#' @param pvals.in.log boolean, are pvalues in input data.frames in -log10(p) format?
#' @param plink.path path to plink 1.9 executable.
#' @param geno.bed prefix of genotype files in plink (bed/bim/fam) format.
#' @param geno.bed.filename prefix of genotype files in plink (bed/bim/fam) format.
#' @param geno.bed.directory directory where genotype files are located.
#' @param plot.r2.thresh minimum LD with qtl snps to plot snps colored by LD. SNPs below are plotted in grey.
#' @param window kilobases on either side of top QTL snp to plot.
#' @param sig.line -log10(p) value to draw line on plot.
#' @param plot.title Optional. Title of plot.
#' @param unique.gwas.model.variable string indicating column name in gwas.res and qtl.df that corresponds to unique gwas models (e.g. phenotype within a location). See [luebbert::get_ld_in_window_multi] for usage in LD calculation.
#' @param sig.hit.color.var Optional. String indicating a column in gwas.res and qtl.df that will color significant hits. (e.g. a phenotype or location)
#' @param sig.hit.shape.var Optional. String indicating a column in gwas.res and qtl.df that will change shape of significant hits. (e.g. a phenotype or location)
#' @param shape.scale Optional. A defined shape scale to keep plotting of shapes consistent when all factors are note present. See another great function I made to do this.
#' @param single.variable.color Optional. a single color to plot qtl.df points as if there is no variable of interest
#' @param include.gene.id boolean, include geneID in gene annotations or not
#' @param plot.effect boolean, include plot of pvalues vs effect size? If TRUE, `gwas.res` and `qtl.df` must have column `EFF` that conatins gwas effect sizes.
#'
#' @returns
#' GGplot of locus zoom style plot with points colored by maximum R2 to markers in the qtl.df. Only results from `unique.gwas.model.variable` present in the qtl are plotted. Alongside are annotations of nearby genes.
#' @export
#'
#' @examples
#' # work in progress
locus_zoom_multi <- function(gwas.res,
                             qtl.df,
                             annotation.table,
                             pvals.in.log = T,
                             plink.path,
                             geno.bed.filename,
                             geno.bed.directory = "/.",
                             plot.r2.thresh = .2,
                             window,
                             sig.line,
                             plot.title = "",
                             unique.gwas.model.variable,
                             sig.hit.color.var = NULL,
                             sig.hit.shape.var = NULL,
                             shape.scale = NULL,
                             single.variable.color = NULL,
                             include.gene.id = FALSE,
                             plot.effect = TRUE){

  # make ld table
  ld.list <- get_ld_in_window_multi(qtl.df,
                                    window,
                                    plink.path,
                                    geno.bed.filename,
                                    geno.bed.directory,
                                    pvals.in.log,
                                    grouping.variable = unique.gwas.model.variable)

  # make manhattan
  man <- make_sideways_manhattan_multi(gwas.res,
                                       qtl.df,
                                       pvals.in.log,
                                       plot.r2.thresh,
                                       ld.list,
                                       window,
                                       sig.line,
                                       unique.gwas.model.variable,
                                       sig.hit.color.var,
                                       sig.hit.shape.var,
                                       single.variable.color,
                                       shape.scale)

  # make annotation
  anno <- make_gene_annotation_plot(annotation.table,
                                    middle.snp = ld.list$key.snp,
                                    window,
                                    include.id = include.gene.id)

  # make effect if you want
  if(plot.effect){
    effect.plot <- make_effect_plot_multi(gwas.res,
                                          qtl.df,
                                          pvals.in.log,
                                          plot.r2.thresh,
                                          ld.list,
                                          window,
                                          sig.line,
                                          unique.gwas.model.variable,
                                          sig.hit.color.var,
                                          sig.hit.shape.var,
                                          shape.scale)

  }

  chrom <- unique(qtl.df$CHR)

  # make final plot
  if(plot.effect){
    out <-
      man + theme(plot.margin = ggplot2::margin(c(6, 6, 6, 120), unit= "points")) +
      patchwork::inset_element(effect.plot, left = .0, bottom = 0, top = .3, right = .3, align_to = "full") +
      anno +
      patchwork::plot_layout(nrow = 1, widths = c(4,6)) +
      patchwork::plot_annotation(title = plot.title,
                      subtitle = paste("Chromosome", chrom))

  } else {
    out <-
      man + theme(plot.margin = ggplot2::margin(6, 6, 6, 6, unit= "points")) +
      anno +
      patchwork::plot_layout(nrow = 1, widths = c(4,6)) +
      patchwork::plot_annotation(title = plot.title,
                      subtitle = paste("Chromosome", chrom))
  }

  return(out)
}
