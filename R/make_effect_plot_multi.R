# make volcano plot

# ------------------------------------------------------------------------\
# unexported functions --------
# ------------------------------------------------------------------------\

get_bp_from_id <- function(marker.ID){
  as.numeric(stringr::str_extract(marker.ID, "-(.*)", group = 1))
}

# ------------------------------------------------------------------------\
# main function --------
# ------------------------------------------------------------------------\

# table should include CHR, POS, PVAL, EFF

#' Make volcano style plot of pvalue and effect size of gwas outputs. Takes gwas output from multiple gwas models.
#'
#' @param gwas.res table of all gwas results, should contain columns (CHR, POS, PVAL, EFF), corresponding to (chromosome, physical position, pvalue, beta/effect size). A column that indicates the unique gwas model a result was generated from indicated via `unique.gwas.model.variable` is required
#' @param qtl.df table with same columns that includes only significant hits in a qtl. QTL are typically defined as hits grouped by LD by something like `plink --clump`. A column that indicates the unique gwas model a result was generated from indicated via `unique.gwas.model.variable` is required
#' @param pvals.in.log boolean, are pvalues in input data.frames in -log10(p)?
#' @param plot.r2.thresh minimum LD with qtl snps to plot snps colored by LD. SNPs below are plotted in grey.
#' @param ld.list output of [luebbert::get_ld_in_window_multi]
#' @param window kilobases on either side of top QTL snp to plot
#' @param sig.line -log10(p) value to draw line on plot
#' @param unique.gwas.model.variable string indicating column name that corresponds to unique gwas models in the qtl.df
#' @param sig.hit.color.var Optional. String indicating a column that will color significant hits.
#' @param sig.hit.shape.var Optional. String indicating a column that will change shape of significant hits.
#' @param shape.scale Optional. A defined shape scale to keep plotting of shapes consistent when all factors are note present. See "another great function" I made to do this.
#' @param include.legend boolean, should legend be included?
#'
#' @returns
#' GGplot of volcano style plot. points colored by maximum R2 to markers in the qtl.df. Only results from gwas models present in the qtl are plotted.
#' @export
#'
#' @examples
#' # work in progress
make_effect_plot_multi <- function(gwas.res,
                                   qtl.df,
                                   pvals.in.log = T,
                                   plot.r2.thresh = .2,
                                   ld.list,
                                   window,
                                   sig.line,
                                   unique.gwas.model.variable = NULL,
                                   sig.hit.color.var = NULL,
                                   sig.hit.shape.var = NULL,
                                   shape.scale = NULL,
                                   include.legend = FALSE){
  # define plot colors
  my.colors <- viridis::viridis_pal(option = "turbo")(5)

  chrom <- unique(qtl.df$CHR)

  gwas.sub_alltraits <- gwas.res %>%
    filter(.data$CHR == chrom) %>%
    mutate(marker.ID = paste(.data$CHR, .data$POS, sep = "-")) %>%
    # filter(between(physical.pos, this.pos - window * 1000, this.pos + window * 1000)) %>%
    left_join(ld.list$table, by = "marker.ID") %>%
    filter(!is.na(.data$R2))

  # Get only hits that are from traits in the clump
  if(!is.null(unique.gwas.model.variable)){
    gwas.sub <- gwas.sub_alltraits %>%
      filter(.data[[unique.gwas.model.variable]] %in% unique(qtl.df[,unique.gwas.model.variable]))
  } else {
    gwas.sub <- gwas.sub_alltraits
  }

  # make manhattan
  plot.df <- gwas.sub %>%
    # alpha scale
    mutate(how.to.plot = case_when(.data$R2 > plot.r2.thresh ~ 1,
                                   TRUE ~ .4)) %>%
    # color scale
    mutate(plot.R2 = case_when(.data$R2 < plot.r2.thresh ~ NA,
                               TRUE ~ .data$R2))

  # change pvalue if needed
  if(!pvals.in.log){
    plot.df <- plot.df %>%
      mutate(PVAL = -log10(.data$PVAL))
  }

  # how far to spread labels past ends, in percentage
  y.spread.expansion <- .1
  y.spread.factors <- c(1 + y.spread.expansion, 1 - y.spread.expansion)
  y.spread.factor.window <- (window * y.spread.expansion) * 1000

  # plot limits
  this.pos <- get_bp_from_id(ld.list$key.snp)
  plot.limits <- c(this.pos + window * 1000, this.pos - window * 1000)
  plot.limits.ex <- c(plot.limits[1] + y.spread.factor.window, plot.limits[2] - y.spread.factor.window)

  # make the plot

  effect.plot <-
    ggplot(aes(x = .data$EFF, y = .data$PVAL), data = plot.df) +
    geom_point(aes(color = .data$plot.R2, alpha = .data$how.to.plot), shape = 16, show.legend = include.legend) +
    # geom_point(aes(x = estim, y = logpval, fill = trait, shape = experiment), data = this.clump.df, size = 4,
    #            show.legend = F) +
    # shapeScale +
    scale_color_stepsn(colors = my.colors, name = "R2") +
    geom_hline(yintercept = sig.line) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(y = bquote(-log[10](p-value)),
         x = "Effect Estimate") +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'gray')

  if(!is.null(sig.hit.color.var) & !is.null(sig.hit.shape.var)){
    effect.plot <- effect.plot +
      geom_point(aes(x = .data$EFF, y = .data$PVAL, fill = .data[[sig.hit.color.var]], shape = .data[[sig.hit.shape.var]]), data = qtl.df, size = 4) +
      guides(fill = guide_legend(override.aes = list(shape=21)))
    if(!is.null(shape.scale)){
      effect.plot <- effect.plot +
        shape.scale
    }
  } else if(!is.null(sig.hit.shape.var)){
    effect.plot <- effect.plot +
      geom_point(aes(x = .data$EFF, y = .data$PVAL, shape = .data[[sig.hit.shape.var]]), data = qtl.df, size = 4)
    if(!is.null(shape.scale)){
      effect.plot <- effect.plot +
        shape.scale
    }
  } else if(!is.null(sig.hit.color.var)){
    effect.plot <- effect.plot +
      geom_point(aes(x = .data$EFF, y = .data$PVAL, fill = .data[[sig.hit.color.var]]), data = qtl.df, size = 4, shape = 21)
  } else {
    effect.plot <- effect.plot +
      geom_point(aes(x = .data$EFF, y = .data$PVAL), fill = 'black', data = qtl.df, size = 4, show.legend = F)

  }

  if(!include.legend){
    effect.plot <- effect.plot +
      theme(legend.position="none")
  }

  return(effect.plot)

}

