# make sideways manhattan plot

# ------------------------------------------------------------------------\
# unexported functions --------
# ------------------------------------------------------------------------\

# extract position from marker.ID in the form "CHR-POS"
get_bp_from_id <- function(marker.ID){
  as.numeric(stringr::str_extract(marker.ID, "-(.*)", group = 1))
}

# ------------------------------------------------------------------------\
# main function --------
# ------------------------------------------------------------------------\

#' Make sideways manhattan plot for building locus zoom. Receives output from a single gwas model.
#'
#' @param gwas.res table of all gwas results, should contain columns (CHR, POS, PVAL), corresponding to (chromosome, physical position, and pvalue).
#' @param qtl.df table with same columns that includes only significant hits in a qtl. QTL are typically defined as hits grouped by LD by something like `plink --clump`
#' @param pvals.in.log boolean, are pvalues in input data.frames in -log10(p)?
#' @param plot.r2.thresh minimum LD with qtl snps to plot snps colored by LD
#' @param ld.list output of [luebbert::get_ld_in_window]
#' @param window kilobases on either side of top QTL snp to plot
#' @param sig.line -log10(p) value to draw line on plot
#'
#' @returns
#' GGplot of manhattan plot with points colored by maximum R2 to markers in the qtl.df.
#' @export
#'
#' @examples
#' # Work in progress
make_sideways_manhattan <- function(gwas.res,
                                    qtl.df,
                                    pvals.in.log = TRUE,
                                    plot.r2.thresh = .2,
                                    ld.list,
                                    window,
                                    sig.line){

  # define plot colors
  my.colors <- viridis::viridis_pal(option = "turbo")(5)

  chrom <- unique(qtl.df$CHR)

  gwas.sub <- gwas.res %>%
    filter(.data$CHR == chrom) %>%
    mutate(marker.ID = paste(.data$CHR, .data$POS, sep = "-")) %>%
    # filter(between(physical.pos, this.pos - window * 1000, this.pos + window * 1000)) %>%
    left_join(ld.list$table, by = "marker.ID") %>%
    filter(!is.na(.data$R2))
  # for rug plot
  marker.list.in.window <- ld.list$table %>%
    mutate(POS = get_bp_from_id(.data$marker.ID)) %>%
    select("marker.ID", "POS") %>%
    distinct()


  # make manhattan
  plot.df <- gwas.sub %>%
    # alpha scale
    mutate(how.to.plot = case_when(.data$R2 > plot.r2.thresh ~ 1,
                                   TRUE ~ .4)) %>%
    # color scale
    mutate(plot.R2 = case_when(.data$R2 < plot.r2.thresh ~ NA,
                               TRUE ~ R2))

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

  # Base plot
  man <-
    ggplot(aes(x = .data$POS, y = .data$PVAL), data = plot.df) +
    geom_point(aes(color = .data$plot.R2, alpha = .data$how.to.plot), shape = 16) +
    scale_alpha(guide = "none") +
    scale_color_stepsn(colors = my.colors, name = "R2") +
    # shape.scale +
    scale_x_reverse(
      limits = plot.limits.ex,
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      name = "-log10pval"
    ) +
    coord_flip() +
    scale_y_reverse() +
    theme(
      panel.background = element_rect(fill = "grey95"),
      axis.title.y = element_blank(),
      legend.position = "left",
      legend.justification = "top",
      panel.grid = element_blank()
    ) +
    labs(y = bquote(-log[10](p - value))) +
    geom_rug(aes(x = .data$POS),
             sides = 't',
             data = marker.list.in.window,
             inherit.aes = F) +
    geom_hline(yintercept = sig.line, linetype = 'dashed')

  return(man)
}
