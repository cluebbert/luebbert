#' Make a manhattan plot from table of gwas results
#'
#' @param df data.frame with chromosome, position and pvalue columns
#' @param chrom.col character - name of chromosome column
#' @param pos.col character - name of position column
#' @param pval.col character - name of pvalue column
#' @param pval.in.log boolean - is pvalue column log10 or not
#' @param bonf.cutoff significance cutoff
#' @param smoothing.cutoff minimum pvalue to plot to reduce size of plot
#' @param highlight.df optional - dataframe in same format as input that has snps to highlight
#'
#' @return manhattan plot
#' @export
#'
#' @examples
make_manhattan <- function(df, chrom.col, pos.col, pval.col, pval.in.log = F, bonf.cutoff, smoothing.cutoff, highlight.df = NULL){
  # format input
  gwas.dat <- df %>%
    select(all_of(c(chrom.col, pos.col, pval.col))) %>%
    setNames(c("chromosome", "physical.pos", "pval"))

  # gwas.dat <- df %>%
  #   rename(chromosome = chrom.col,
  #          physical.pos = pos.col,
  #          pval = pval.col)
  if(!pval.in.log){
    gwas.dat <- gwas.dat %>%
      mutate(logpval = -log10(.data$pval))
  } else {
    gwas.dat <- gwas.dat %>%
      mutate(logpval = .data$pval)
  }

  offset.table <- make_manhattan_offsets(gwas.dat, "chromosome", "physical.pos")

  # cutoff yaxis here
  # smoothing_cutoff <- 4

  # offset xaxis values and add alternating color indicator
  plot.df <- left_join(gwas.dat, offset.table,
                       by = "chromosome") %>%
    filter(.data$logpval > smoothing.cutoff) %>%
    # mutate(x.value = (physical.pos - 1e4) + (offset - 1e4)) %>%
    mutate(x.value = .data$physical.pos + .data$offset) %>%
    arrange(.data$x.value) %>%
    mutate(pt.color = case_when(.data$chromosome %% 2 == 1 ~ "odd",
                                TRUE ~ "even")) %>%
    tibble::rownames_to_column(var = "num") %>%
    mutate(num = as.numeric(.data$num))

  # make chromosome names on axis
  xbreaks.df <- plot.df %>%
    group_by(.data$chromosome) %>%
    mutate(x.break = mean(.data$x.value, na.rm = T))
  xbreaks <- sort(unique(xbreaks.df$x.break))
  xlabs <- sort(unique(plot.df$chromosome))

  # make nice yaxis label
  logobj <- expression(paste("-log"[10], plain("(p-value)")))

  p <-
    ggplot(plot.df, aes(x = .data$x.value, y = .data$logpval)) +
    geom_point(aes(color = .data$pt.color)) +
    geom_hline(yintercept = bonf.cutoff,
               linetype = "dashed",
               color = "red",
               size = 1.5) +
    scale_color_manual(values = c("gray38", "gray61")) +
    theme_bw() +
    scale_x_continuous(name = "Chromosome",
                       breaks = xbreaks,
                       labels = xlabs) +
    theme(legend.position = "none",
          axis.text = element_text(size = 30),
          axis.title.y = element_text(size = 32),
          axis.title.x = element_text(size = 32),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank()) +
    labs(y = logobj)


  return(list(plot = p, table = plot.df))
}
