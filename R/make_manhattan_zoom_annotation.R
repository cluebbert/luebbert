#' Title Make a manhattan in a region and label genes and LD for a table of top gwas hits
#'
#' @param out.dir where to output plots
#' @param plot.type png or pdf?
#' @param trait.name character string of name of trait
#' @param top.hits table of top hits with columns CHROM, POS, and PVAL column, p-value is log transformed in function so should not be pre-transformed
#' @param gwas.res table of all gwas results, whats plotted
#' @param window window in KB to plot in
#' @param geno.bed path to bed file, no .bed extension, character string
#' @param plot.r2.thresh minimum LD to plot
#' @param file.prefix optional label for output files
#' @param annotation.table table with gene descriptions and start/stop, currently only supports a specific table from maize
#' @param gene.highlight name of column that indicates which genes to highlight, column should be boolean, T is displayed red
#'
#' @return outputs 1 plot per row of top.hits table
#'
#' @examples
#' # no example
make_manhattan_zoom_annotation <- function(out.dir = "./",
                                           plot.type = c("png", "pdf"),
                                           trait.name,
                                           top.hits,
                                           gwas.res,
                                           window = 1000,
                                           geno.bed,
                                           plot.r2.thresh,
                                           file.prefix = "",
                                           annotation.table,
                                           gene.highlight){
  if(nrow(top.hits) > 0){
    if(plot.type == "pdf"){
      pdf(paste0(out.dir, trait.name, "_GwasGeneAnnotation.pdf"),
          height = 10,
          width = 15)
    }

    for(i in 1:nrow(top.hits)){
      this.chr <-top.hits$CHROM[i]
      this.pos <- top.hits$POS[i]
      this.snp.name <- paste0(this.chr, "-", this.pos)

      print(paste("Starting", trait.name, this.snp.name))

      # get info for hit from farmcpu output
      this.top.hit <- filter(top.hits, .data$SNP == this.snp.name)

      luebbert::make_ld(this.snp.name, window, geno.bed)

      ld.table <- read.table("./ld_out_temp.ld", header = T)
      ld.table_sub <- ld.table %>%
        select("SNP" = "SNP_B", "R2")

      # if plot window is different than the window that you calculated the LD for
      gwas.sub <- gwas.res %>%
        filter(.data$CHROM == this.chr) %>%
        filter(between(.data$POS, this.pos - window * 1000, this.pos + window * 1000)) %>%
        filter(.data$SNP != this.top.hit$SNP) %>%
        left_join(ld.table_sub, by = "SNP") %>%
        mutate(logpval = -log10(.data[[paste0(trait.name, ".GLM")]]))

      # make manhattan
      plot.df <- gwas.sub %>%
        # alpha scale
        mutate(how.to.plot = case_when(R2 > plot.r2.thresh ~ 1,
                                       TRUE ~ .4)) %>%
        # color scale
        mutate(plot.R2 = case_when(R2 < plot.r2.thresh ~ NA,
                                   TRUE ~ R2))

      # how far to spread labels past ends, in percentage
      y.spread.expansion <- .1
      y.spread.factors <- c(1 + y.spread.expansion, 1 - y.spread.expansion)
      y.spread.factor.window <- (window * y.spread.expansion) * 1000

      # plot limits
      plot.limits <- c(this.pos + window * 1000, this.pos - window * 1000)
      plot.limits.ex <- c(plot.limits[1] + y.spread.factor.window, plot.limits[2] - y.spread.factor.window)

      # my.colors <- viridis::viridis_pal(option = "turbo")(5)
      my.colors <- viridis::viridis_pal(option = "plasma")(5)

      man <-
        ggplot(aes(x = .data$POS, y = .data$logpval), data = plot.df) +
        geom_point(aes(color = .data$plot.R2, alpha = .data$how.to.plot), shape = 16) +
        scale_alpha(guide = "none") +
        geom_point(aes(x = .data$POS, y = -log10(.data[["PVAL"]])),
                   data = this.top.hit,
                   fill = "red", size = 4, shape = 23) +
        scale_color_stepsn(colors = my.colors, name = "R2") +
        theme(panel.grid = element_blank()) +
        scale_x_reverse(limits = plot.limits.ex,
                        labels = label_number(scale_cut = cut_short_scale()),
                        name = "-log10pval") +
        coord_flip() +
        scale_y_reverse() +
        theme(panel.background = element_rect(fill = "grey95"),
              axis.title.y = element_blank(),
              legend.position = "left") +
        labs(y = bquote(-log[10](p-value))) +
        geom_rug(sides = 't')

      # make smoothed R2 plot
      ldsmooth <-
        ggplot(aes(x = .data$POS, y = .data$R2), data = gwas.sub) +
        geom_smooth(se = F, color = 'black') +
        coord_flip() +
        scale_x_reverse(limits = plot.limits.ex) +
        scale_y_reverse() +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              panel.grid = element_blank())



      # make gene plot
      if(missing(gene.highlight)){
        x <- annotation.table %>%
          select("GeneID", "CHROM", "start", "end",
                 "Description", "GeneNameShort")
      } else {
        x <- annotation.table %>%
          select("GeneID", "CHROM", "start", "end",
                 "Description", "GeneNameShort", .data[[gene.highlight]])
      }

      anno.sub <- x %>%
        mutate(plot_name = case_when(is.na(.data$GeneNameShort) ~ .data$Description,
                                     TRUE ~ paste0(.data$GeneNameShort, ", ", .data$Description))) %>%
        filter(.data$CHROM == this.chr) %>%
        rowwise() %>%
        mutate(dist.from.snp = luebbert::get_gene_dist_from_snp(this.pos, .data$start, .data$end)) %>%
        filter(.data$dist.from.snp <= window * 1000)

      anno.spread <- anno.sub %>%
        add_column(y.pos = luebbert::get_gene_y_pos(min(plot.limits.ex),
                                                    max(plot.limits.ex),
                                                    nrow(anno.sub)))

      mid <- this.top.hit$POS
      breaks.anno <- seq(from = mid - window * 1000,
                         to = mid + window * 1000,
                         length.out = 9)

      if(nrow(anno.spread) > 100){
        text.min.size <- 3
      } else {
        text.min.size <- 4
      }

      anno.labels.fun <- function(x){
        paste0((x - this.top.hit$POS) / 1000, " KB")
      }

      anno <-
        ggplot(data = anno.spread) +
        geom_segment(aes(y = .data$start, yend = .data$end, x = .5, xend = .5), linewidth = 2, color = "red") +
        scale_y_reverse(limits = plot.limits.ex,
                        labels = anno.labels.fun,
                        breaks = breaks.anno) +
        xlim(.5, .875) +
        # geom_text(aes(y = y.pos, x = .5, label = paste0(" ", plot_name)),
        #           nudge_x = .05,
        #           hjust = 0,
        #           size = 2) +
        # ggfittext::geom_fit_text(aes(xmin = .55, xmax = .85, y = .data$y.pos, label = paste0(.data$plot_name)),
        #                          place = "left",
        #                          #grow = TRUE,
        #                          hjust = 0,
        #                          padding.y = grid::unit(.1, "lines"),
        #                          min.size = text.min.size) +
        geom_segment(aes(x = .5, xend = .55, y = .data$start, yend = .data$y.pos)) +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              #axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              panel.grid = element_blank())

      if(missing(gene.highlight)){
        anno <- anno +
          ggfittext::geom_fit_text(aes(xmin = .55, xmax = .85, y = .data$y.pos, label = paste0(.data$plot_name)),
                                   place = "left",
                                   #grow = TRUE,
                                   hjust = 0,
                                   padding.y = grid::unit(.1, "lines"),
                                   min.size = text.min.size)
      } else {
        anno <- anno +
          ggfittext::geom_fit_text(aes(xmin = .55, xmax = .85, y = .data$y.pos, label = paste0(.data$plot_name),
                                       color = .data[[gene.highlight]]),
                                   place = "left",
                                   #grow = TRUE,
                                   hjust = 0,
                                   padding.y = grid::unit(.1, "lines"),
                                   min.size = text.min.size) +
          scale_color_manual(guide="none", values = c("black", "red"))
      }

      out.plot <- man + ldsmooth + anno +
        patchwork::plot_layout(nrow = 1, widths = c(3,1,6)) +
        patchwork::plot_annotation(title = paste(this.snp.name, trait.name, sep = ", "))

      if(plot.type == "png"){
        memesave(paste0(out.dir, "/", file.prefix, trait.name, "_", this.snp.name, "_GwasGeneAnnotation.png"),
                 plot = out.plot,
                 height = 10 * 1.2, width = 15 * 1.2)
      }

      if(plot.type == "pdf"){
        print(out.plot)
      }
    }


    if(plot.type == "pdf"){
      dev.off()
    }
  } else {
    write.table("NO SIGNALS", file = paste0(out.dir, "/NO_SIGNIFICANT_SIGNALS_for_", trait.name, ".txt"))
  }
  # clean up temp ld file
  system("rm ./ld_out_temp*")
}
