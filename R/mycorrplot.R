

#' Make visuzlization of correlation matrix
#'
#' @param dat dataframe of data to compute correlations, columns are traits
#' @param use how to deal with missing data in correlation, default is pairwise complete
#' @param axis.text.size size of trait labels in plot
#' @param coef.text.size size of correlation coefficient labels in plot
#' @param color.scheme which color scheme to use, diverging schemes from khroma package
#' @param output.png output corrplot as png?
#' @param out.file.png where to output png file
#'
#' @return prints visualization of correlation plot
#' @export
#'
#' @examples
mycorrplot <- function(dat,
                       use = "pairwise.complete",
                       axis.text.size = 1,
                       coef.text.size = 1,
                       color.scheme = c("sunset", "BuRd", "nightfall", "PRGn"),
                       output.png = F,
                       out.file.png = NULL,
                       pvalthresh = .05){

  color.scheme <- match.arg(color.scheme)
  corr.colors <- khroma::color(color.scheme)
  m <- cor(dat, use = use)

  p.mat <- corrplot::cor.mtest(dat, use = use)$p

  if(output.png){
    png(out.file.png,
        units = "in",
        height = 8, width = 8,
        res = 100)
    corrplot::corrplot(m,
                       method="color",
                       col=corr.colors(8),
                       type="full",
                       order="original",
                       addCoef.col = "black", # Add coefficient of correlation
                       tl.col="black", tl.srt=45, #Text label color and rotation
                       # Combine with significance
                       p.mat = p.mat, sig.level = pvalthresh, insig = "blank",
                       # hide correlation coefficient on the principal diagonal
                       diag=FALSE,
                       tl.cex = axis.text.size,
                       number.cex = coef.text.size,
                       bg = "light grey",
                       addgrid.col = "black")
    dev.off()
  } else {
    corrplot::corrplot(m,
                       method="color",
                       col=corr.colors(8),
                       type="full",
                       order="original",
                       addCoef.col = "black", # Add coefficient of correlation
                       tl.col="black", tl.srt=45, #Text label color and rotation
                       # Combine with significance
                       p.mat = p.mat, sig.level = pvalthresh, insig = "blank",
                       # hide correlation coefficient on the principal diagonal
                       diag=FALSE,
                       tl.cex = axis.text.size,
                       number.cex = coef.text.size,
                       bg = "light grey",
                       addgrid.col = "black")
  }
}

