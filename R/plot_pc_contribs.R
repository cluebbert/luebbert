#' Plot PC contributions
#'
#' @param pca_obj pca object created from "princomp()"
#' @param plot.title title for plot
#'
#' @return GGplot of variable contributions for each PC
#' @export
#'
#' @examples
#' pcobj <- princomp(data)
#' plot_pc_contribs(pcobj)
plot_pc_contribs <- function(pca_obj, plot.title = ""){

  var.contrib <- factoextra::get_pca_var(pca_obj)$contrib
  var.contrib.long <- reshape2::melt(var.contrib)

  component.contribs <- round(pca_obj$sdev^2 / sum(pca_obj$sdev^2) * 100, 2)
  names(component.contribs) <- str_replace(names(component.contribs), "Comp", "Dim")
  component.contribs <- data.frame(Var2 = names(component.contribs), contrib_pct = component.contribs)
  var.contrib.long <- left_join(var.contrib.long, component.contribs, by = "Var2")

  p <-
    ggplot2::ggplot(filter(var.contrib.long, Var2 %in% c("Dim.1", "Dim.2", "Dim.3", "Dim.4", "Dim.5", "Dim.6")),
                    aes(x = Var1, y = value, fill = Var1)) +
    geom_bar(stat = "identity") +
    geom_label(aes(x = length(unique(Var1)) / 2, y = max(value) * .9, label = contrib_pct), fill = "white") +
    facet_wrap(~Var2) +
    theme_bw() +
    coord_flip() +
    theme(axis.text.x = element_text(hjust = 0),
          legend.position = "none",
          axis.title.y = element_blank()) +
    labs(title = plot.title,
         y = "Percent Contribution")
  if(length(unique(var.contrib.long$Var1)) > 30){
    p <- p + theme(axis.text.y = element_text(size = 3))
  }

  print(p)
}
