#' Plot PC contributions
#'
#' @param pca_obj pca object created from "princomp()"
#' @param plot.title title for plot
#' @param num.pcs number pcs to plot either as maximum or as vector of desired pc's
#' @param label.offset percentage to offset percent contribution label from side of plot
#'
#' @return GGplot of variable contributions for each PC
#' @export
#'
#' @examples
#' #pcobj <- princomp(data)
#' #plot_pc_contribs(pcobj)
plot_pc_contribs <- function(pca_obj, plot.title = "", num.pcs = 6, label.offset){

  if(length(num.pcs) == 1){
    get.these.pcs <- c(1:num.pcs)
  } else {
    get.these.pcs <- num.pcs
  }

  var.contrib <- factoextra::get_pca_var(pca_obj)$contrib
  var.contrib.long <- reshape2::melt(var.contrib)

  component.contribs <- round(pca_obj$sdev^2 / sum(pca_obj$sdev^2) * 100, 2)
  names(component.contribs) <- stringr::str_replace(names(component.contribs), "Comp", "Dim")
  component.contribs <- data.frame(Var2 = names(component.contribs), contrib_pct = component.contribs)
  var.contrib.long <- dplyr::left_join(var.contrib.long, component.contribs, by = "Var2") %>%
    rename("Var" = "Var1", "PC.dim" = "Var2") %>%
    mutate(pc.num = str_extract(.data$PC.dim, "\\d+")) %>%
    filter(.data$pc.num %in% get.these.pcs)

  contrib.pcts <- var.contrib.long %>%
    select("PC.dim", "contrib_pct") %>%
    distinct()

  label.y <- max(var.contrib.long$value) * label.offset
  label.x <- length(unique(var.contrib.long$Var)) / 2

  p <-
    ggplot2::ggplot(var.contrib.long,
                    ggplot2::aes(x = .data$Var, y = .data$value, fill = .data$Var)) +
    geom_bar(stat = "identity") +
    geom_label(aes(x = length(unique(.data$Var)) / 2,
                   y = label.y,
                   label = .data$contrib_pct),
               fill = "white",
               data = contrib.pcts) +
    facet_wrap(ggplot2::vars(.data$PC.dim)) +
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
