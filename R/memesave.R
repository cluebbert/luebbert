#' GGsave for me! ggsave with some defaults
#'
#' @param myfilename filename to output
#' @param device device (can change this)
#' @param units defaults to inches
#' @param height some sensible and large default
#' @param width see height
#' @param bg no transparent background pngs
#'
#' @return write plot
#' @export
#'
#' @examples
memesave <- function(myfilename,
                     device = "png",
                     units = "in",
                     height = 10,
                     width = 15,
                     bg = "white"){
  if(device == "png"){
    ggplot2::ggsave(filename = myfilename,
                    device = "png",
                    dpi = "retina",
                    height = height,
                    width = width,
                    units = units,
                    bg = bg)
  } else {
    ggplot2::ggsave(filename = myfilename,
                    device = device,
                    height = height,
                    width = width,
                    units = units,
                    bg = bg)
  }
}
