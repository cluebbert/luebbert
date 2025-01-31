#' Make table of offsets for plotting manhattans
#'
#' @param df dataframe of gwas results with minimally chromosome and position columns
#' @param chrom.col character - name of chromosome column
#' @param pos.col character - name of position column
#'
#' @return table of offsets to add for each chromosome to get genome in continuous line
#' @export
#'
#' @examples
make_manhattan_offsets <- function(df, chrom.col, pos.col){
  df.sub <- df %>%
    select(one_of(chrom.col, pos.col)) %>%
    setNames(c("chrom.col", "pos.col"))

  offset.table <- df.sub %>%
    group_by(.data$chrom.col) %>%
    mutate(max.pos = max(.data$pos.col)) %>%
    select(-"pos.col") %>%
    distinct() %>%
    arrange(.data$chrom.col) %>%
    ungroup() %>%
    mutate(lag_max = dplyr::lag(.data$max.pos, 1)) %>% # need to account for the fact that chrom 1 does not have offset
    mutate(lag_max = case_when(is.na(.data$lag_max) ~ 0,
                               TRUE ~ as.numeric(.data$lag_max))) %>%
    mutate(offset = cumsum(.data$lag_max)) %>%
    select(-"max.pos", -"lag_max")

  names(offset.table)[1] <- "chromosome"

  return(offset.table)
}
