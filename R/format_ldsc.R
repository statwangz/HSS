#' Format the LD score data.
#'
#' @param file_ldsc name of the LD score file.
#'
#' @return formatted LD score data
#' @export
#'
format_ldsc <- function(file_ldsc) {
  message("Begin reading in LD scores data... ")

  if (is.null(file_ldsc)) {
    stop("please provide the information on LD scores!")
  }

  ldsc_path <- paste0(file_ldsc, "/", c(1:22), ".l2.ldscore.gz")

  ldsc <- suppressMessages(readr::read_delim(ldsc_path[1], "\t",
    escape_double = F, trim_ws = T, progress = F
  ))
  for (i in 2:22) {
    ldsc <- bind_rows(ldsc, suppressMessages(readr::read_delim(ldsc_path[i], "\t",
      escape_double = F, trim_ws = T, progress = F
    )))
  }

  message("The LD score data set has ", nrow(ldsc), " lines.")

  return(ldsc)
}
