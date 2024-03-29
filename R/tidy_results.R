#' Prints results from moveSIM() or energySIM() in an easier-to-read table.
#' @import kableExtra
#' @param data The output from moveSIM() or energySIM() -- a list of two dataframes.
#' @param type "run_params" or "results", corresponding to which component of
#' your moveSIM() or energySIM() output you'd like to print out. Default
#' "results", which contains the movement data.
#' @param nrows The number of rows to print.
#' @details missing_pct and mortality_pct are not function paramet?ers,
#' but are nonetheless computed and returned here for your convenience
#'  \itemize{
#'  \item{missing_pct}{: What percent of rows exhibiting a missing
#'  location value (due to agent death or agent stopping)}
#'  \item{mortality_pct}{: What percent of simulated agents experienced death?}
#'  }
#' @return Prints a cleaned table of moveSIM() or energySIM() results.
#' @export

tidy_results <- function(data, type = "results", nrows = NULL) {
  if (type == "run_params") {
    datapar <- data$run_params
    my_params <- colnames(datapar)
    my_data <- as.list(datapar)
    my_data <- unlist(my_data)
    filler <- rep(NA, length(my_params))
    test2 <- as.data.frame(filler)
    test2[, 1] <- my_data
    rownames(test2) <- my_params
    colnames(test2) <- "Param value"
    my_output <- kable(test2, digits = 3, format = "html", row.names = TRUE) %>%
      kable_styling(
        bootstrap_options = c("bordered", "hover", "striped", "condensed"),
        full_width = F,
        font_size = 12,
        position = "left"
      )
    return(my_output)
  }
  else if (type == "results") {
    data.res <- data$results
    if (!is.null(nrows)) {
      data.res <- head(data.res, nrows)
    }
    my_output <- kable(data.res, digits = 3, format = "html", row.names = TRUE) %>%
      kable_styling(
        bootstrap_options = c("bordered", "hover", "striped", "condensed"),
        full_width = F,
        font_size = 12,
        position = "left"
      )
    return(my_output)
  }
}
