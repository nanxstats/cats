#' Line plots for cats object
#'
#' @param x description
#' @param type description
#' @param col description
#' @param lty description
#' @param ... description
#'
#' @export lines.cats
lines.cats <- function(x, type = "Replication", col = NULL, lty = 2, ...) {
  if (type == "Joint") {
    power <- x$powerList.J
  } else if (type == "One") {
    power <- x$powerList.O
  } else if (type == "Replication") {
    power <- x$powerList.R
  } else if (type == "First") {
    power <- x$powerList.F
  }

  if (x$curve) {
    if (is.null(col)) {
      col <- 1:length(x$MAF)
    }
    for (nmaf in 1:length(x$MAF)) {
      lines(x$RR, power[, nmaf], col = col[nmaf], lwd = 2, lty = lty)
    }
  }
  else {
    cat("only for curves \n")
  }
}
