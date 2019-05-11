#' Plot supercats
#'
#' @param x description
#' @param type description
#' @param file description
#' @param col description
#' @param main description
#' @param ... description
#'
#' @method plot supercats
#'
#' @export
plot.supercats <- function(x, type = "Joint", file = "power.pdf", col = NULL, main = paste("POWER N=", x$ncases, ":", x$ncontrols, ",", x$ncases2, ":", x$ncontrols2, " rep=", x$rep, sep = ""), ...) {
  if (type == "Joint") {
    power <- x$powerList.J
  } else if (type == "One") {
    power <- x$powerList.O
  } else if (type == "Replication") {
    power <- x$powerList.R
  } else if (type == "First") {
    power <- x$powerList.F
  }

  if (!is.null(file)) {
    pdf(file)
  }
  # curve
  if (x$curve) {
    if (is.null(col)) {
      col <- 1:length(x$MAF)
    }
    plot(x$RR, power[, 1], ylim = c(0, 1), main = main, col = "transparent", ...)
    for (nmaf in 1:length(x$MAF)) {
      lines(x$RR, power[, nmaf], col = col[nmaf], lwd = 2)
    }
    legend(min(x$RR), 1, paste("MAF=", x$MAF), col = col, lwd = 2, bty = "n")
  }
  else {
    # image
    if (is.null(col)) {
      col <- heat.colors(80)
    }

    image(x$RR, x$MAF, power, col = col, main = main, log = "y", ylim = c(0.005, .5), ylab = "MAF", xlab = "RR", ...)
    legend("topright", paste(1:10 * 10, "%"), fill = col[1:10 * 8], bty = "n")
  }
  if (!is.null(file)) {
    dev.off()
  }
}
