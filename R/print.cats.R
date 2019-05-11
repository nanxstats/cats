#' Print the calculated power
#'
#' Print the calculated power
#'
#' @param x a cats object
#' @param ... Dynamic parameter for the values of additional parameters for the
#' summary method
#'
#' @method print cats
#'
#' @seealso cats
#'
#' @export
print.cats <-
  function(x, ...) {
    if (!inherits(x, "cats")) {
      stop("Not an object of class cats!")
    }
    cat("Expected Power is;\n")
    print(structure(list(
      "For a one-stage study" = signif(x$P.one.study, 3),
      "For first stage in two-stage study" = signif(x$P.first.stage, 3),
      "For second stage in replication analysis" = signif(x$P.rep.study, 3),
      "For second stage in a joint analysis" = signif(x$P.joint, 3),
      "pi" = signif(x$pi, 3)
    ),
    class = "power.htest"
    ))
  }
