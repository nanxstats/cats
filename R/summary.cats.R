#' Summary of a cats object
#'
#' Summary of a cats object
#'
#' @param object a cats object
#' @param ... Dynamic parameter for the values of additional parameters
#' for the print method
#'
#' @method summary cats
#'
#' @seealso cats
#'
#' @export
#'
#' @examples
#' x <- cats()
#' summary(x)
summary.cats <-
  function(object, ...) {
    if (!inherits(object, "cats")) {
      stop("Not an object of class cats!")
    }
    cat("Options \n")
    ob <- t(object$options)
    colnames(ob) <- "chosen"
    print(ob)

    cat("Recommended thresholds:")
    print(structure(list(
      "One stage Design" = object$T.one.study,
      "Stage 1 Threshold" = object$T.first.stage,
      "Replication Threshold" = object$T.second.stage.rep,
      "Joint Analysis Threshold" = object$T.second.stage.joint
    ),
    class = "power.htest"
    ))

    cat("Eobjectpected disesase allele frequencies")
    print(structure(list(
      "Cases in stage 1" = object$E.Disease.freq.cases1,
      "Controls in stage 1 " = object$E.Disease.freq.controls1,
      "Cases in stage 2" = object$E.Disease.freq.cases2,
      "Controls in stage 2" = object$E.Disease.freq.controls2
    ),
    class = "power.htest"
    ))

    cat("Expected Power is:")
    print(structure(list("For a one-stage study" = signif(
      object$P.one.study,
      3
    ), "For first stage in two-stage study" = signif(
      object$P.first.stage,
      3
    ), "For second stage in replication analysis" = signif(
      object$P.rep.study,
      3
    ), "For second stage in a joint analysis" = signif(
      object$P.joint,
      3
    ), pi = signif(object$pi, 3)), class = "power.htest"))
  }
