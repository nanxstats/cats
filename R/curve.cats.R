#' Curves plot for cats object
#'
#' @param RR description
#' @param MAF description
#' @param rep description
#' @param SNPs description
#' @param ncases description
#' @param ncontrols description
#' @param ncases2 description
#' @param ncontrols2 description
#' @param alpha description
#' @param ... description
#'
#' @export curve.cats
curve.cats <- function(
  RR, MAF, rep = 1536, SNPs = 1E6,
  ncases, ncontrols, ncases2, ncontrols2,
  alpha = 0.05 / SNPs, ...
) {
  powerList.O <- c()
  powerList.J <- c()
  powerList.R <- c()
  powerList.F <- c()
  power.O <- rep(0, length(RR))
  power.F <- rep(0, length(RR))
  power.J <- rep(0, length(RR))
  power.R <- rep(0, length(RR))
  for (nmaf in 1:length(MAF)) {
    for (tal in 1:length(RR)) {
      if (power.F[tal] > 0.995 & power.R[tal] > 0.995) {
        power.O[tal] <- 1
        power.R[tal] <- 1
        power.J[tal] <- 1
        power.F[tal] <- 1
        break
      }
      temp <- cats(risk = RR[tal], freq = MAF[nmaf], ncases = ncases, ncontrols = ncontrols, ncases2 = ncases2, ncontrols2 = ncontrols2, alpha = alpha, pimarkers = rep / SNPs, ...)
      power.O[tal] <- temp$P.one.study
      power.J[tal] <- temp$P.joint
      power.R[tal] <- temp$P.rep.study
      power.F[tal] <- temp$P.first.stage
    }
    powerList.O <- cbind(powerList.O, power.O)
    powerList.J <- cbind(powerList.J, power.J)
    powerList.R <- cbind(powerList.R, power.R)
    powerList.F <- cbind(powerList.F, power.F)
    cat(nmaf, " ")
  }
  cat("\n")

  obs <- list(powerList.O = powerList.O, powerList.J = powerList.J, powerList.R = powerList.R, powerList.F = powerList.F, RR = RR, MAF = MAF, ncases = ncases, ncontrols = ncontrols, ncases2 = ncases2, ncontrols2 = ncontrols2, rep = rep, curve = T)

  class(obs) <- "supercats"
  return(obs)
}
