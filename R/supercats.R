#' Generate supercats
#'
#' @param RR description
#' @param MAFmax description
#' @param MAFmin description
#' @param by description
#' @param rep description
#' @param SNPs description
#' @param ncases description
#' @param ncontrols description
#' @param ncases2 description
#' @param ncontrols2 description
#' @param alpha description
#' @param ... description
#'
#' @export super.cats

# functions for plots
super.cats <- function(RR, MAFmax = 0.5, MAFmin = 0.005, by = 50, rep = 1536, SNPs = 1E6, ncases, ncontrols, ncases2, ncontrols2, alpha = 0.05 / SNPs, ...) {
  powerList.O <- c()
  powerList.J <- c()
  powerList.R <- c()
  powerList.F <- c()
  power.O <- rep(0, length(RR))
  power.F <- rep(0, length(RR))
  power.J <- rep(0, length(RR))
  power.R <- rep(0, length(RR))

  MAF <- exp(seq(log(MAFmin), log(MAFmax), by = (log(MAFmax) - log(MAFmin)) / by))
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

  obs <- list(powerList.O = powerList.O, powerList.J = powerList.J, powerList.R = powerList.R, powerList.F = powerList.F, RR = RR, MAF = MAF, ncases = ncases, ncontrols = ncontrols, ncases2 = ncases2, ncontrols2 = ncontrols2, rep = rep, curve = F)

  class(obs) <- "supercats"
  return(obs)
}

if (FALSE) {
  # heat plot
  rr <- seq(1, 2, by = 0.025)
  c <- super.cats(rr, by = length(rr), ncases = 765, ncontrols = 1274, ncases2 = 100, ncontrols2 = 100, alpha = 0.001, prevalence = 0.01)
  plot(c, main = "power", file = NULL)

  # curves
  rr <- seq(1, 3, by = 0.05)
  maf <- c(0.01, 0.05, 0.2, 0.5)
  c2 <- curve.cats(rr, maf, ncases = 765, ncontrols = 1274, ncases2 = 100, ncontrols2 = 100, alpha = 0.001, prevalence = 0.01)
  plot(c2, main = "power2", ylab = "Power", xlab = "RR", file = NULL, col = 1:4)
}
