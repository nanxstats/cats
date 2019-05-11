#' Power calculation for a joint analysis of a two-stage case control design
#' for SNP data
#'
#' Power calculation for a joint analysis of a two-stage case control design
#' for SNP data.
#'
#' These power analysis are based on Skol et al. 2006, But are generized so
#' that the ratio between cases and controls may vary between stages. Also the
#' allele frequencies, disease prevalence and relative risk are also allowed to
#' vary. The joint statistic $z_joint=z_1\\sqrt\\pi+z_2\\sqrt1-\\pi$ where $z_1$ is
#' the z-score for the first stage and the weight $\\pi$ is calculated as
#' $\\pi=1/var(\\hatp'_1-\\hatp_1)*(1/var(\\hatp'_1-\\hatp_1)+1/var(\\hatp'_2-\\hatp_2))^-1$,
#' where $\\hatp'_1$ is the estimate of the allele frequency of the cases in the
#' first stage. This is consistent with Skol et al 2006 when the ratios of
#' cases and controls are the same in both stages. When this is not the case
#' the weight $\\pi$ may vary slightly with different allele frequencies and
#' different relative risks. For power calculations I would recommend
#' calculating the weight at a likely scenario where there is about 80-90\%
#' power and fixing the weights at other scenarios (and the testing of the real
#' data) to this weight. This can be done by assigning pisample to a value. In
#' practice this will hardly affect the power.
#'
#' @param freq numeric. The minor allele frequency (MAF) in the first stage
#' @param freq2 numeric. The MAF in the second stage, Optional, if -1 the same
#' value as for the first stage is given
#' @param ncases integer. The number of cases in the first stage
#' @param ncontrols integer. The number of controls in the first stage
#' @param ncases2 integer. The number of cases in the second stage
#' @param ncontrols2 integer. The number of controls in the second stage
#' @param risk numeric. The relative risk in the first stage
#' @param risk2 numeric. The relative risk in the second stage, Optional, if -1
#' the same value as for the first stage is given
#' @param pisamples numeric. The weights used for the joint statistic.
#' Optional. see details
#' @param prevalence numeric. The prevalence of the disease in the population
#' for the first stage
#' @param prevalence2 numeric. The prevalence of the disease in the population
#' for the second stag, Optional, if -1 the same value as for the first stage
#' is given
#' @param additive boolean. if 1 an additive model is assumed
#' @param recessive boolean. if 1 a recessive model is assumed
#' @param dominant boolean. if 1 a dominant model is assumed
#' @param multiplicative boolean. if 1 a multiplicative model is assumed
#' @param alpha numeric. The significance threshold. Often the a threshold of
#' 0.05 divided by the number of markers is chosen
#' @param pimarkers numeric. The fraction of markers genotyped in the second
#' stage
#'
#' @return \item{P.one.study}{The power if only one study was performed, NB!
#' This is only a valid estimate if the relative risk and allele frequency is
#' the same for both stages} \item{P.first.stage}{The power for a marker to
#' proceed the the second stage} \item{P.rep.study }{The power of the study if
#' based on replication and not a joint analysis} \item{P.joint.min }{The power
#' of the joint analysis tp detect at least one susceptibility SNP assuming
#' that five susceptibility SNPs exits} \item{P.joint }{The power of the joint
#' analysis} \item{pi }{The weight used to calculate the joint statistic}
#' \item{T.one.study}{Recommended thresholds for a one-stage study}
#' \item{T.first.stage }{Recommended thresholds for the first stage in
#' two-stage study} \item{T.second.stage.rep }{Recommended thresholds for the
#' second stage in replication analysis} \item{T.second.stage.joint
#' }{Recommended thresholds for the second stage in a joint analysis}
#' \item{E.Disease.freq.cases1 }{The expected disease allele frequency in stage
#' 1 for cases} \item{E.Disease.freq.controls1 }{The expected disease allele
#' frequency in stage 1 for controls} \item{E.Disease.freq.cases2 }{The
#' expected disease allele frequency in stage 2 for cases}
#' \item{E.Disease.freq.controls2 }{The expected disease allele frequency in
#' stage 2 for controls}
#'
#' @author Anders Albrechtsen
#'
#' @seealso \url{http://www.sph.umich.edu/csg/abecasis/CaTS/}
#'
#' @references Skol AD, Scott LJ, Abecasis GR, Boehnke M: Joint analysis is
#' more efficient than replication-based analysis for two-stage genome-wide
#' association studies. Nat Genet 38: 209-213, 2006.
#'
#' @keywords htest
#'
#' @examples
#' # calculate the power under a multiplicative model using a two stage design
#' # and assuming a relative risk of 1.5
#' cats(
#'   freq = 0.2,
#'   ncases = 500, ncases2 = 500,
#'   ncontrols = 1000, ncontrols2 = 1000,
#'   risk = 1.5, multiplicative = 1
#' )
#'
#' power.J <- c()
#' power.R <- c()
#' power.O <- c()
#' RR <- 23:32 / 20
#' for (tal in 1:length(RR)) {
#'   temp <- cats(risk = RR[tal])
#'   power.J[tal] <- temp$P.joint
#'   power.R[tal] <- temp$P.rep.study
#'   power.O[tal] <- temp$P.one.study
#' }
#' plot(RR, power.J, type = "b", lwd = 2, ylab = "Power")
#' lines(RR, power.R, lwd = 2, col = 2, type = "b")
#' lines(RR, power.O, lwd = 2, col = 3, type = "b")
#' legend(1.4, 0.4, c(
#'   "joint analysis", "replication design",
#'   "one stage design"
#' ), col = 1:3, lwd = 3, bty = "n")
#' @export cats
cats <-
  function(
             freq = 0.5, freq2 = -1,
             ncases = 500, ncontrols = 500,
             ncases2 = 500, ncontrols2 = 500,
             risk = 1.5, risk2 = -1,
             pisamples = -1,
             prevalence = 0.1, prevalence2 = -1,
             additive = 0, recessive = 0, dominant = 0, multiplicative = 1,
             alpha = 0.0000001,
             pimarkers = 0.00316) {
    model <- c(additive, recessive, dominant, multiplicative)

    if (sum(model == 1) != 1) {
      stop("chose only one model. e.i. one model must be 1 the others 0")
    }
    if (sum(model == 0) != 3) {
      stop("chose only one model. e.i. one model must be 1 the others 0")
    }


    if (freq < 0 | freq > 1) {
      stop("freq must be between 0 and 1")
    }
    if ((freq2 < 0 | freq2 > 1) & freq2 != -1) {
      stop("freq2 must be between 0 and 1 (or undefined as -1)")
    }
    if ((pisamples < 0 | pisamples > 1) & pisamples != -1) {
      stop("pisamples must be between 0 and 1")
    }
    if ((prevalence2 < 0 | prevalence2 > 1) & prevalence2 != -1) {
      stop("prevalence2 must be between 0 and 1 (or undefined as -1)")
    }
    if (alpha < 0 | alpha > 1) {
      stop("alpha must be between 0 and 1")
    }
    if (prevalence < 0 | prevalence > 1) {
      stop("prevalence must be between 0 and 1")
    }
    if (pimarkers < 0 | pimarkers > 1) {
      stop("pimarkers must be between 0 and 1")
    }
    if (ncases != as.integer(ncases) | ncases < 0) {
      stop("ncases must be a positive integer")
    }
    if (ncases2 != as.integer(ncases2) | ncases2 < 0) {
      stop("ncases2 must be a positive integer")
    }
    if (ncontrols != as.integer(ncontrols) | ncontrols < 0) {
      stop("ncontrols must be a positive integer")
    }
    if (ncontrols2 != as.integer(ncontrols2) | ncontrols2 < 0) {
      stop("ncontrols2 must be a positive integer")
    }
    if (risk < 0) {
      stop("risk must be positive")
    }
    if (risk2 < 0 & risk2 != -1) {
      stop("risk2 must be positive(or undefined as -1)")
    }

    res <- cats_c(
      freq, freq2, ncases, ncontrols,
      ncases2, ncontrols2, risk, risk2,
      pisamples, prevalence, prevalence2,
      additive, recessive, dominant,
      multiplicative, alpha, pimarkers
    )

    options <- cbind(freq, freq2, ncases = ncases, ncontrols = ncontrols, ncases2 = ncases2, ncontrols2 = ncontrols2, risk, risk2, pisamples, prevalence, prevalence2, additive, recessive, dominant, multiplicative, alpha, pimarkers)

    result <- list(P.one.study = res[1, 1], P.first.stage = res[2, 1], P.rep.study = res[3, 1], P.joint.min = res[4, 1], P.joint = res[5, 1], pi = res[6, 1], T.one.study = res[7, 1], T.first.stage = res[8, 1], T.second.stage.rep = res[9, 1], T.second.stage.joint = res[10, 1], E.Disease.freq.cases1 = res[11, 1], E.Disease.freq.controls1 = res[12, 1], E.Disease.freq.cases2 = res[13, 1], E.Disease.freq.controls2 = res[14, 1], options = options)
    class(result) <- "cats"
    return(result)
  }
