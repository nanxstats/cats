#' Power calculation for a joint analysis of a two-stage case control design
#' for SNP data
#'
#' Power calculation for a joint analysis of a two-stage case control design
#' for SNP data.
#'
#' @name cats-package
#'
#' @docType package
#'
#' @author Anders Albrechtsen <albrecht@@binf.ku.dk>
#'
#' @seealso \url{http://www.sph.umich.edu/csg/abecasis/CaTS/}
#'
#' @references Skol AD, Scott LJ, Abecasis GR, Boehnke M: Joint analysis is
#' more efficient than replication-based analysis for two-stage genome-wide
#' association studies. Nat Genet 38: 209-213, 2006.
#'
#' @importFrom grDevices dev.off heat.colors pdf
#' @importFrom graphics image legend lines plot
#'
#' @examples
#' # calculate the power under a multiplicative model using a
#' # two stage design and assuming a relative risk of 1.5
#' cats(
#'   freq = 0.2,
#'   ncases = 500, ncases2 = 500,
#'   ncontrols = 1000, ncontrols2 = 1000,
#'   risk = 1.5, multiplicative = 1
#' )
NULL
