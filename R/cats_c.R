#' @useDynLib cats,.registration = TRUE
#' @useDynLib cats cats_
cats_c <- function(
  freq, freq2, ncases, ncontrols,
  ncases2, ncontrols2, risk, risk2,
  pisamples, prevalence, prevalence2,
  additive, recessive, dominant,
  multiplicative, alpha, pimarkers) {
  .Call(
    cats_,
    as.double(freq), as.double(freq2), as.integer(ncases), as.integer(ncontrols),
    as.integer(ncases2), as.integer(ncontrols2), as.double(risk), as.double(risk2),
    as.double(pisamples), as.double(prevalence), as.double(prevalence2),
    as.integer(additive), as.integer(recessive), as.integer(dominant),
    as.integer(multiplicative), as.double(alpha), as.double(pimarkers),
    PACKAGE = "cats"
  )
}
