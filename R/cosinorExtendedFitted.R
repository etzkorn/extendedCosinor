#' cosinorExtendedFitted
#'
#' @description
#' Get fitted values from [extendedCosinor::cosinorExtendedModel].
#'
#' @param t time in hours since midnight.
#' @param par_ext named vector of length 5 with \code{mu_ext}, \code{gamma_ext}, \code{alpha_ext}, \code{beta_ext}, and \code{phi_ext}.
#'
#' @returns
#' A vector of fitted values for extended cosinor model.
#'
#' @export
cosinorExtendedFitted <- function(par_ext, t){
    ct = cos((t - par_ext["phi_ext"]) * 2 * pi / 24)
    lct = expit(par_ext["beta_ext"] * (ct - par_ext["alpha_ext"]))
    par_ext["mu_ext"] + par_ext["gamma_ext"] * lct
}
