#' cosinorExtendedTransformParameters
#'
#' @description
#' Transform parameters from those used in optimization routine to those specified in model (see [cosinorExtendedModel]).
#'
#' @param par_optim vector of length 5 corresponding to transformed parameter values.
#' @param phi_cos a centering reference value for the transformation of phi.
#'
#' @returns
#' Named vector of length 5 with \code{min_ext}, \code{amp_ext}, \code{alpha_ext}, \code{beta_ext}, and \code{phi_ext}.
#'
#' @export
cosinorExtendedTransformParameters <- function(par_optim, phi_cos){
    c(
        "min_ext" = par_optim[1],
        "amp_ext" = exp(par_optim[2]),
        "alpha_ext" = atan(par_optim[3])*2/pi,
        "beta_ext" = exp(par_optim[4]),
        "phi_ext" = (phi_cos + atan(par_optim[5])*24/2/pi)%%24
    )
}
