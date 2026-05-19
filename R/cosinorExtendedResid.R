#' cosinorExtendedResid
#'
#' @description
#' Get residual values. For use in optimization routine for model fit.
#'
#' @param par_optim vector of length corresponding to transformed parameter values.
#' See [cosinorExtendedTransformParameters].
#' @param phi_cos a centering reference value for the transformation of phi.
#' @param Y observed data timeseries.
#' @param t time (in hours since midnight) corresponding to values of Y.
#'
#' @returns
#' A vector of residual values for extended cosinor model.
#'
#' @export
cosinorExtendedResid = function(par_optim, phi_cos, Y, t){
    par_ext = cosinorExtendedTransformParameters(par_optim, phi_cos)
    rt = cosinorExtendedFitted(par_ext, t)
    Y - rt
}
