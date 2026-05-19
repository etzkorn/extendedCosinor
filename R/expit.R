#' Logit and Expit
#'
#' @description
#' Helper functions.
#'
#' @param x numeric.
#' @param p numeric in (0,1).
#'
#' @export
logit <- function(p){
    log(p/(1-p))
}
#' @export
expit <- function(x){
    1/(1 + exp(-x))
}
