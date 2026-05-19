#' Logit Function
#'
#' @description
#' Helper functions.
#'
#' @param p numeric in (0,1).
#'
#' @export
logit <- function(p){
    log(p/(1-p))
}
#' Expit Function
#'
#'
#' @param x numeric
#' @export
expit <- function(x){
    1/(1 + exp(-x))
}
#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`
