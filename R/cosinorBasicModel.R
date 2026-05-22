#' cosinorBasicModel
#'
#' @description
#' Fit the three parameter cosinor model.
#' This function is a modified from the function ActExtendCosinor in ActCR.
#' Updates include performing the optimization on a transformed version of the
#' parameter space to avoid issues with boundaries.
#'
#' @details The cosinor model is:
#' \deqn{y(t) = \mesor_{\mathrm{cos}} + \mathrm{gam}_{\mathrm{cos}}\cdot \cos(\left[ t - \phi_{\mathrm{ext}} \right]\cdot 2\pi/24) + \epsilon(t)}
#'
#' @param x input timeseries with missing epochs represented as "NA" values. Note
#' \code{length(x)*window/1440} must be an integer.
#' @param window epoch duration in minutes.
#' @param export_ts Should the original data and fitted values be exported?
#'
#' @returns
#' A list with up to two elements:
#' \itemize{
#'  \item \code{estimates}: a one-row tibble of cosinor parameters and model fit summaries.
#'      \itemize{
#'          \item The estimates of the three primary cosinor parameters \code{mesor_cos}, \code{gam_cos}, \code{phi_cos}.
#'          \item Secondary estimates from the extended cosinor model fit including the amplitude and model resuduals (\code{rss_cos} and model R^2 (\code{R2_cos})).
#'      }
#'  \item \code{cosinor_ts}: a tibble containing epoch-level cosinor model fits and observed data (see \code{export_ts}).
#' }
#'
#' @examples
#' t <- seq(0, 1440*7,by = 1)[-1]/60
#' x <- 2 + 3*cos((t- 5)/24*2*pi) + rnorm(length(t),1)
#' cosOut <- cosinorBasicModel(x = x, window = 1,export_ts = TRUE)
#' cosOut$estimates
#' plot(
#'     data = cosOut$cosinor_ts,
#'     y~time_across_days
#' )
#' lines(
#'     data = cosOut$cosinor_ts,
#'     y_cos~time_across_days, col = "red"
#' )
#'
#'
#' @export
cosinorBasicModel <- function(
    x,
    window = 1,
    export_ts = FALSE
){
    t <- seq(window, length(x*window), by = window)/60
    cosOut = lm(
        x ~ I(cos(t*2*pi/24)) + I(sin(t*2*pi/24))
    )

    e_coef <- coef(cosOut)
    names(e_coef) <- NULL
    e_phi = atan2(e_coef[3], e_coef[2])*24/2/pi
    e_amp = e_coef[2]/cos(e_phi*2*pi/24)

    ret = list(
        estimates = tibble::tibble(
            mesor_cos = e_coef[1],
            phi_cos = e_phi,
            gam_cos = e_coef[2]/cos(e_phi*2*pi/24),
            amplitude_cos = 2*e_coef[2]/cos(e_phi*2*pi/24),
            rss_cos = sum((x - predict(cosOut))^2),
            tss_y = sum((x - mean(x))^2),
            R2_cos = (tss_y - rss_cos)/tss_y
        )
    )

    if(export_ts){
        ret$cosinor_ts = tibble::tibble(
            time_across_days = t,
            y = x,
            y_cos = predict(cosOut)
        )
    }
    return(ret)
}
