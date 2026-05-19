#' cosinorExtendedModel
#'
#' @description
#' Fit the five parameter cosinor model from Marler et al. (2006).
#' This function is a modified from the function ActExtendCosinor in ActCR.
#' Updates include performing the optimization on a transformed version of the
#' parameter space to avoid issues with boundaries.
#'
#' @details The model specified by Marler et al uses the following equations:
#' \deqn{y(t) = \min_{\mathrm{ext}} + \mathrm{amp}_{\mathrm{ext}}\cdot \mathrm{expit}\bigg(\beta_{\mathrm{ext}} \Big(\cos(\left[t - \phi_{\mathrm{ext}}\right]\cdot 2\pi/24)-\alpha_{\mathrm{ext}} \Big) \bigg) + \epsilon(t)}
#'
#' @param x input timeseries with missing epochs represented as "NA" values. Note
#' \code{length(x)*window/1440} must be an integer.
#' @param window epoch duration in minutes.
#' @param export_ts Should the original data and fitted values be exported?
#' @param export_nls_outcome Should the message and convergence code be added to the output from
#' minpack.lm::nls.lm.
#' @param export_cosinor_param Should estimates from the basic three-parameter
#' cosinor model be added to the output?
#' @param maxiter Passed to [minpack.lm::nls.lm.control].
#' @param maxfev Passed to [minpack.lm::nls.lm.control].
#'
#' @returns
#' A list with up to two elements:
#' \itemize{
#'  \item \code{estimates}: a one-row tibble of extended cosinor parameters, derived estimates, model fit diagnostics (see \code{export_nls_outcome}), and basic cosinor parameters/estimates (see \code{export_cosinor_param}).
#'      \itemize{
#'          \item The estimates of the five primary extended cosinor parameters \code{min_ext}, \code{amp_ext}, \code{alpha_ext}, \code{beta_ext}, and \code{phi_ext} (see "details").
#'          \item Secondary estimates from the extended cosinor model fit including the maximum, minimum, range (amplitude), midpoint (mesor), time of activity rises above and falls below the mesor (\code{up_mesor_ext} and \code{down_mesor_ext}), total daily time the curve stays above the mesor (\code{hours_above_mesor_ext}), model resuduals (\code{rss_ext} and \code{F_pseudo_ext}).
#'          \item Convergence codes (\code{niter} and \code{info}) and messages (\code{message}) from minpack.lm::nls.lm.
#'          \item Parameters and secondary estimates from the three-parameter cosinor model: midpoint (\code{mesor_cos}), range (\code{amp_cos}), acrotime (\code{acrotime_cos}),  minimum(\code{minimum_cos}), and residual sum of squares (\code{rss_cos}).
#'          \item \code{ndays}: the number of unique days of data.
#'      }
#'  \item \code{cosinor_ts}: a tibble containing epoch-level cosinor and extended cosinor model fits and observed data (see \code{export_ts}).
#' }
#'
#' @references
#' Marler MR, Gehrman P, Martin JL, Ancoli‐Israel S.
#' The sigmoidally transformed cosine curve: a mathematical model for circadian
#' rhythms with symmetric non‐sinusoidal shapes. Statistics in medicine.
#' 2006 Nov 30;25(22):3893-904.
#'
#' @export
cosinorExtendedModel <- function(
        x,
        window = 1,
        export_ts = FALSE,
        export_cosinor_param = TRUE,
        export_nls_outcome = FALSE,
        maxiter = 1000,
        maxfev = 600
){
    dim = 1440/window
    ndays = length(x)/dim
    tmp.dat = data.frame(
        day = rep(1:ndays, each = dim),
        time = rep(1:dim, ndays)/(60/window),
        Y = x
    )
    cosOut = cosinor::cosinor.lm(
        Y ~ time(time) + 1,
        period = 24,
        data = tmp.dat[!is.na(tmp.dat$Y), ]
    )
    e_cos <- cosOut$coefficients
    names(e_cos) = NULL
    exCosOut = minpack.lm::nls.lm(
        c( # starting parameters
            max(e_cos[1] - e_cos[2], 0),
            log(2 * e_cos[2]),
            0,
            log(2),
            0
        ),
        fn = cosinorExtendedResid,
        theta0 = (-1) * cosinor2:::correct.acrophase(cosOut) * 24/(2 * pi),
        Y=tmp.dat$Y,
        t=tmp.dat$t,
        control = minpack.lm::nls.lm.control(
            maxiter = maxiter,
            maxfev = maxfev
        )
    )
    e_ext <- cosinorExtendedTransformParameters(
        exCosOut$par,
        phi_cos = (-1) * cosinor2:::correct.acrophase(cosOut) * 24/(2 * pi)
    )
    e_24 = cosinorExtendedFitted(e_ext,t = seq(0,24, length = 1441)[-1])
    estimates_out = tibble(
        t(e_ext)
    ) %>% mutate(
        ndays = ndays,
        acrotime_ext = e_ext["phi_ext"],
        maximum_ext = e_ext["min_ext"] + e_ext["amp_ext"] * expit(e_ext["beta_ext"]*(1-e_ext["alpha_ext"])),
        minimum_ext = e_ext["min_ext"] + e_ext["amp_ext"] * expit(-e_ext["beta_ext"]*(1+e_ext["alpha_ext"])),
        amplitude_ext = (maximum_ext - minimum_ext),
        mesor_ext = (maximum_ext + minimum_ext)/2,
        hours_above_mesor_ext = 24*mean(e_24 > mesor_ext),
        up_mesor_ext = ( e_ext["phi_ext"] - hours_above_mesor_ext/2)%%24,
        down_mesor_ext = ( e_ext["phi_ext"] + hours_above_mesor_ext/2)%%24,
        rss_cos = sum((cosOut$fit$residuals)^2),
        rss_ext = sum(residuals(exCosOut)^2),
        F_pseudo_ext = ((rss_cos - rss_ext)/2)/
            (rss_ext/(nrow(tmp.dat) - 5))
    )
    if(export_cosinor_param){
        estimates_out = estimates_out %>% mutate(
            mesor_cos = e_cos[1],
            amp_cos = e_cos[2],
            acrotime_cos = e_cos[3],
            minimum_cos = e_cos[1] - e_cos[2]/2
        )
    }
    if(export_nls_outcome){
        estimates_out = estimates_out %>% tibble(
            niter = exCosOut$niter,
            info = exCosOut$info,
            message = exCosOut$message
        )
    }
    ret = list(
        estimates = estimates_out
        #cosinor_ts = cosinor_ts
    )
    if(export_ts){
        ret$cosinor_ts = tibble(
            time = tmp.dat$time,
            time_across_days = tmp.dat$time + (tmp.dat$day-1)*24,
            y = tmp.dat$Y,
        ) %>% mutate(
            y_cos = cosOut$fit$fitted.values,
            y_ext = cosinorExtendedFitted(
                e_ext,
                t = time
            )
        )
    }
    return(ret)
}
