#' cosinorExtendedModel0
#'
#' @description
#' Fit the five parameter cosinor model from Marler et al. (2006).
#' Using the current approach of ActExtendCosinor in ActCR.
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
cosinorExtendedModel0 <- function(
        x,
        window = 1,
        export_ts = FALSE,
        export_cosinor_param = TRUE,
        export_nls_outcome = FALSE,
        maxiter = 1000,
        maxfev = 600
){
    dim = 1440/window
    n.days = length(x)/dim
    tmp.dat = data.frame(time = rep(1:dim, n.days)/(60/window),
                         Y = x)
    fit = cosinor.lm(Y ~ time(time) + 1, data = tmp.dat, period = 24)
    mesor = fit$coefficients[1]
    amp = fit$coefficients[2]
    acr = correct.acrophase(fit)
    acrotime = (-1) * acr * 24/(2 * pi)
    names(mesor) = names(amp) = names(acr) = names(acrotime) = NULL
    e_min0 = max(mesor - amp, 0)
    e_amp0 = 2 * amp
    e_phi0 = acrotime
    e_par0 = c(e_min0, e_amp0, 0, 2, e_phi0)
    tmp.dat = tmp.dat[!is.na(tmp.dat$Y), ]
    fit_nls = nls.lm(e_par0, fn = fn_obj, lower = lower, upper = upper,
                     tmp.dat = tmp.dat, control = nls.lm.control(maxiter = 1000))
    if (export_ts == TRUE) {
        fittedYext = tmp.dat$Y - residuals(fit_nls)
        fittedY = fitted(fit$fit)
        original = tmp.dat$Y
        time = tmp.dat$time
        time2 = time
        drops = which(diff(time) < 0) + 1
        for (k in drops) {
            time2[k:length(time)] = time2[k:length(time)] + 24
        }
        time_across_days = time2
        cosinor_ts = as.data.frame(cbind(time, time_across_days,
                                         original, fittedY, fittedYext))
    }
    else {
        cosinor_ts = NULL
    }
    coef.nls = coef(fit_nls)
    e_min = coef.nls[1]
    e_amp = coef.nls[2]
    e_alpha = coef.nls[3]
    e_beta = coef.nls[4]
    e_acrotime = coef.nls[5]
    RSS_cos = sum((fit$fit$residuals)^2)
    RSS_ext = sum(residuals(fit_nls)^2)
    F_pseudo = ((RSS_cos - RSS_ext)/2)/(RSS_ext/(nrow(tmp.dat) -
                                                     5))
    UpMesor = -acos(e_alpha)/(2 * pi/24) + e_acrotime
    DownMesor = acos(e_alpha)/(2 * pi/24) + e_acrotime
    MESOR = e_min + e_amp/2
    params = tibble(
        minimum0 = e_min,
        amp0 = e_amp,
        alpha0 = e_alpha,
        beta0 = e_beta,
        acrotime0 = e_acrotime,
        F_pseudo0 = F_pseudo,
        UpMesor0 = UpMesor,
        DownMesor0 = DownMesor,
        MESOR0 = MESOR,
        ndays0 = n.days
    )

    if (export_nls_outcome == TRUE) {
        params = bind_cols(
            params,
            tibble(
                niter0 = fit_nls$niter,
                info0 = fit_nls$info,
                message0 = fit_nls$message
            )
        )
    }
    ret = list(params = params, cosinor_ts = cosinor_ts)
    return(ret)
}
