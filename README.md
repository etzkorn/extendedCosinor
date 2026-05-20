extendedCosinorVignette
================

``` r
library(extendedCosinor)
library(tidyverse)
#> ── Attaching core tidyverse packages ───────────────────────────────────────────────────────────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.1.4     ✔ readr     2.1.5
#> ✔ forcats   1.0.0     ✔ stringr   1.5.1
#> ✔ ggplot2   3.5.2     ✔ tibble    3.2.1
#> ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
#> ✔ purrr     1.0.4     
#> ── Conflicts ─────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
library(cluster)
#> Warning: package 'cluster' was built under R version 4.5.2
```

## Data Requirements

We’ll simulate seven days of minute-level data `x` for one person for
this demonstration. `cosinorExtendedModel` and `cosinorBasicModel` both
take a vector of epoch-level data as input, as well as the argument
`window` which specifies the epoch duration in minutes. The vector must
have a length that is an integer multiple of `1440/60`. The model
assumes the first element of x corresponds to the first epoch after
midnight.

``` r
t_hours <- seq(0, 24*7,by = 1/60)[-1]
x <- sqrt(exp(2 + 3*cos((t_hours- 5)/24*2*pi) + rnorm(length(t_hours),1)))
plot(x~t_hours)
```

![](../README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
cosOut <- cosinorExtendedModel(
    x = x, 
    window = 1, 
    export_ts = TRUE,
    export_cosinor_param = TRUE,
    export_nls_outcome = TRUE
)
cosOut$estimates
#>     mu_ext gamma_ext alpha_ext beta_ext  phi_ext ndays acrotime_ext maximum_ext minimum_ext amplitude_ext mesor_ext
#> 1 1.001684  43.00079 0.9999951 2.346702 5.026552     7     5.026552    22.50221     1.39181       21.1104  11.94701
#>   hours_above_mesor_ext up_mesor_ext down_mesor_ext up_mesor_deriv_ext  rss_ext    tss_y F_pseudo_ext    R2_ext mesor_cos  phi_cos
#> 1              7.633333     1.209885       8.843218           4.212311 371515.8 952019.9     862.7549 0.6097605  8.464432 5.012861
#>    amp_cos rss_cos    R2_cos niter info                                                 message
#> 1 10.12693  435144 0.5429255    44    1 Relative error in the sum of squares is at most `ftol'.
```

``` r
cosOut$cosinor_ts
#> # A tibble: 10,080 × 4
#>    time_across_days     y y_cos y_ext
#>               <dbl> <dbl> <dbl> <dbl>
#>  1           0.0167  2.04  11.1  7.39
#>  2           0.0333  3.89  11.1  7.45
#>  3           0.05    2.22  11.2  7.50
#>  4           0.0667  8.45  11.2  7.56
#>  5           0.0833  6.21  11.3  7.61
#>  6           0.1    11.8   11.3  7.67
#>  7           0.117   2.21  11.4  7.72
#>  8           0.133  14.3   11.4  7.78
#>  9           0.15    3.96  11.4  7.83
#> 10           0.167   4.90  11.5  7.89
#> # ℹ 10,070 more rows
```

## Interpreting the Output

### Primary Model Parameter Estimates

The five primary model parameter estimates are the following:

``` r
cosOut$estimates %>% select(mu_ext:phi_ext)
#>     mu_ext gamma_ext alpha_ext beta_ext  phi_ext
#> 1 1.001684  43.00079 0.9999951 2.346702 5.026552
```

These correspond to the model outlined by Marler et al (2006). We’ve
renamed some of the parameters ($\mu_ext$ was “min” and $\gamma_ext$ was
“amp”) from the original version to avoid confusion with the amplitude
and minimum of the fitted function.

$$
    y(t) = \mu_{\mathrm{ext}} + \gamma_{\mathrm{ext}}\cdot \mathrm{expit}\bigg(\beta_{\mathrm{ext}} \Big(\cos(\left[t - \phi_{\mathrm{ext}}\right]\cdot 2\pi/24)-\alpha_{\mathrm{ext}} \Big) \bigg) + \epsilon(t)
$$

The plot below shows that $mu_ext$ (formerly “min”) can be close to the
function minumum, but it is not exact. Additionally, $gam_ext$ (formerly
“amp”) influences the function amplitude, but it is not the amplitude of
the fitted function.

``` r
cosOut$cosinor_ts %>% filter(
    time_across_days <= 24
) %>% ggplot() + geom_point(
    aes(
        x = time_across_days, y = y
    ), alpha = 0.2, color = "blue"
)  + geom_line(
    aes(
        x = time_across_days, y = y_ext
    )
) + scale_x_continuous(
    "Time of Day (hours since midnight)",
    breaks = 4*(0:6)
) + geom_hline(
    "hline",
    yintercept = c(
        cosOut$estimates$mu_ext
    ), color = "red"
) + geom_segment(
    aes(
        xend = cosOut$estimates$acrotime_ext+12,
        x = cosOut$estimates$acrotime_ext+12,
        y = cosOut$estimates$mu_ext,
        yend = cosOut$estimates$mu_ext + 
            cosOut$estimates$gamma_ext
    ), color = "purple", arrow = arrow(ends = "both")
) + geom_label(
    data = tibble(
        name = c("mu_ext", "mu_ext + gamma_ext"),
        x = c(cosOut$estimates$acrotime_ext, cosOut$estimates$acrotime_ext + 12),
        y = c(
            cosOut$estimates$mu_ext, 
            cosOut$estimates$mu_ext + cosOut$estimates$gamma_ext + 3
        )
    ),
    aes(
        x = x,
        y = y,
        label = name
    )
) + ggtitle(
    "Primary Estimates from extendedCosinor"
)
```

![](../README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

### Secondary Estimates

Secondary estimates are more appropriate for analysis, and are derived
from the fitted function and residuals.

``` r
cosOut$estimates %>% select(acrotime_ext:up_mesor_deriv_ext,F_pseudo_ext,R2_ext)
#>   acrotime_ext maximum_ext minimum_ext amplitude_ext mesor_ext hours_above_mesor_ext up_mesor_ext down_mesor_ext up_mesor_deriv_ext
#> 1     5.026552    22.50221     1.39181       21.1104  11.94701              7.633333     1.209885       8.843218           4.212311
#>   F_pseudo_ext    R2_ext
#> 1     862.7549 0.6097605
```

The connection of each estimate to the fitted function is depicted
below.

``` r
label.df = cosOut$estimates %>% select(
    acrotime_ext:up_mesor_deriv_ext
) %>% pivot_longer(
    acrotime_ext:up_mesor_deriv_ext
) %>% mutate(
    x = ifelse(
        name %in% c(
            "minimum_ext",
            "mesor_ext",
            "maximum_ext"
        ), 20, value
    ), 
    y = ifelse(
        !name %in% c(
            "minimum_ext",
            "mesor_ext",
            "maximum_ext"
        ), 60 + c(0,-15,3, 6), value
    ),
    y = ifelse(
        name == "amplitude_ext",
        cosOut$estimates$maximum_ext*0.75,
        y
    ),
    x = ifelse(
        name == "amplitude_ext",
        cosOut$estimates$acrotime_ext+11.5,
        x
    ),
    y = ifelse(
        name == "hours_above_mesor_ext",
        38,
        y
    ),
    x = ifelse(
        name == "hours_above_mesor_ext",
        cosOut$estimates$acrotime_ext,
        x
    ),
    y = ifelse(
        name == "up_mesor_deriv_ext",
        cosOut$estimates$mesor_ext + (15*cosOut$estimates$up_mesor_deriv_ext),
        y
    ),
    x = ifelse(
        name == "up_mesor_deriv_ext",
        cosOut$estimates$up_mesor_ext + 15,
        x
    )
)

cosOut$cosinor_ts %>% filter(
    time_across_days <= 24
) %>% ggplot() + geom_point(
    aes(
        x = time_across_days, y = y
    ), alpha = 0.2, color = "blue"
)  + geom_line(
    aes(
        x = time_across_days, y = y_ext
    )
) + scale_x_continuous(
    "Time of Day (hours since midnight)",
    breaks = 4*(0:6)
) + geom_hline(
    "hline",
    yintercept = c(
        cosOut$estimates$minimum_ext,
        cosOut$estimates$mesor_ext,
        cosOut$estimates$maximum_ext
    ), color = "red"
) + geom_vline(
    xintercept = c(
        cosOut$estimates$up_mesor_ext,
        cosOut$estimates$down_mesor_ext,
        cosOut$estimates$acrotime_ext
    ), color = "orange"
) + geom_segment(
    aes(
        xend = cosOut$estimates$down_mesor_ext,
        x = cosOut$estimates$up_mesor_ext,
        y = 40, yend = 40
    ), color = "green", arrow = arrow(ends = "both")
) + geom_segment(
    aes(
        xend = cosOut$estimates$acrotime_ext+12,
        x = cosOut$estimates$acrotime_ext+12,
        y = cosOut$estimates$minimum_ext,
        yend = cosOut$estimates$minimum_ext +
            cosOut$estimates$amplitude_ext
    ), color = "purple", arrow = arrow(ends = "both")
) + geom_abline(
    aes(
        intercept = cosOut$estimates$mesor_ext - 
            cosOut$estimates$up_mesor_deriv_ext*
            cosOut$estimates$up_mesor_ext,
        slope = cosOut$estimates$up_mesor_deriv_ext
    ), color = "cyan"
) + geom_label(
    data = label.df,
    aes(
        x = x,
        y = y,
        label = name
    )
) + ggtitle(
    "Secondary Estimates from extendedCosinor"
)
```

![](../README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

### Comparison to ActCR::ActExtendCosinor

``` r
actcrOut <- ActCR::ActExtendCosinor(
    x = x, 
    window = 1, 
    export_ts = TRUE
)
actcrOut$params %>% as.data.frame
#>    minimum      amp alpha     beta acrotime F_pseudo  UpMesor DownMesor   MESOR ndays
#> 1 1.001678 43.00105     1 2.346694 5.026552 862.7552 5.026552  5.026552 22.5022     7
```

The plot below shows how certain estimates from the
`ActCR::ActExtendCosinor` output map to the fitted function. The plot
below confirms that the estimates labelled `amp` and `minimum` do not
correspond to the amplitude and minimum of the fitted function.
Additionally, the estimate of `MESOR`, `UpMesor`, and `DownMesor` do not
seem to capture the intended construct–the midline of rhythm.

``` r
actcrOut$cosinor_ts %>% filter(
    time_across_days <= 24
) %>% ggplot() + geom_point(
    aes(
        x = time_across_days, y = original
    ), alpha = 0.2, color = "blue"
)  + geom_line(
    aes(
        x = time_across_days, y = fittedYext
    )
) + scale_x_continuous(
    "Time of Day (hours since midnight)",
    breaks = 4*(0:6)
) + geom_hline(
    "hline",
    yintercept = c(
        actcrOut$params$minimum,
        actcrOut$params$MESOR,
        actcrOut$params$minimum + actcrOut$params$amp
    ), color = "red"
) + geom_vline(
    xintercept = c(
        actcrOut$params$acrotime,
        actcrOut$params$UpMesor,
        actcrOut$params$DownMesor
    ), color = "orange"
) + geom_segment(
    aes(
        xend = actcrOut$params$acrotime + 12,
        x = actcrOut$params$acrotime + 12,
        y = actcrOut$params$minimum,
        yend = actcrOut$params$minimum + actcrOut$params$amp
    ), color = "purple", arrow = arrow(ends = "both")
)  + geom_label(
    data = tibble(
        name = c("minimum", "amp", "acrotime", "MESOR", "UpMesor", "DownMesor"),
        x = c(
            20, actcrOut$params$acrotime + 12, actcrOut$params$acrotime, 
            20, actcrOut$params$UpMesor, actcrOut$params$DownMesor
        ),
        y = c(
            actcrOut$params$minimum, 
            actcrOut$params$minimum + actcrOut$params$amp + 3,
            actcrOut$params$MESOR*2 + 15,
            actcrOut$params$MESOR,
            actcrOut$params$MESOR*2+10,
            actcrOut$params$MESOR*2 + 5
        )
    ),
    aes(
        x = x,
        y = y,
        label = name
    )
) + ggtitle(
    "Estimates from ActCR::ActExtendCosinor"
)
```

![](../README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

### Cosinor Model Estimates

By setting the option `export_cosinor_param = TRUE`, the function
includes basic cosinor parameter estimates in the output. These
estimates are mapped onto the basic cosinor fit below.

``` r
cosOut$cosinor_ts %>% filter(
    time_across_days <= 24
) %>% ggplot() + geom_point(
    aes(
        x = time_across_days, y = y
    ), alpha = 0.2, color = "blue"
)  + geom_line(
    aes(
        x = time_across_days, y = y_ext
    ), linetype = 2
) + geom_line(
    aes(
        x = time_across_days, y = y_cos
    ), linewidth = 1
) + scale_x_continuous(
    "Time of Day (hours since midnight)",
    breaks = 4*(0:6)
) + geom_hline(
    "hline",
    yintercept = cosOut$estimates$mesor_cos, color = "red"
) + geom_vline(
    xintercept = cosOut$estimates$phi_cos, color = "orange"
) + geom_segment(
    aes(
        yend = cosOut$estimates$mesor_cos - cosOut$estimates$amp_cos,
        y = cosOut$estimates$mesor_cos + cosOut$estimates$amp_cos,
        x = cosOut$estimates$phi_cos + 12, 
        yend = cosOut$estimates$phi_cos + 12
    ), color = "green", arrow = arrow(ends = "both")
) + geom_label(
    data = tibble(
        name = c("mesor_cos", "2*amp_cos", "phi_cos"),
        x = c(
            20,
            cosOut$estimates$phi_cos + 12, 
            cosOut$estimates$phi_cos
        ),
        y = c(
            cosOut$estimates$mesor_cos, 
            cosOut$estimates$mesor_cos + cosOut$estimates$amp_cos + 3,
            cosOut$estimates$mesor_cos + cosOut$estimates$amp_cos + 10
        )
    ),
    aes(
        x = x,
        y = y,
        label = name
    )
) + ggtitle(
    "Basic Cosinor Estimates from extendedCosinorModel"
)
```

![](../README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

### Comparing Basic and Extended Cosinor Model Estimates

The figure above shows that the extended cosinor model (dashed line)
appears to fit the observed data better than the basic cosinor model
(solid line). The estimates of the fitted amplitude, acrotime, and mesor
can be compared between both models.

``` r
cosOut$estimates %>% transmute(
    amplitude_ext, 2*amp_cos, 
    acrotime_ext, phi_cos,
    mesor_ext, mesor_cos
)
#>   amplitude_ext 2 * amp_cos acrotime_ext  phi_cos mesor_ext mesor_cos
#> 1       21.1104    20.25386     5.026552 5.012861  11.94701  8.464432
```

The pseudo F-statistic indicates the degree to which the extended
cosinor model fits the data better than the basic cosinor model.
`R2_ext` tells us that 61% of the variation in the epoch level data can
be predicted by the extended cosinor fit. `R2_cos` tells us that 54.3%
of the variation in the epoch level data can be predicted by the basic
cosinor fit.

``` r
cosOut$estimates %>% select(
    F_pseudo_ext, R2_ext, R2_cos
)
#>   F_pseudo_ext    R2_ext    R2_cos
#> 1     862.7549 0.6097605 0.5429255
```
