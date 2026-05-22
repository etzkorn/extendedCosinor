## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>",
    fig.height = 6,
    fig.width = 10
)

## ----setup--------------------------------------------------------------------
library(extendedCosinor)
library(tidyverse)
library(cluster)

## -----------------------------------------------------------------------------
t_hours <- seq(0, 24*7,by = 1/60)[-1]
x <- sqrt(exp(2 + 3*cos((t_hours- 5)/24*2*pi) + rnorm(length(t_hours),1)))
plot(x~t_hours)

## -----------------------------------------------------------------------------
cosOut <- cosinorExtendedModel(
    x = x, 
    window = 1, 
    export_ts = TRUE,
    export_cosinor_param = TRUE,
    export_nls_outcome = TRUE
)
cosOut$estimates

## -----------------------------------------------------------------------------
cosOut$cosinor_ts

## -----------------------------------------------------------------------------
cosOut$estimates %>% select(mu_ext:phi_ext)

## ----warning = FALSE, message = FALSE-----------------------------------------
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

## -----------------------------------------------------------------------------
cosOut$estimates %>% select(acrotime_ext:up_mesor_deriv_ext,F_pseudo_ext,R2_ext)

## ----warning = FALSE, message = FALSE-----------------------------------------
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

## -----------------------------------------------------------------------------
actcrOut <- ActCR::ActExtendCosinor(
    x = x, 
    window = 1, 
    export_ts = TRUE
)
actcrOut$params %>% as.data.frame

## ----warning = FALSE, message = FALSE-----------------------------------------
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

## ----warning = FALSE, message = FALSE-----------------------------------------
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
        yend = cosOut$estimates$mesor_cos - 0.5*cosOut$estimates$amplitude_cos,
        y = cosOut$estimates$mesor_cos + 0.5*cosOut$estimates$amplitude_cos,
        x = cosOut$estimates$phi_cos + 12, 
        yend = cosOut$estimates$phi_cos + 12
    ), color = "green", arrow = arrow(ends = "both")
) + geom_label(
    data = tibble(
        name = c("mesor_cos", "amplitude_cos", "phi_cos"),
        x = c(
            20,
            cosOut$estimates$phi_cos + 12, 
            cosOut$estimates$phi_cos
        ),
        y = c(
            cosOut$estimates$mesor_cos, 
            cosOut$estimates$mesor_cos + 0.5*cosOut$estimates$amplitude_cos + 3,
            cosOut$estimates$mesor_cos + 0.5*cosOut$estimates$amplitude_cos + 10
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

## -----------------------------------------------------------------------------
cosOut$estimates %>% transmute(
    amplitude_ext, amplitude_cos, 
    acrotime_ext, phi_cos,
    mesor_ext, mesor_cos
)

## -----------------------------------------------------------------------------
cosOut$estimates %>% select(
    F_pseudo_ext, R2_ext, R2_cos
)

