#' Adjusted Serfling
#'
#' @description Adjusted Serfling regression for periodic disease surveillance, automating epidemic baseline
#' estimation through iterative threshold optimization. Enhances traditional Serfling models by objectively
#' determining epidemic periods and improving peak detection accuracy.
#'
#' @param data A data frame containing the warning indicator columns, arranged in time-based order.
#' @param col_name A column name for the warning indicator (character).
#' @param cycles A numeric vector of disease cycles (e.g., c(52,26) for weekly annual + semi-annual patterns)
#'
#' @details
#' Implements an iterative periodic regression for time series with at least 2 full cycles. Key features:
#'
#' \enumerate{
#'   \item Dynamic Epidemic Filtering:
#'   \itemize{
#'     \item Automatically excludes outbreak points via iterative prediction-CI comparison
#'     \item Terminates when adjusted R-squared stabilizes (maximized model fit)
#'   }
#'
#'   \item Flexible Seasonality Modeling:
#'   \deqn{Y = \beta_0 + \beta_1 t + \beta_2 t^2 + \sum_{k=1}^K \left[\gamma_k \sin\left(\frac{2\pi t}{C_k}\right) + \delta_k \cos\left(\frac{2\pi t}{C_k}\right)\right] + \epsilon}
#'   \itemize{
#'     \item Supports multiple cycles via \code{cycles} parameter (e.g., c(52,26) for weekly annual + semi-annual patterns)
#'     \item Self-adapts to pathogen seasonality shifts
#'   }
#'
#'   \item Peak-Centric Alerting:
#'   \itemize{
#'     \item Flags peaks via optimized threshold (final model's 95\% CI upper bound)
#'     \item Avoids subjective epidemic-onset definitions
#'   }
#' }
#'
#'
#' @returns A list containing:
#' \itemize{
#'   \item output: Full dataset with warning flags (1=alert, 0=normal)
#'   \item best_fit: Final lm model object
#'   \item fit_times: Iteration count for convergence
#'   \item cycles: Input cycle parameters
#' }
#'
#' @export
#' @references Wang X, Wu S, MacIntyre CR, et al. Using an adjusted Serfling regression model to improve the early warning at the arrival of peak timing of influenza in Beijing. PLoS One, 2015,10(3):e0119923.
#'
#' @examples
#' ## modeling
#' data(sample_ili)
#' sf <- aSerfling(data = sample_ili, 'case', cycles = c(52, 26))
#' sf
#'
#' ## visualize alerts
#' output <- sf$output
#' plot(output$date, output$case, type = "l")
#' points(output$date[output$warning == 1],
#'        output$case[output$warning == 1], col = "red")
#'
#' @importFrom stats as.formula lm predict rnorm
aSerfling <- function(data, col_name, cycles) {
  dt <- as.data.frame(data)
  dt$T <- seq(1, nrow(dt), 1)
  dt$T2 <- dt$T^2

  for (cycle in cycles) {
    dt[[paste0("sin_", cycle)]] <- sin(2 * pi * dt$T / cycle)
    dt[[paste0("cos_", cycle)]] <- cos(2 * pi * dt$T / cycle)
  }

  formula_terms <- paste(c("T", "T2", sapply(cycles, function(cycle) c(paste0("sin_", cycle), paste0("cos_", cycle)))), collapse = "+")
  formula <- as.formula(paste(col_name, "~", formula_terms))

  # 1
  fit1 <- lm(formula, data = dt)
  adj_Rsquare_1 <- summary(fit1)$adj.r.squared
  predict1 <- predict(fit1, dt, interval = "prediction", level = 0.95)
  dt1 <- cbind(dt, predict1)
  dt1$diff <- dt1[[col_name]] - dt1$fit
  dt1$case1 <- ifelse(dt1$diff > 0, NA, dt1[[col_name]])

  # 2
  fit2 <- lm(formula, data = dt1, subset = !is.na(dt1$case1))
  adj_Rsquare_2 <- summary(fit2)$adj.r.squared
  predict2 <- predict(fit2, dt, interval = "prediction", level = 0.95)
  dt2 <- cbind(dt, predict2)
  dt2$diff <- dt2[[col_name]] - dt2$upr
  dt2$case1 <- ifelse(dt2$diff > 0, NA, dt2[[col_name]])
  fit_times <- 2
  output <- dt1
  last_fit <- fit1

  # 3
  while (adj_Rsquare_2 > adj_Rsquare_1) {
    fit3 <- lm(formula, data = dt2, subset = !is.na(dt2$case1))
    adj_Rsquare_1 <- adj_Rsquare_2
    adj_Rsquare_2 <- summary(fit3)$adj.r.squared
    last_fit <- fit2
    fit2 <- fit3
    output <- dt2
    fit_times <- fit_times + 1
    predict2 <- predict(fit2, dt, interval = "prediction", level = 0.95)
    dt2 <- cbind(dt, predict2)
    dt2$diff <- dt2[[col_name]] - dt2$upr
    dt2$case1 <- ifelse(dt2$diff > 0, NA, dt2[[col_name]])
  }

  #
  output$warning <- ifelse(output$diff > 0, 1, 0)

  #
  list(output = output, best_fit = last_fit, fit_times = fit_times, cycles = cycles)
}
