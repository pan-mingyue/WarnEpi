#' Exponentially Weighted Moving Average
#' @description
#' Detects anomalies in infectious disease surveillance data using an
#' Exponentially Weighted Moving Average (EWMA) algorithm.
#' Designed for time series data, it flags potential outbreaks by smoothing
#' past observations with decayed weights and comparing against control thresholds.
#'
#' @param data A data frame containing the warning indicator columns, arranged in time-based order.
#' @param column A column name or column number, used to specify the warning indicator.
#' @param lambda The weight factor \eqn{\lambda}, ranging from 0 to 1(higher values prioritize recent observations).
#' @param k The standard deviation coefficient \eqn{k}.
#' @param move_t The moving period \eqn{t_{move}}.
#' @param ignore_t The number of nearest time units to be ignored by the model, \eqn{t_{ignore}}.
#'
#' @returns A data frame containing warning results. The value of the warning column is 1 for warning and 0 for no warning.
#' @details
#' Let \eqn{\mathbf{X} = (X_1,\ldots,X_T)^\top} be an observed time series of disease case counts,
#' where \eqn{X_t} represents the aggregated counts at time \eqn{t} (e.g., daily, weekly, or monthly observations).
#' We assume \eqn{X_t \sim N(\mu, \sigma^2)} for the underlying distribution.
#'
#' The EWMA (Exponentially Weighted Moving Average) model is defined as:
#' \deqn{Z_1 = X_1}
#' \deqn{Z_t = \lambda X_t + (1-\lambda)Z_{t-1}}
#' \deqn{UCL_t = \hat{\mu}_t + k\hat{\sigma}_t\sqrt{\frac{\lambda}{2-\lambda}}}
#'
#' where:
#' \itemize{
#'   \item \eqn{Z_t}: The EWMA statistic at time \eqn{t}, representing an exponentially weighted average of current and past observations.
#'   \item \eqn{\lambda}: Weight factor (\eqn{0 < \lambda < 1}), higher values prioritize recent observations
#'   \item \eqn{k}: Standard deviation coefficient (typically 2-3)
#'   \item \eqn{UCL_t}: Upper Control Limit at time \eqn{t}, forming a dynamic threshold for anomaly detection.
#'   \item \eqn{\hat{\mu}_t, \hat{\sigma}_t}: Estimated from moving window \eqn{(X_{t-t_{move}-t_{ignore}},\ldots,X_{t-1-t_{ignore}})}
#' }
#'
#' An alarm is triggered when \eqn{Z_t > UCL_t}, with the alarm set defined as:
#' \deqn{\mathcal{T} = \{t: Z_t > UCL_t\}}
#'
#' @export
#' @references Wang X, Zeng D, Seale H, et al. Comparing early outbreak detection algorithms based on their optimized parameter values. J Biomed Inform, 2010,43(1):97-103.
#'
#' @examples
#' ## simulate reported cases
#' set.seed(123)
#' cases <- c(round(rnorm(10, 10, 1)), seq(12,21,3), seq(15,5,-5))
#' dates <- seq(as.Date("2025-01-01"), by = "7 days", length.out = length(cases))
#' data_frame <- data.frame(date = dates, case = cases)
#'
#' ## modeling
#' output <- EWMA(data_frame,'case',lambda = 0.5, k = 3, move_t = 4, ignore_t = 2)
#' output
#'
#' ## visualize alerts
#' plot(output$date, output$case, type = "l")
#' points(output$date[output$warning == 1],
#'        output$case[output$warning == 1], col = "red")
#'
#' @importFrom stats sd rnorm
EWMA <- function(data,column,lambda=0.5,k=3,move_t,ignore_t=2){
  dt <- as.data.frame(data)
  dt$z <- dt[,column]
  dt <- dt[which(is.na(dt$z)==FALSE),]
  for(i in 2:nrow(dt)){
    dt[i,'z'] <- lambda*dt[i,column]+(1-lambda)*dt[i-1,'z']
  }
  for(i in (1+move_t+ignore_t):nrow(dt)){
    dt[i,'mu'] <- mean(dt[(i-move_t-ignore_t):(i-ignore_t-1),column])
    dt[i,'sigma'] <- sd(dt[(i-move_t-ignore_t):(i-ignore_t-1),column])
    dt[i,'UCL'] <- dt[i,'mu']+k*dt[i,'sigma']*sqrt(lambda/(2-lambda))
    dt[i,'subtract'] <- dt[i,'z']-dt[i,'UCL']
    dt[i,'warning'] <- ifelse(dt[i,'subtract']>0,1,0)
  }
  return(dt)
}
