#' Modified Cumulative Sum
#'
#' @description Modified CUSUM method for outbreak detection in infectious disease surveillance data.
#' Implements three variants (C1', C2', C3') with dynamic thresholds for time series analysis.
#'
#'
#' @param data A data frame containing the warning indicator columns, arranged in time-based order.
#' @param column A column name or column number, used to specify the warning indicator.
#' @param k The standard deviation coefficient \eqn{k}.
#' @param h The threshold coefficient \eqn{h}.
#' @param move_t The moving period \eqn{t_{move}}.
#'
#' @details
#' Let \eqn{\mathbf{X} = (X_1,\ldots,X_T)^\top} be an observed time series of disease case counts,
#' where \eqn{X_t} represents the aggregated counts at time \eqn{t} (e.g., daily, weekly, or monthly observations).
#' We assume \eqn{X_t \sim N(\mu, \sigma^2)} for the underlying distribution.
#'
#' The modified CUSUM models accumulate excess cases beyond control limits:
#'
#' \deqn{C1'_0 = C2'_0 = 0}
#' \deqn{C1'_t = \max\left(0, X_t - (\hat{\mu}_t + k\hat{\sigma}_t) + C1'_{t-1}\right)}
#' \deqn{C2'_t = \max\left(0, X_t - (\hat{\mu}_t + k\hat{\sigma}_t) + C2'_{t-1}\right)}
#' \deqn{C3'_t = C2'_t + C2'_{t-1} + C2'_{t-2}}
#' \deqn{H_t = h\hat{\sigma}_t}
#'
#' where:
#' \itemize{
#'   \item \eqn{k}: Standard deviation coefficient (typical range 0.5–1.5), adjusts sensitivity to deviations
#'   \item \eqn{h}: Threshold coefficient (typical range 2–5), controls alarm stringency
#'   \item \eqn{H}: Threshold
#'   }
#'
#' Model specifications:
#' \itemize{
#'   \item \strong{C1'}: Baseline \eqn{\hat{\mu}_t, \hat{\sigma}_t} estimated from \eqn{(X_{t-t_{move}},...,X_{t-1})}
#'   \item \strong{C2'}: Baseline \eqn{\hat{\mu}_t, \hat{\sigma}_t} estimated from \eqn{(X_{t-2-t_{move}},...,X_{t-3})} to avoid recent outbreaks
#'   \item \strong{C3'}: 3-day cumulative sum of C2' values
#'   \item Alarms trigger when \eqn{Cx'_t > H_t} for each model (x = 1,2,3)
#' }
#'
#'
#' @returns A data frame containing C1', C2' and C3' warning results. The value of the warning column is 1 for warning and 0 for no warning.
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
#' output <- mCUSUM(data_frame, 'case', k = 1, h = 2.5, move_t = 4)
#' output
#'
#' ## visualize alerts
#' ### C1'
#' plot(output$date, output$case, type = "l")
#' points(output$date[output$C1_prime_warning == 1],
#'        output$case[output$C1_prime_warning == 1], col = "red")
#'
#' ### C2'
#' plot(output$date, output$case, type = "l")
#' points(output$date[output$C2_prime_warning == 1],
#'        output$case[output$C2_prime_warning == 1], col = "red")
#'
#' ### C3'
#' plot(output$date, output$case, type = "l")
#' points(output$date[output$C3_prime_warning == 1],
#'        output$case[output$C3_prime_warning == 1], col = "red")
#'
#' @importFrom stats sd rnorm
mCUSUM <- function(data,column,k=1,h=2,move_t){
  dt <- as.data.frame(data)
  dt <- dt[which(is.na(dt[,column])==FALSE),]
  for(i in (move_t+1):nrow(dt)){
    dt[i,'mu'] <- mean(dt[(i-move_t):(i-1),column])
    dt[i,'sigma'] <- sd(dt[(i-move_t):(i-1),column])
    dt[,'threshold'] <- h*dt[,'sigma']
  }
  #
  dt[move_t+1,'C1_prime'] <- max(0,dt[move_t+1,column]-dt[move_t+1,'mu']-k*dt[move_t+1,'sigma'])
  for(i in (move_t+2):nrow(dt)){
    dt[i,'C1_prime'] <- max(0,dt[i,column]-dt[i,'mu']-k*dt[i,'sigma']+dt[i-1,'C1_prime'])
  }
  dt[,'C1_prime_subtract'] <- dt[,'C1_prime']-h*dt[,'sigma']
  dt[,'C1_prime_warning'] <- ifelse(dt[,'C1_prime_subtract']>0,1,0)
  #
  dt[move_t+3,'C2_prime'] <- max(0,dt[move_t+3,column]-dt[move_t+1,'mu']-k*dt[move_t+1,'sigma'])
  for(i in (move_t+4):nrow(dt)){
    dt[i,'C2_prime'] <- max(0,dt[i,column]-dt[i-2,'mu']-k*dt[i-2,'sigma']+dt[i-1,'C2_prime'])
  }
  dt[,'C2_prime_subtract'] <- dt[,'C2_prime']-h*dt[,'sigma']
  dt[,'C2_prime_warning'] <- ifelse(dt[,'C2_prime_subtract']>0,1,0)
  #
  for(i in (move_t+5):nrow(dt)){
    dt[i,'C3_prime'] <- sum(dt[(i-2):i,'C2_prime'])
  }
  dt[,'C3_prime_subtract'] <- dt[,'C3_prime']-h*dt[,'sigma']
  dt[,'C3_prime_warning'] <- ifelse(dt[,'C3_prime_subtract']>0,1,0)

  return(dt)
}
