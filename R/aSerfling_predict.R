#' Apply Adjusted Serfling Model to Subsequent Time Periods
#'
#' @description
#' Projects an existing Serfling model onto new temporally contiguous data to detect
#' epidemic signals. Requires test data to immediately follow training data chronologically
#' to maintain periodicity.
#'
#' @param sf Model object from \code{\link{aSerfling}} (must contain \code{best_fit},
#'           \code{output}, and \code{cycles} components)
#' @param df_test New data frame with identical structure to training data, containing
#'                subsequent time points. Must include the response variable column used
#'                in original modeling.
#' @details
#' This function extends the surveillance capability of an established \code{aSerfling} model by:
#' \itemize{
#'   \item Automatically generating time indices continuing from the training set
#'   \item Preserving all terms from the original model fit
#'   \item Calculating prediction intervals using the trained coefficients
#'   \item Flagging values exceeding the 95\% upper prediction bound as warnings
#' }
#'
#' Critical requirements:
#' \enumerate{
#'   \item Test data must maintain the same time resolution (weekly/monthly) as training data
#'   \item The first test observation must be the immediate next time point after the last training observation
#'   \item Column names and cycle parameters must match the original model specification
#' }
#'
#' @returns A data frame containing warning results. The value of the warning column is 1 for warning and 0 for no warning.
#' @export
#' @references Wang X, Wu S, MacIntyre CR, et al. Using an adjusted Serfling regression model to improve the early warning at the arrival of peak timing of influenza in Beijing. PLoS One, 2015,10(3):e0119923.
#'
#' @examples
#' data(sample_ili)
#'
#' ## Split into sequential training/test sets
#' df_train <- sample_ili[1:150,]
#' df_test <- sample_ili[151:200,]
#'
#' ## modeling
#' sf <- aSerfling(df_train, 'case', cycles = c(52, 26))
#'
#' ## apply the model to test set
#' pre <- aSerfling_predict(sf, df_test)
#'
#' ## visualize alerts
#' plot(pre$date, pre$case, type = "l")
#' points(pre$date[pre$warning == 1],
#'        pre$case[pre$warning == 1], col = "red")
#'
#' @importFrom stats as.formula lm predict rnorm
aSerfling_predict <- function(sf,df_test){
  output_train <- sf$output

  df_test$T <- seq(nrow(output_train)+1,nrow(output_train)+nrow(df_test),1)
  df_test$T2 <- df_test$T^2

  cycles <- sf$cycles
  for (cycle in cycles) {
    df_test[[paste0("sin_", cycle)]] <- sin(2 * pi * df_test$T / cycle)
    df_test[[paste0("cos_", cycle)]] <- cos(2 * pi * df_test$T / cycle)
  }

  predict_test <- predict(sf$best_fit,df_test,interval="prediction",level=0.95)
  df_test <- cbind(df_test,predict_test)
  df_test$diff <- df_test$case-df_test$upr
  df_test$warning <- ifelse(df_test$diff>0,1,0)

  return(df_test)
}
