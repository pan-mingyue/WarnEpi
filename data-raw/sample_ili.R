#' Sample Influenza-like Illness (ILI) Data
#'
#' A simulated dataset containing weekly ILI case counts and dates.
#' @format A data frame with 200 rows and 2 variables:
#' \describe{
#'   \item{date}{Date in YYYY-MM-DD format}
#'   \item{case}{Integer count of reported cases}
#' }
"sample_ili"

#
set.seed(1)
T <- c(1:200)
cycle1 <- 52
cycle2 <- 26
sin1 <- sin(2*pi*T/cycle1)
cos1 <- cos(2*pi*T/cycle1)
sin2 <- sin(2*pi*T/cycle2)
cos2 <- cos(2*pi*T/cycle2)
error <- stats::rnorm(length(T), 0, 1)
y <- 200 + 0.1*T - 0.0002*T^2 + 26*sin1 + 21*cos1 + 5*sin1 + 28*cos1 + error
epi <- c(stats::rnorm(51, 0, 10),
         seq(50, 170, 30), seq(150, 30, -30),
         stats::rnorm(49, 0, 15),
         seq(100, 300, 50), seq(350, 70, -70),
         stats::rnorm(35, 0, 10),
         seq(30, 110, 20), seq(140, 60, -20),
         stats::rnorm(35, 0, 10))
y_new <- round(y + epi, 0)
cases <- round(y_new, 0)
dates <- seq(as.Date("2021-01-01"), by = "7 days", length.out = length(cases))
sample_ili <- data.frame(date = dates, case = cases)

# save
save(sample_ili, file = "data/sample_ili.rda", compress = "xz")

#
rm(list = setdiff(ls(), "sample_ili"))
