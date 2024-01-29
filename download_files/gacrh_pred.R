library(dplyr)

# Setting up log returns of financial time-series:
log_returns <- function(path_to_csv_file) {

  library(xts)

  close_prices <- read.csv(path_to_csv_file, header=TRUE, sep=",")
  df$timestamp <- as.POSIXct(df$timestamp, format = "%Y-%m-%d %H:%M:%S")
  df <- na.omit(df)
  df$log_returns <- append(diff(log(df$close_price)), 0, 0)
  xts_log_returns <- xts(df$log_returns, order.by = df$timestamp)

  return(xts_log_returns)
}

####################################################################################
####################################################################################

# Selecting conditional volatility model to get PDF statistics for each time stamp:

model_selection <- function(model = "sGARCH",
                             ts = log_returns,
                             garchOrder = c(1, 1),
                             armaOrder = c(0, 0),
                             distribution.model = "std") {
  library(rugarch)

  garch_model <- ugarchspec(
      variance.model = list(model = model, garchOrder = garchOrder),
      mean.model = list(armaOrder),
      distribution.model = distribution.model
    )

  garchroll <- ugarchroll(spec = garch_model,
                          data = ts[1:(nrow(ts)/4 + 2016)],
                          n.start = 2016,
                          refit.window = "moving", refit.every = 7*288)

  preds <- as.data.frame(garchroll)
  e  <- preds$Realized - preds$Mu
  d  <- e^2 - preds$Sigma^2

  print(paste("MSE of", model, "model:", mean(d^2)))

  return(preds)
}

####################################################################################
####################################################################################

# Getting statistics from conditional volatility model grouped by mean sigma:

statsGarch <- functional(preds, precision_bar = 0.0005) {
    precision_bar <- 0.0005

    sigma.df <- data.frame(sigma_bar = double(), ave_mean = double(), ave_sigma = double(),
                       ave_shape = double(), ave_skew = double())

    for (sigma_bar in seq(from = precision_bar, to = round(max(preds$Sigma), 3), by = precision_bar)) {

        l <- sigma_bar - precision_bar

        ave_mean <- mean(preds$Mu[preds$Sigma > l & preds$Sigma < sigma_bar])
        ave_sigma <- mean(preds$Sigma[preds$Sigma > l & preds$Sigma < sigma_bar])
        ave_shape <- mean(preds$Shape[preds$Sigma > l & preds$Sigma < sigma_bar])
        ave_skew <- mean(preds$Skew[preds$Sigma > l & preds$Sigma < sigma_bar])

        sigma.df[nrow(sigma.df) + 1,] = c(sigma_bar, ave_mean, ave_sigma, ave_shape, ave_skew)
    }

    return(sigma.df)
}

####################################################################################
####################################################################################

library(MASS)
library(fGarch)
library(VarianceGamma)
library(goftest)

# Example target distribution (replace with your desired distribution)
# target_distribution <- rsstd(10000, mean = mu_mean, sd = sigma_mean, nu = shape_mean, xi = skew_mean)

ks_test <- function(params, target_distribution, num_steps, get_distr = FALSE, get_paths = FALSE) {

  num_paths <- 10000
  num_steps <- num_steps

  paths <- replicate(num_paths, {

    mu <- params[1]
    sigma <- params[2]
    nu <- params[3]
    theta <- params[4]


    if (length(params) > 2) {
      log_returns <- rvg(n = num_steps, vgC = 0, sigma = sigma, nu = nu, theta = theta)
    } else {
      log_returns <- rnorm(num_steps, mean = mu, sd = sigma)
    }

    cumsum(log_returns)
  })

  if (get_paths) {
    return(paths)
  }

  # Extract final values
  final_values <- paths[num_steps,]
  highest_points <- apply(paths, MARGIN = 2, FUN = max)

  if (get_distr) {
    return(final_values)
  }

  ks_statistic <- ks.test(target_distribution, ecdf(final_values))$statistic
  return(ks_statistic)
}