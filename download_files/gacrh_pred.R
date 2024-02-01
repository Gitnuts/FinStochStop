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

ks_test <- function(params, target_distribution, MJD = FALSE, get_distr = FALSE, get_paths = FALSE, num_steps, nsim = 1) {

  num_paths <- 10000
  num_steps <- num_steps
  ks_statistic_list <- list()

  for (i in 1:nsim) {
  paths <- replicate(num_paths, {

    mu <- params[1]
    sigma <- params[2]
    nu <- params[3]
    theta <- params[4]
    lambda <- params[5]


    if (length(params) > 2) {
      log_returns <- rvg(n = num_steps, vgC = mu, sigma = sigma, nu = nu, theta = theta)
    } else {
      log_returns <- rnorm(num_steps, mean = mu, sd = sigma)
    }

    if (MJD & length(params) > 4) {
      jumps <- rpois(num_steps, lambda = lambda)
      jumps[as.logical(jumps)] <- rnorm(length(which(jumps != 0)), mean = mu, sd = sigma)
      cumsum(jumps) + cumsum(log_returns)

    } else {
      cumsum(log_returns)
    }

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
  ks_statistic_list[i] <- ks_statistic
  }

  return(mean(unlist(ks_statistic_list)))
}}

estimator_seq <- function(sigma.df, estimator, num_steps = 10) {

  sigma_upper_bound <- 5 / num_steps
  if (estimator == "mean") {
    estimator_x <- seq(from = sigma.df$ave_mean, to = 0, length.out = 100)
  } else if (estimator == "sigma") {
    estimator_x <- seq(from = 0.1*sigma.df$ave_sigma, to = 0.2 * sigma.df$ave_sigma , length.out = 100)
  } else if (estimator == "shape") {
    estimator_x <- seq(from = sigma.df$ave_shape -2, to = sigma.df$ave_shape +2, length.out = 100)
  } else if (estimator == "skew") {
    estimator_x <- seq(from = (sigma.df$ave_skew-1)/1000, to = abs(sigma.df$ave_skew-1)/1000, length.out = 100)
  } else if (estimator == "lambda") {
    estimator_x <- seq(from = 0, to = 0.1, length.out = 100)
  }

  return(estimator_x)
}

estimator_local_minimum <- function(estimator_seq, estimator_list, polynomial = 2, get_score = FALSE) {

  estimator.df <- as.data.frame(cbind(estimator_seq, estimator_list))
  colnames(estimator.df) <- c("estimator", "ks_error")
  estimator.df$estimator <- as.numeric(estimator.df$estimator)
  estimator.df$ks_error <- as.numeric(estimator.df$ks_error)
  estimator.fit <- lm(ks_error ~ poly(estimator, polynomial), data = estimator.df)
  estimator_fitted_line <- predict(estimator.fit, data = estimator.df)
  min_estimator <- estimator_seq[which.min(estimator_fitted_line)]

  if (get_score) {
    return(list(min_estimator, min(estimator_fitted_line)))
  }
  return(min_estimator)
}

single_fit <- function(sigma.df, num_steps = 30, MJD = FALSE) {


split.df <- data.frame(sigma_bar = double(), number_of_splits = integer(),
                       mean = double(), sigma = double(),
                       shape = double(), skew = double(),
                       lambda = double())


for (row in seq(1, nrow(sigma.df), by = 1)) {

  target_distribution <- rsstd(10000, mean = sigma.df[row,]$ave_mean,
                               sd = sigma.df[row,]$ave_sigma,
                               nu = sigma.df[row,]$ave_shape,
                               xi = sigma.df[row,]$ave_skew)


  mean_x <- estimator_seq(sigma.df, "mean")
  mean_list <- sapply(mean_x, function(mean_val) ks_test(c(mean_val, 0.1), target_distribution, num_steps = num_steps))
  min_mean <- mean_x[which.min(mean_list)]

  sigma_x <- estimator_seq(sigma.df, "sigma", num_steps = num_steps)
  sigma_list <- sapply(sigma_x, function(sigma_val) ks_test(c(min_mean, sigma_val), target_distribution, num_steps = num_steps))
  min_sigma <- estimator_local_minimum(sigma_x, sigma_list, polynomial = 3)

  skew_x <- estimator_seq(sigma.df, "skew")
  skew_list <- sapply(skew_x, function(skew_val) ks_test(c(min_mean, min_sigma, sigma.df[row,]$ave_shape, skew_val), target_distribution, num_steps = num_steps))
  min_skew <- estimator_local_minimum(skew_x, skew_list, polynomial = 3)

  nu_x <- estimator_seq(sigma.df, "shape")
  nu_list <- sapply(nu_x, function(nu_val) ks_test(c(min_mean, min_sigma, nu_val, min_skew), target_distribution, num_steps = num_steps))
  min_nu <- estimator_local_minimum(nu_x, nu_list, polynomial = 3)

  if (MJD) {
    lambda_x <- estimator_seq(sigma.df, "lambda")
    lambda_list <- sapply(lambda_x, function(lambda_val) ks_test(c(min_mean, min_sigma, min_nu, min_skew, lambda_val), target_distribution, num_steps = num_steps, MJD = TRUE))
    min_lambda <- estimator_local_minimum(lambda_x, lambda_list, polynomial = 3)

  } else {
    min_lambda <- NA_real_
  }

  split.df[nrow(split.df) + 1,] = c(sigma.df[row,]$sigma_bar, num_steps,
                                    min_mean, min_sigma, min_nu, min_skew,
                                    min_lambda)
  #print(paste("row processed for sigma bar =", sigma.df[row,]$sigma_bar))
}

  return(split.df)
}

####################################################################################
####################################################################################

stop_loss <- function(params, paths, trigger = FALSE, summary = FALSE) {

  take_profit_condition <- params[1]
  stop_loss_condition <- params[2]
  upper_trigger <- params[3]

  # Identify the time steps where conditions are met
  num_paths <- 1000
  sampled_indices <- sample(seq(ncol(paths)), num_paths)


  tp_indices <- which(paths[,sampled_indices] >= take_profit_condition, arr.ind = TRUE)
  sl_indices <- which(paths[,sampled_indices] <= stop_loss_condition, arr.ind = TRUE)
  tp_first_incident <- tapply(tp_indices[,1], tp_indices[,2], min)
  sl_first_incident <- tapply(sl_indices[,1], sl_indices[,2], min)

  fixed_returns <- merge(as.matrix(tp_first_incident),
                         as.matrix(sl_first_incident),
                         by=0, all = TRUE)

  fixed_returns <- merge(as.matrix(seq(1, num_paths, length.out = num_paths)),
                         fixed_returns, by.x = 0, by.y = "Row.names", all = TRUE)
  fixed_returns <- subset(fixed_returns, select = -c(V1))
  colnames(fixed_returns) <- c("row", "take_profit_index", "stop_loss_index")

  # if number of parameters is grearer than 2
  # This snippet finds indeces whether, after triggering a level at timestamp T,
  # a timesries bottoms down bellow zero at timestamp T + k.
  # default is long position stop-loss, so to switch to short sell, inequality signs
  # in ´ut_indices´ and ´ut_sl_indices´ must be turned.

  if (trigger) {
  ut_indices <- which(paths[,sampled_indices] >= upper_trigger, arr.ind = TRUE)
  ut_sl_indices <- which(paths[,sampled_indices] < 0, arr.ind = TRUE)
  ut_first_incident <- tapply(ut_indices[,1], ut_indices[,2], min)

  d <- merge(as.matrix(ut_first_incident), ut_sl_indices, by.x = 0, by.y = "col")
  d <- d %>%
    filter(row > V1) %>%
    group_by(Row.names) %>%
    summarize(min_column1 = min(row)) %>%
    left_join(d %>% filter(row > V1), multiple = "any", by = "Row.names")

  d <- merge(as.matrix(ut_first_incident),
             subset(d, select=-c(row, V1)), by.x = 0, by.y = "Row.names", all = TRUE)

  colnames(d) <- c("row", "upper_trigger_level", "upper_trigger_stop_loss")
  fixed_returns <- merge(fixed_returns, d, by ='row', all = TRUE)


  fixed_returns <- fixed_returns %>%
    mutate(return = case_when(

      (upper_trigger_level < stop_loss_index | (is.na(stop_loss_index) & !is.na(upper_trigger_level))) & (is.na(take_profit_index) & is.na(upper_trigger_stop_loss)) ~ take_profit_condition/2,
      (upper_trigger_level < stop_loss_index | (is.na(stop_loss_index) & !is.na(upper_trigger_level))) & (is.na(take_profit_index) & !is.na(upper_trigger_stop_loss)) ~ 0,
      (upper_trigger_level < stop_loss_index | (is.na(stop_loss_index) & !is.na(upper_trigger_level))) & (!is.na(take_profit_index) & is.na(upper_trigger_stop_loss)) ~ take_profit_condition,
      (upper_trigger_level < stop_loss_index | (is.na(stop_loss_index) & !is.na(upper_trigger_level))) & (take_profit_index < upper_trigger_stop_loss) ~ take_profit_condition,
      (upper_trigger_level < stop_loss_index | (is.na(stop_loss_index) & !is.na(upper_trigger_level))) & (take_profit_index > upper_trigger_stop_loss) ~ 0,
      upper_trigger_level > stop_loss_index | (!is.na(stop_loss_index) & is.na(upper_trigger_level)) ~ stop_loss_condition,
      TRUE ~ NA_real_
    ))
  } else {
    fixed_returns <- fixed_returns %>%
    mutate(return = case_when(
            is.na(take_profit_index) & is.na(stop_loss_index) ~ NA_real_,
            is.na(take_profit_index) | take_profit_index > stop_loss_index ~ stop_loss_condition,
            is.na(stop_loss_index) | take_profit_index < stop_loss_index ~ take_profit_condition,
            TRUE ~ NA_real_
    ))
  }

  counts <- table(fixed_returns$return)
  take_profit_counts <- as.numeric(counts[as.character(take_profit_condition)])
  stop_loss_counts <- as.numeric(counts[as.character(stop_loss_condition)])
  upper_trigger_counts <- as.numeric(counts[as.character(take_profit_condition/2)])
  upper_trigger_stop_loss_counts <- as.numeric(counts[as.character(0)])

  if (is.na(upper_trigger_counts) | is.na(upper_trigger_stop_loss_counts)) {
    upper_trigger_counts <- 0
    upper_trigger_stop_loss_counts <- 0
  }

  if (is.na(take_profit_counts)) {
    take_profit_counts <- 0
  }

  penalty_counts <- num_paths - sum(take_profit_counts, stop_loss_counts, upper_trigger_counts, upper_trigger_stop_loss_counts)

  if (penalty_counts > stop_loss_counts) {
    stop_loss_counts <- penalty_counts
  }

  score <- sum(take_profit_counts * take_profit_condition,
               stop_loss_counts * stop_loss_condition,
               upper_trigger_counts * take_profit_condition/2,
               upper_trigger_stop_loss_counts * 0)

  if (summary) {
    print(paste("stop_loss_counts:", stop_loss_counts))
    print(paste("take_profit_counts:", take_profit_counts))
    print(paste("upper_trigger_counts:", upper_trigger_counts))
    print(paste("upper_trigger_stop_loss_counts:", upper_trigger_stop_loss_counts))
    print(paste("No hits:", penalty_counts))
    print(paste("return per trade:", score / (num_paths - penalty_counts)))
  }

  return(score*(-1))
}

####################################################################################
####################################################################################

jump_diffusion_test <- function(params, num_steps, target_df, sigma_bar, get_paths = FALSE, nsim = 1) {

  # extracting actual highs and lows for sigma_bar from data frame given positive (negative) returns.

  real.lows.positive_returns <- target_df$low_close_ratio[target_df$Sigma > sigma_bar - 0.0005 & target_df$Sigma < sigma_bar & target_df$close_log_return >= 0]
  real.highs.positive_returns <- target_df$high_close_ratio[target_df$Sigma > sigma_bar - 0.0005 & target_df$Sigma < sigma_bar & target_df$close_log_return >= 0]

  real.lows.negative_returns <- target_df$low_close_ratio[target_df$Sigma > sigma_bar - 0.0005 & target_df$Sigma < sigma_bar & target_df$close_log_return < 0]
  real.highs.negative_returns <- target_df$high_close_ratio[target_df$Sigma > sigma_bar - 0.0005 & target_df$Sigma < sigma_bar & target_df$close_log_return < 0]

  num_paths <- 10000
  num_steps <- num_steps

  mu <- params[1]
  sigma <- params[2]
  nu <- params[3]
  theta <- params[4]

  # params vectors contains parameters in order to exctract paths
  lambda <- params[5]
  delta <- params[6]
  tau <- params[7]

  # setting a grid where:
  #   - lambda := jump intensity, specified by Poisson
  #   - delta  := jump amplitide, i.e a multiplier for scaled jumps
  #   - tau    := threshold for returns to be scaled with s.t. global minimum (maximum) of a path
  #               with positive (negative) end value (i.e. value at time t = T) is scaled by a tau iff
  #               value of global minimum (maximum) is positive (negative).

  grid_values <- expand.grid(lambda = seq(0, 1, length.out = 10),
                             delta = seq(0, 1, length.out = 10),
                             tau = 0.0005
                             #tau = seq(0, sigma/5, length.out = 10)
  )

  subtract_min_conditionally <- function(path, tau) {
    if (min(path) > 0 && tail(path, 1) > 0) {
      min_index <- which.min(path)
      path[min_index] <- path[min_index] - tau
    }
    path
  }

  add_max_conditionally <- function(path, tau) {
    if (max(path) < 0 && tail(path, 1) < 0) {
      max_index <- which.max(path)
      path[max_index] <- path[max_index] + tau
    }
    path
  }

  grid_search <- function(grid_value, get_paths = FALSE) {

    lambda <- as.numeric(grid_value["lambda"])
    delta <- as.numeric(grid_value["delta"])
    tau <- as.numeric(grid_value["tau"])

    for (i in 1:nsim) {
      paths <- replicate(num_paths, {

        log_returns <- rvg(n = num_steps, vgC = mu, sigma = sigma, nu = nu, theta = theta)
        jumps <- rpois(num_steps, lambda = lambda)
        jumps[as.logical(jumps)] <- rnorm(length(which(jumps != 0)), mean = mu, sd = sigma)
        cumulative_jumps <- cumsum(jumps)

        cumsum(log_returns) + delta * (cumulative_jumps - mean(cumulative_jumps))

      })
    }

    paths <- apply(paths, 2, subtract_min_conditionally, tau = tau)
    paths <- apply(paths, 2, add_max_conditionally, tau = tau)

    # add first row for time step t = 0
    paths <- rbind(numeric(ncol(paths)), paths)

    if (get_paths) {
      return(paths)
    }

    final_values <- paths[num_steps+1,]

    highs.positive_returns <- apply(paths[, which(final_values >= 0)], MARGIN = 2, FUN = max)
    lows.positive_returns <- apply(paths[, which(final_values >= 0)], MARGIN = 2, FUN = min)
    highs.negative_returns <- apply(paths[, which(final_values < 0)], MARGIN = 2, FUN = max)
    lows.negative_returns <- apply(paths[, which(final_values < 0)], MARGIN = 2, FUN = min)

    # Due to non-continiuos nature of distributions, i.e. distribution contains ties (repeated values),
    # k-s test is not robust solution for this task. Instead, Wasserstein distance is calculated between
    # two samples.

    ks_res <- c(
         wasserstein1d(real.highs.positive_returns, highs.positive_returns),
         wasserstein1d(real.lows.positive_returns, lows.positive_returns),
         wasserstein1d(real.highs.negative_returns, highs.negative_returns),
         wasserstein1d(real.lows.negative_returns, lows.negative_returns)
    )

    return(mean(ks_res))
  }

  if (get_paths && length(params) == 7) {
    paths <- grid_search(data.frame(lambda=lambda, delta=delta, tau=tau), get_paths = TRUE)
    return(paths)
  }

  result <- apply(grid_values, 1, grid_search)
  grid_results <- cbind(grid_values, result)

  # returns a table with wasserstain scores given parameters
  return(grid_results)

}