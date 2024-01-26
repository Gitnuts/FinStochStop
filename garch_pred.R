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
model_selection <-  function(model = "sGARCH",
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




