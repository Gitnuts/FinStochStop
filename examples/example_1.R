library(FinStochStopLoss)

# Load the data
log_returns <- logReturns("close_price.csv")

# Fit the GARCH model to the time series
preds <- modelPrediction(log_returns, model = "gjrGARCH", garchOrder = c(1, 1), armaOrder = c(0, 0), distribution.model = "sstd")

# Get average for each estimator of probability distribution grouped by sigma
sigma.df <- getSigma(preds)

# Retreive an examlpe row from sigma.df (e.g. 7th row)
row <- sigma.df[7,]

# Set a target distribution
target_distribution <- fGarch::rsstd(100000, mean = row$ave_mean, sd = row$ave_sigma, xi = row$ave_skew, nu = row$ave_shape)

# Optimize the parameters
params.df <- seqmodelOptim(row)

# Get paths for the optimized parameters
params <- params.df[,c("mean", "sigma", "shape", "skew")]
num_steps <- params.df[,"num_steps"]
paths <- ksdiffScore(params, target_distribution, num_steps = num_steps, get_paths = TRUE)

# Visually compare end value of paths with the target distribution
end_values <- paths[,num_steps]
hist(end_values, breaks = 100, col = "lightblue", border = "pink", probability = TRUE)
lines(stats::density(end_values), col = "red", lwd = 2)
lines(stats::density(target_distribution), lwd = 2)



