# Forecasts each buoy covariate forward in time as an ensemble
# One quantile-regression RF per covariate (today's value and day of year predict tomorrow's value)
# Applied recursively, each ensemble member drawing its own random quantile at every step
library(ranger)
library(lubridate)
library(dplyr)

# Trains one RF for one covariate (value + day of year -> next day's value)
train_covariate_model <- function(data, cov_name, quantreg = TRUE) {
  d <- data.frame(
    value = data[[cov_name]],
    doy   = lubridate::yday(data$date)
  )
  d$next_value <- dplyr::lead(d$value)
  d <- d[stats::complete.cases(d), ]
  if (nrow(d) < 10) return(NULL)
  ranger::ranger(next_value ~ value + doy, data = d, num.trees = 200,
                 min.node.size = 5, num.threads = 1, quantreg = quantreg, seed = 42)
}

# Recursively forecasts one covariate for every ensemble member
# u: n_members x horizon matrix of random quantile levels
forecast_one_covariate_ens <- function(model, ref_date, seed_value, horizon, u) {
  n_members <- nrow(u)
  values <- matrix(NA_real_, n_members, horizon)
  current <- rep(seed_value, n_members)
  for (h in seq_len(horizon)) {
    newdata <- data.frame(value = current, doy = lubridate::yday(ref_date + h))
    current <- draw_quantile(model, newdata, u[, h])
    values[, h] <- current
  }
  values
}

# Trains and forecasts every covariate in a buoy's covariate_combo forward from ref_date
# Returns an array [member, day, covariate]
forecast_covariates_ens <- function(cfg, buoy_id, data, ref_date, horizon, n_members, u_cov = NULL) {
  buoy <- cfg$buoys[[buoy_id]]
  covs <- buoy$covariate_combo
  if (length(covs) == 0) covs <- names(buoy$covariates)
  
  if (is.null(u_cov)) {
    u_cov <- array(runif(n_members * horizon * length(covs)), c(n_members, horizon, length(covs)))
  }
  result <- array(NA_real_, c(n_members, horizon, length(covs)), dimnames = list(NULL, NULL, covs))
  
  train_data <- data[data$date <= ref_date, ]
  for (k in seq_along(covs)) {
    model <- train_covariate_model(train_data, covs[k])
    if (is.null(model)) next
    
    seed_value <- lookup_with_fallback(data, covs[k], ref_date)
    if (is.na(seed_value)) next
    
    result[, , k] <- forecast_one_covariate_ens(model, ref_date, seed_value, horizon, u_cov[, , k])
  }
  result
}
