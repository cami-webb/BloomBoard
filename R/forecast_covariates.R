# Forecasts each buoy covariate forward in time
# One RF per covariate: today's value and day of year predict tomorrow's value
# Applied recursively so day 2 uses the model's own day 1 prediction, etc 
library(ranger)
library(lubridate)
library(dplyr)

# Trains one RF for one covariate (value + day of year -> next day's value)
train_covariate_model <- function(data, cov_name) {
  d <- data.frame(
    value = data[[cov_name]],
    doy   = lubridate::yday(data$date)
  )
  d$next_value <- dplyr::lead(d$value)
  d <- d[stats::complete.cases(d), ]
  if (nrow(d) < 10) return(NULL)
  ranger::ranger(next_value ~ value + doy, data = d, num.trees = 200,
                 min.node.size = 5, num.threads = 1, seed = 42)
}

# Recursively forecasts one covariate forward from seed value
forecast_one_covariate <- function(model, ref_date, seed_value, horizon) {
  values <- numeric(horizon)
  current_value <- seed_value
  for (h in seq_len(horizon)) {
    target_date <- ref_date + h
    pred <- predict(model, data = data.frame(value = current_value, doy = lubridate::yday(target_date)))$predictions
    values[h] <- pred
    current_value <- pred
  }
  data.frame(date = ref_date + seq_len(horizon), value = values)
}

# Trains and forecasts every covariate in a buoy's covariate_combo forward from ref_date
# Falls back to all configured covariates if none chosen yet
forecast_covariates <- function(cfg, buoy_id, data, ref_date, horizon) {
  buoy <- cfg$buoys[[buoy_id]]
  covs <- buoy$covariate_combo
  if (length(covs) == 0) covs <- names(buoy$covariates)

  result <- data.frame(date = ref_date + seq_len(horizon))
  for (cov_name in covs) {
    train_data <- data[data$date <= ref_date, ]
    model <- train_covariate_model(train_data, cov_name)
    if (is.null(model)) next

    seed_value <- lookup_with_fallback(data, cov_name, ref_date)
    if (is.na(seed_value)) next

    fc <- forecast_one_covariate(model, ref_date, seed_value, horizon)
    result[[cov_name]] <- fc$value
  }
  result
}
