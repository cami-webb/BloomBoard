# Forecasts chlorophyll forward in time every day, differently depending
# on whether today is a "bloom" or "calm" day. See README for full design

library(ranger)
library(dplyr)
library(lubridate)

# Predict this upper quantile (not the mean) of peak/duration/decline on bloom days (there was meaningful underestimation when
# predicting the mean and a 0.7 quantile produced the best results)
bloom_peak_quantile <- 0.7
bloom_duration_quantile <- 0.7
bloom_decline_quantile <- 0.7

# Fixed bloom threshold (µg/L) (see README for why not dynamic)
compute_bloom_threshold <- function(data) 5

# Empirical median fraction of an event's duration from start to peak
compute_bloom_rise_fraction <- function(data, threshold) {
  d <- data[order(data$date), ]
  is_bloom <- !is.na(d$chlorophyll) & d$chlorophyll > threshold
  rle_obj <- rle(is_bloom)
  event_id_all <- rep(seq_along(rle_obj$lengths), rle_obj$lengths)

  bloom_idx <- which(is_bloom)
  bloom_rows <- d[bloom_idx, ]
  bloom_event_id <- event_id_all[bloom_idx]

  chl_by_event <- split(bloom_rows$chlorophyll, bloom_event_id)
  fractions <- vapply(chl_by_event, function(vals) which.max(vals) / length(vals), numeric(1))
  median(fractions, na.rm = TRUE)
}

# Median per-day decline rate (peak to event end) across historical bloom events with an observed decline phase
compute_bloom_decline_rate <- function(data, threshold) {
  d <- data[order(data$date), ]
  is_bloom <- !is.na(d$chlorophyll) & d$chlorophyll > threshold
  rle_obj <- rle(is_bloom)
  event_id_all <- rep(seq_along(rle_obj$lengths), rle_obj$lengths)

  bloom_idx <- which(is_bloom)
  bloom_rows <- d[bloom_idx, ]
  bloom_event_id <- event_id_all[bloom_idx]

  chl_by_event <- split(bloom_rows$chlorophyll, bloom_event_id)
  rates <- vapply(chl_by_event, function(vals) {
    peak_idx <- which.max(vals)
    n <- length(vals)
    if (peak_idx >= n) return(NA_real_)   # peak was the last observed day, no decline phase seen
    (vals[n] - vals[peak_idx]) / (n - peak_idx)
  }, numeric(1))
  median(rates, na.rm = TRUE)
}

# Today's actual chlorophyll minus a naive extrapolation of yesterday's trend
compute_recent_forecast_error <- function(data) {
  d <- data$date
  chl <- data$chlorophyll
  val_at <- function(offset) chl[match(d - offset, d)]
  yesterday <- val_at(1)
  day_before <- val_at(2)
  naive_expected <- yesterday + (yesterday - day_before)
  chl - naive_expected
}

# One row per day inside a detected bloom event: 
# days_since_onset
# rise_rate_since_onset
# remaining_peak
# remaining_duration
# recent_forecast_error
# each covariate + its 7-day trend
# this event's own decline_rate
compute_bloom_training_panel <- function(data, covs, threshold) {
  d <- data[order(data$date), ]
  d$recent_forecast_error <- compute_recent_forecast_error(d)
  is_bloom <- !is.na(d$chlorophyll) & d$chlorophyll > threshold
  rle_obj <- rle(is_bloom)
  event_id_all <- rep(seq_along(rle_obj$lengths), rle_obj$lengths)

  bloom_idx <- which(is_bloom)
  bloom_rows <- d[bloom_idx, ]
  bloom_event_id <- event_id_all[bloom_idx]

  onset_value <- ave(bloom_rows$chlorophyll, bloom_event_id, FUN = function(x) x[1])
  bloom_rows$days_since_onset <- ave(seq_along(bloom_event_id), bloom_event_id, FUN = seq_along)
  bloom_rows$rise_rate_since_onset <- (bloom_rows$chlorophyll - onset_value) / bloom_rows$days_since_onset
  bloom_rows$doy <- lubridate::yday(bloom_rows$date)   # seasonality -- see README

  bloom_rows$remaining_peak <- ave(bloom_rows$chlorophyll, bloom_event_id,
                                   FUN = function(x) rev(cummax(rev(x))))
  event_len <- ave(bloom_rows$days_since_onset, bloom_event_id, FUN = max)
  bloom_rows$remaining_duration <- event_len - bloom_rows$days_since_onset + 1

  # this event's own peak-to-end decline rate (NA if peak was the last day)
  bloom_rows$decline_rate <- ave(bloom_rows$chlorophyll, bloom_event_id, FUN = function(vals) {
    peak_idx <- which.max(vals)
    n <- length(vals)
    if (peak_idx >= n) return(rep(NA_real_, n))
    rep((vals[n] - vals[peak_idx]) / (n - peak_idx), n)
  })

  dts_full <- d$date
  for (cov_name in covs) {
    val_full <- d[[cov_name]]
    now_val  <- val_full[match(bloom_rows$date, dts_full)]
    past_val <- val_full[match(bloom_rows$date - 7, dts_full)]
    bloom_rows[[paste0(cov_name, "_trend7")]] <- now_val - past_val
  }

  bloom_rows
}

# Trains the three bloom outcome models (peak magnitude, remaining duration, decline rate) on bloom-only data)
train_bloom_models <- function(data, covs, bloom_threshold, importance = "none") {
  panel <- compute_bloom_training_panel(data, covs, bloom_threshold)
  trend_cols <- paste0(covs, "_trend7")
  hourly_cols <- c("chlora_hourly_sd", "chlora_hourly_range", "chlora_hourly_trend")
  predictor_cols <- c("chlora_state", "days_since_onset", "rise_rate_since_onset", "doy",
                      "recent_forecast_error", hourly_cols, covs, trend_cols)

  base_d <- data.frame(chlora_state = panel$chlorophyll,
                       days_since_onset = panel$days_since_onset,
                       rise_rate_since_onset = panel$rise_rate_since_onset,
                       doy = panel$doy,
                       recent_forecast_error = panel$recent_forecast_error)
  for (col in hourly_cols) base_d[[col]] <- panel[[col]]
  for (col in c(covs, trend_cols)) base_d[[col]] <- panel[[col]]

  d <- base_d
  d$remaining_peak <- panel$remaining_peak
  d$remaining_duration <- panel$remaining_duration
  d <- d[stats::complete.cases(d), ]
  if (nrow(d) < 10) return(NULL)

  # filtered separately so NA decline_rate rows don't also shrink peak_model/duration_model's training data
  d_decline <- base_d
  d_decline$decline_rate <- panel$decline_rate
  d_decline <- d_decline[stats::complete.cases(d_decline), ]

  form_rhs <- paste0("`", predictor_cols, "`", collapse = " + ")
  # fixed seed (ranger's bootstrap sampling is random otherwise)
  peak_model <- ranger::ranger(as.formula(paste("remaining_peak ~", form_rhs)), data = d,
                               num.trees = 200, min.node.size = 5, num.threads = 1, quantreg = TRUE,
                               importance = importance, seed = 42)
  duration_model <- ranger::ranger(as.formula(paste("remaining_duration ~", form_rhs)), data = d,
                                   num.trees = 200, min.node.size = 5, num.threads = 1, quantreg = TRUE,
                                   importance = importance, seed = 42)
  decline_rate_model <- if (nrow(d_decline) >= 10) {
    ranger::ranger(as.formula(paste("decline_rate ~", form_rhs)), data = d_decline,
                   num.trees = 200, min.node.size = 5, num.threads = 1, quantreg = TRUE,
                   importance = importance, seed = 42)
  } else NULL

  # attached for diagnostic use (e.g. rRMSE) see README
  attr(peak_model, "target_mean") <- mean(d$remaining_peak)
  attr(duration_model, "target_mean") <- mean(d$remaining_duration)
  if (!is.null(decline_rate_model)) attr(decline_rate_model, "target_mean") <- mean(d_decline$decline_rate)

  list(peak_model = peak_model, duration_model = duration_model,
      decline_rate_model = decline_rate_model, predictor_cols = predictor_cols)
}

# ref_date's day-index within its bloom event (1 = onset day) and the chlorophyll value on onset day
# Both NA if ref_date itself isn't currently above threshold
find_bloom_onset_info <- function(data, ref_date, threshold) {
  bloom_dates <- data$date[!is.na(data$chlorophyll) & data$chlorophyll > threshold]
  if (!(ref_date %in% bloom_dates)) return(list(days_since_onset = NA_integer_, onset_value = NA_real_))
  count <- 1L
  check_date <- ref_date - 1
  while (check_date %in% bloom_dates) {
    count <- count + 1L
    check_date <- check_date - 1
  }
  onset_date <- ref_date - (count - 1)
  list(days_since_onset = count, onset_value = data$chlorophyll[data$date == onset_date][1])
}

# Bloom-day forecast: predicts peak magnitude and remaining duration directly, then builds the trajectory from those two numbers
forecast_bloom <- function(train_data, ref_row, covs, ref_date, horizon, bloom_threshold) {
  models <- train_bloom_models(train_data, covs, bloom_threshold)
  if (is.null(models)) {
    message("Not enough bloom-period training data, skipping")
    return(NULL)
  }

  onset_info <- find_bloom_onset_info(train_data, ref_date, bloom_threshold)
  onset_day <- onset_info$days_since_onset
  onset_value <- onset_info$onset_value
  current_state <- ref_row$chlorophyll[1]
  rise_rate <- (current_state - onset_value) / onset_day

  yesterday <- train_data$chlorophyll[train_data$date == ref_date - 1]
  day_before <- train_data$chlorophyll[train_data$date == ref_date - 2]
  recent_forecast_error <- if (length(yesterday) == 1 && length(day_before) == 1 &&
                              !is.na(yesterday) && !is.na(day_before)) {
    current_state - (yesterday + (yesterday - day_before))
  } else 0

  newdata <- data.frame(chlora_state = current_state, days_since_onset = onset_day,
                        rise_rate_since_onset = rise_rate, doy = lubridate::yday(ref_date),
                        recent_forecast_error = recent_forecast_error)
  # within-day chlorophyll stats, fall back to 0 (no known volatility) if this day's hourly readings weren't available
  for (col in c("chlora_hourly_sd", "chlora_hourly_range", "chlora_hourly_trend")) {
    v <- ref_row[[col]][1]
    newdata[[col]] <- if (!is.null(v) && !is.na(v)) v else 0
  }
  for (cov_name in covs) {
    now_val <- lookup_with_fallback(train_data, cov_name, ref_date)
    newdata[[cov_name]] <- now_val
    past_val <- lookup_with_fallback(train_data, cov_name, ref_date - 7)
    # fall back to 0 (no known change) if that sensor has a gap that a few days of lookback still can't fill (a single gappy covariate
    # shouldn't block the whole day's forecast)
    newdata[[paste0(cov_name, "_trend7")]] <- if (!is.na(now_val) && !is.na(past_val)) now_val - past_val else 0
  }

  if (anyNA(newdata)) {
    message("Missing predictor at ", ref_date, ", skipping")
    return(NULL)
  }

  pred_peak <- predict(models$peak_model, data = newdata, type = "quantiles",
                       quantiles = bloom_peak_quantile)$predictions[, 1]
  pred_peak <- max(pred_peak, current_state)   # can't be below what's already been observed
  pred_duration_raw <- predict(models$duration_model, data = newdata, type = "quantiles",
                               quantiles = bloom_duration_quantile)$predictions[, 1]
  pred_remaining_duration <- max(1, round(pred_duration_raw))

  # Build the trajectory from the two predicted numbers:
    # rise to pred_peak over the empirical rise fraction, then decline at decline_rate with no forced endpoint
  rise_fraction <- compute_bloom_rise_fraction(train_data, bloom_threshold)
  # collapse rise phase to 0 if there's little predicted rise left (avoids a flat plateau before the decline)
  rise_gap <- pred_peak - current_state
  rise_days <- if (rise_gap <= 0.2) 0 else max(1, round(pred_remaining_duration * rise_fraction))
  decline_rate <- if (!is.null(models$decline_rate_model)) {
    predict(models$decline_rate_model, data = newdata, type = "quantiles",
           quantiles = bloom_decline_quantile)$predictions[, 1]
  } else {
    compute_bloom_decline_rate(train_data, bloom_threshold)
  }
  if (is.na(decline_rate)) decline_rate <- 0

  values <- numeric(horizon)
  for (h in seq_len(horizon)) {
    if (h <= rise_days) {
      frac <- h / rise_days
      values[h] <- current_state + frac * (pred_peak - current_state)
    } else {
      days_past_peak <- h - rise_days
      values[h] <- pred_peak + days_past_peak * decline_rate
    }
  }

  data.frame(date = ref_date + seq_len(horizon), chlorophyll_forecast = pmax(0, values))
}

# Trains a simple RF on calm-only day-to-day transitions (today's value and covariates predict tomorrow's change)
train_calm_model <- function(data, covs, bloom_threshold, importance = "none") {
  d <- data[order(data$date), ]
  full_delta <- dplyr::lead(d$chlorophyll) - d$chlorophyll
  is_calm <- !is.na(d$chlorophyll) & d$chlorophyll <= bloom_threshold

  panel <- data.frame(chlora_state = d$chlorophyll, chlora_delta_target = full_delta)
  for (cov_name in covs) panel[[cov_name]] <- d[[cov_name]]
  panel <- panel[is_calm, ]
  panel <- panel[stats::complete.cases(panel), ]
  if (nrow(panel) < 10) return(NULL)

  form <- as.formula(paste("chlora_delta_target ~ chlora_state +",
                           paste0("`", covs, "`", collapse = " + ")))
  fit <- ranger::ranger(form, data = panel, num.trees = 200, min.node.size = 5, num.threads = 1,
                        importance = importance, seed = 42)
  attr(fit, "target_mean") <- mean(panel$chlora_delta_target)
  fit
}

# Calm-day forecast: simple recursive using forecasted covariates, where flattening toward a stable
# baseline is realistic
forecast_calm <- function(train_data, ref_row, covariate_forecast, covs, ref_date, horizon, bloom_threshold) {
  model <- train_calm_model(train_data, covs, bloom_threshold)
  if (is.null(model)) {
    message("Not enough calm-period training data, skipping")
    return(NULL)
  }

  predictions <- numeric(horizon)
  current_state <- ref_row$chlorophyll[1]

  for (h in seq_len(horizon)) {
    target_date <- ref_date + h
    if (h == 1) {
      # real covariates on ref_date w/ gap-tolerant fallback
      cov_row <- as.data.frame(lapply(covs, function(cn) lookup_with_fallback(train_data, cn, ref_date)))
      names(cov_row) <- covs
    } else {
      cov_row <- covariate_forecast[covariate_forecast$date == target_date, covs, drop = FALSE]
    }
    if (nrow(cov_row) == 0 || anyNA(cov_row)) {
      predictions[h:horizon] <- NA_real_
      break
    }

    newdata <- cov_row
    newdata$chlora_state <- current_state
    pred_delta <- predict(model, data = newdata)$predictions
    pred <- pmax(0, current_state + pred_delta)
    predictions[h] <- pred
    current_state <- pred
  }

  data.frame(date = ref_date + seq_len(horizon), chlorophyll_forecast = predictions)
}

# data: the buoy's full daily training table (date, chlorophyll, covariates)
# covariate_forecast: output of forecast_covariates() for the same ref_date/horizon
forecast_chlorophyll <- function(cfg, buoy_id, data, covariate_forecast, ref_date, horizon,
                                 bloom_threshold = NULL) {
  buoy <- cfg$buoys[[buoy_id]]
  covs <- buoy$covariate_combo
  if (length(covs) == 0) covs <- names(buoy$covariates)

  train_data <- data[data$date <= ref_date, ]
  if (is.null(bloom_threshold)) bloom_threshold <- compute_bloom_threshold(train_data)

  ref_row <- data[data$date == ref_date, ]
  if (nrow(ref_row) == 0 || is.na(ref_row$chlorophyll[1])) {
    message("No real chlorophyll seed value for ", ref_date, ", skipping")
    return(NULL)
  }

  if (ref_row$chlorophyll[1] > bloom_threshold) {
    forecast_bloom(train_data, ref_row, covs, ref_date, horizon, bloom_threshold)
  } else {
    forecast_calm(train_data, ref_row, covariate_forecast, covs, ref_date, horizon, bloom_threshold)
  }
}
