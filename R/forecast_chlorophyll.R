# Forecasts chlorophyll forward in time every day, differently depending
# on whether today is a "bloom" or "calm" day

library(ranger)
library(dplyr)
library(lubridate)

# Predict the 0.7 quantile (not the mean) of peak/duration/rise/decline on bloom days
bloom_peak_quantile <- 0.7
bloom_duration_quantile <- 0.7
bloom_rise_quantile <- 0.7
bloom_decline_quantile <- 0.7

# Bloom threshold in µg/L
compute_bloom_threshold <- function(data) 5

# Median rise fraction across past bloom events (fallback for rise_fraction_model)
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

# Median decline rate across past bloom events
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
  bloom_rows$doy <- lubridate::yday(bloom_rows$date)   # seasonality

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

  # this event's own fraction of total duration spent rising to peak
  bloom_rows$rise_fraction <- ave(bloom_rows$chlorophyll, bloom_event_id, FUN = function(vals) {
    peak_idx <- which.max(vals)
    n <- length(vals)
    rep(peak_idx / n, n)
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

# Trains the four bloom outcome models (peak magnitude, remaining duration, decline rate, rise fraction) on bloom-only data)
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

  # filtered separately to avoid dropping peak/duration rows
  d_decline <- base_d
  d_decline$decline_rate <- panel$decline_rate
  d_decline <- d_decline[stats::complete.cases(d_decline), ]

  d_rise <- base_d
  d_rise$rise_fraction <- panel$rise_fraction
  d_rise <- d_rise[stats::complete.cases(d_rise), ]

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
  rise_fraction_model <- if (nrow(d_rise) >= 10) {
    ranger::ranger(as.formula(paste("rise_fraction ~", form_rhs)), data = d_rise,
                   num.trees = 200, min.node.size = 5, num.threads = 1, quantreg = TRUE,
                   importance = importance, seed = 42)
  } else NULL

  # oob_bias: mean signed OOB residual, for rRMSE/rBias diagnostics
  attr(peak_model, "target_mean") <- mean(d$remaining_peak)
  attr(peak_model, "oob_bias") <- mean(peak_model$predictions - d$remaining_peak, na.rm = TRUE)
  attr(duration_model, "target_mean") <- mean(d$remaining_duration)
  attr(duration_model, "oob_bias") <- mean(duration_model$predictions - d$remaining_duration, na.rm = TRUE)
  if (!is.null(decline_rate_model)) {
    attr(decline_rate_model, "target_mean") <- mean(d_decline$decline_rate)
    attr(decline_rate_model, "oob_bias") <- mean(decline_rate_model$predictions - d_decline$decline_rate, na.rm = TRUE)
  }
  if (!is.null(rise_fraction_model)) {
    attr(rise_fraction_model, "target_mean") <- mean(d_rise$rise_fraction)
    attr(rise_fraction_model, "oob_bias") <- mean(rise_fraction_model$predictions - d_rise$rise_fraction, na.rm = TRUE)
  }

  list(peak_model = peak_model, duration_model = duration_model,
      decline_rate_model = decline_rate_model, rise_fraction_model = rise_fraction_model,
      predictor_cols = predictor_cols)
}

# Day-index within the bloom event and chlorophyll value on onset day (NA if not in a bloom)
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

# Builds chlorophyll trajectories (one row per member) from predicted peak, rise length, and decline rate
bloom_trajectories <- function(current_state, pred_peak, rise_days, decline_rate, horizon) {
  n_members <- length(pred_peak)
  values <- matrix(NA_real_, n_members, horizon)
  for (h in seq_len(horizon)) {
    rising <- h <= rise_days
    values[, h] <- ifelse(rising,
                          current_state + (h / pmax(rise_days, 1)) * (pred_peak - current_state),
                          pred_peak + (h - rise_days) * decline_rate)
  }
  pmax(values, 0)
}

# Bloom-day ensemble forecast: each member draws its own peak magnitude, remaining duration,
# rise fraction, and decline rate, then builds its trajectory from those numbers
# u_bloom: n_members x 4 matrix of random quantile levels (peak, duration, rise, decline)
# shift_q: quantile the draws are centered on (NULL for the median)
# noise_sd: daily noise standard deviation, one value or one per lead time
forecast_bloom_ens <- function(train_data, ref_row, covs, ref_date, horizon, bloom_threshold,
                               n_members, u_bloom = NULL, shift_q = NULL, noise_sd = 0) {
  models <- train_bloom_models(train_data, covs, bloom_threshold)
  if (is.null(models)) {
    message("Not enough bloom-period training data, skipping")
    return(NULL)
  }
  if (is.null(u_bloom)) {
    u_bloom <- matrix(runif(n_members * 4), n_members, 4,
                      dimnames = list(NULL, c("peak", "duration", "rise", "decline")))
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
  # hourly chlorophyll stats, default to 0 if unavailable
  for (col in c("chlora_hourly_sd", "chlora_hourly_range", "chlora_hourly_trend")) {
    v <- ref_row[[col]][1]
    newdata[[col]] <- if (!is.null(v) && !is.na(v)) v else 0
  }
  for (cov_name in covs) {
    now_val <- lookup_with_fallback(train_data, cov_name, ref_date)
    newdata[[cov_name]] <- now_val
    past_val <- lookup_with_fallback(train_data, cov_name, ref_date - 7)
    # fall back to 0 if the sensor has a gap
    newdata[[paste0(cov_name, "_trend7")]] <- if (!is.na(now_val) && !is.na(past_val)) now_val - past_val else 0
  }
  
  if (anyNA(newdata)) {
    message("Missing predictor at ", ref_date, ", skipping")
    return(NULL)
  }
  
  # shift each model's draws so their center moves from the median to the shift_q quantile (from config)
  qshift <- function(model) {
    if (is.null(shift_q) || is.null(model)) return(0)
    q <- predict(model, data = newdata, type = "quantiles", quantiles = c(0.5, shift_q))$predictions
    q[1, 2] - q[1, 1]
  }
  
  # one random draw per member from each model's predicted distribution
  pred_peak <- pmax(draw_quantile(models$peak_model, newdata, u_bloom[, "peak"]) + qshift(models$peak_model),
                    current_state)
  pred_remaining_duration <- pmax(1, round(draw_quantile(models$duration_model, newdata, u_bloom[, "duration"]) +
                                             qshift(models$duration_model)))
  
  rise_fraction <- if (!is.null(models$rise_fraction_model)) {
    draw_quantile(models$rise_fraction_model, newdata, u_bloom[, "rise"]) + qshift(models$rise_fraction_model)
  } else {
    rep(compute_bloom_rise_fraction(train_data, bloom_threshold), n_members)
  }
  rise_fraction[is.na(rise_fraction)] <- compute_bloom_rise_fraction(train_data, bloom_threshold)
  rise_fraction <- pmin(pmax(rise_fraction, 0.05), 0.95)
  
  # skip rise phase if there's little predicted rise left
  rise_gap <- pred_peak - current_state
  rise_days <- ifelse(rise_gap <= 0.2, 0, pmax(1, round(pred_remaining_duration * rise_fraction)))
  
  decline_rate <- if (!is.null(models$decline_rate_model)) {
    draw_quantile(models$decline_rate_model, newdata, u_bloom[, "decline"]) + qshift(models$decline_rate_model)
  } else {
    rep(compute_bloom_decline_rate(train_data, bloom_threshold), n_members)
  }
  decline_rate[is.na(decline_rate)] <- 0
  
  trajectories <- bloom_trajectories(current_state, pred_peak, rise_days, decline_rate, horizon)
  
  # daily noise (multiplicative, mean 1) with one standard deviation per lead time
  sd_matrix <- matrix(rep_len(noise_sd, horizon), n_members, horizon, byrow = TRUE)
  trajectories * exp(matrix(rnorm(n_members * horizon), n_members, horizon) * sd_matrix - sd_matrix^2 / 2)
}

# Trains a simple RF on calm-only day-to-day transitions (today's value and covariates predict tomorrow's change)
train_calm_model <- function(data, covs, bloom_threshold, importance = "none", quantreg = TRUE) {
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
                        quantreg = quantreg, importance = importance, seed = 42)
  attr(fit, "target_mean") <- mean(panel$chlora_delta_target)
  fit
}

# Calm-day ensemble forecast: recursive, each member uses its own forecasted covariates
# cov_ens: array [member, day, covariate] from forecast_covariates_ens()
# u_chl: n_members x horizon matrix of random quantile levels
forecast_calm_ens <- function(train_data, ref_row, cov_ens, covs, ref_date, horizon, bloom_threshold,
                              n_members, u_chl = NULL) {
  model <- train_calm_model(train_data, covs, bloom_threshold)
  if (is.null(model)) {
    message("Not enough calm-period training data, skipping")
    return(NULL)
  }
  if (is.null(u_chl)) u_chl <- matrix(runif(n_members * horizon), n_members, horizon)
  
  members <- matrix(NA_real_, n_members, horizon)
  current <- rep(ref_row$chlorophyll[1], n_members)
  
  for (h in seq_len(horizon)) {
    if (h == 1) {
      # real covariates on ref_date w/ gap-tolerant fallback
      real_cov <- as.data.frame(lapply(covs, function(cn) lookup_with_fallback(train_data, cn, ref_date)))
      names(real_cov) <- covs
      cov_now <- real_cov[rep(1, n_members), , drop = FALSE]
    } else {
      cov_now <- as.data.frame(matrix(cov_ens[, h, ], nrow = n_members, dimnames = list(NULL, covs)))
    }
    if (anyNA(cov_now)) break
    
    newdata <- cov_now
    newdata$chlora_state <- current
    delta <- draw_quantile(model, newdata, u_chl[, h])
    current <- pmax(0, current + delta)
    members[, h] <- current
  }
  members
}

# data: the buoy's full daily training table (date, chlorophyll, covariates)
# Returns an n_members x horizon matrix of chlorophyll forecasts (NULL if no forecast possible)
forecast_chlorophyll_ens <- function(cfg, buoy_id, data, ref_date, horizon, n_members,
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
    noise <- cfg$ensemble$bloom_noise_sd
    noise_sd <- noise$day1 + (noise$last_day - noise$day1) * (seq_len(horizon) - 1) / max(horizon - 1, 1)
    forecast_bloom_ens(train_data, ref_row, covs, ref_date, horizon, bloom_threshold, n_members,
                       shift_q = cfg$ensemble$bloom_shift_quantile, noise_sd = noise_sd)
    } else {
    cov_ens <- forecast_covariates_ens(cfg, buoy_id, data, ref_date, horizon, n_members)
    forecast_calm_ens(train_data, ref_row, cov_ens, covs, ref_date, horizon, bloom_threshold, n_members)
  }
}