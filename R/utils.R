# Shared helpers for BloomBoard

library(yaml)
library(httr)
library(readr)
library(dplyr)
library(lubridate)
library(aws.s3)

load_config <- function(path = "config.yml") {
  yaml::read_yaml(path)
}

# Looks up data[[col]] on target_date, falling back to the most recent
# value within max_lookback days, then that day-of-year's climatological
# mean
lookup_with_fallback <- function(data, col, target_date, max_lookback = 5) {
  for (offset in 0:max_lookback) {
    v <- data[[col]][data$date == (target_date - offset)]
    if (length(v) == 1 && !is.na(v)) return(v)
  }
  doy <- lubridate::yday(target_date)
  clim_vals <- data[[col]][!is.na(data[[col]]) & lubridate::yday(data$date) == doy]
  if (length(clim_vals) > 0) return(mean(clim_vals))
  NA_real_
}

# Historical files use a pre-generated static CSV 
# realtime files use NetCDF-to-CSV conversion

# Download of historical CSV 
fetch_direct_csv <- function(path, timeout_s = 300) {
  resp <- httr::GET(
    paste0("http://gyre.umeoce.maine.edu/data/gomoos/buoy/", path),
    httr::timeout(timeout_s)
  )
  if (httr::http_error(resp)) {
    stop("Download failed for ", path, ": ", httr::http_status(resp)$message)
  }
  tmp <- tempfile(fileext = ".csv")
  writeBin(httr::content(resp, as = "raw"), tmp)
  d <- readr::read_csv(tmp, show_col_types = FALSE)
  unlink(tmp)
  names(d)[1] <- "time"
  d
}

# NetCDF to CSV conversion
fetch_converted_csv <- function(ncfile, timeout_s = 170) {
  resp <- httr::GET(
    "http://gyre.umeoce.maine.edu/data/gomoos/buoy/php/view_csv_file.php",
    query = list(ncfile = ncfile),
    httr::timeout(timeout_s)
  )
  if (httr::http_error(resp)) {
    stop("Download failed for ", ncfile, ": ", httr::http_status(resp)$message)
  }
  tmp <- tempfile(fileext = ".csv")
  writeBin(httr::content(resp, as = "raw"), tmp)
  d <- readr::read_csv(tmp, show_col_types = FALSE)
  unlink(tmp)
  names(d)[1] <- "time"
  d
}

# Pulls one variable from a sensor spec, combines historical+realtime, and collapses to daily (mean of sub-daily readings)
# With include_hourly_stats = TRUE, also returns value_sd/value_range/value_trend computed from the underlying hourly readings
load_sensor_daily <- function(spec, include_hourly_stats = FALSE) {
  hist <- fetch_direct_csv(spec$historical_url)
  d <- hist
  if (isTRUE(spec$use_realtime)) {
    rt <- fetch_converted_csv(spec$realtime_url)
    d <- dplyr::bind_rows(hist, rt)
  }
  d <- d %>%
    dplyr::mutate(date = as.Date(time), val = .data[[spec$variable]]) %>%
    dplyr::filter(!is.na(val))

  if (!include_hourly_stats) {
    return(d %>%
      dplyr::group_by(date) %>%
      dplyr::summarise(value = mean(val), .groups = "drop") %>%
      dplyr::filter(!is.nan(value)))
  }

  d <- d %>% dplyr::mutate(hour = lubridate::hour(time))
  d %>%
    dplyr::group_by(date) %>%
    dplyr::summarise(
      value = mean(val),
      value_sd = if (dplyr::n() >= 2) stats::sd(val) else 0,
      value_range = max(val) - min(val),
      value_trend = if (dplyr::n() >= 2) stats::coef(stats::lm(val ~ hour, data = data.frame(val, hour)))[2] else 0,
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.nan(value))
}

# Loads a daily chlorophyll series already stored on S3 for use in date_range/backtest mode
load_calibrated_chlorophyll <- function(cfg, fname) {
  d <- read_s3_csv(cfg, cfg$s3$read_path, fname)
  if (is.null(d)) stop("Calibrated chlorophyll file not found on S3: ", fname)
  d$date <- as.Date(d$date)
  d
}

s3_key <- function(cfg, path, fname) {
  paste0(path, "/", fname)
}

read_s3_csv <- function(cfg, path, fname) {
  key <- s3_key(cfg, path, fname)
  exists <- tryCatch(
    aws.s3::object_exists(object = key, bucket = cfg$s3$bucket, base_url = cfg$s3$base_url,
                          use_https = TRUE, region = "",
                          key = Sys.getenv("OSN_KEY"), secret = Sys.getenv("OSN_SECRET")),
    error = function(e) FALSE
  )
  if (!isTRUE(exists)) return(NULL)
  tmp <- tempfile(fileext = ".csv")
  aws.s3::save_object(object = key, bucket = cfg$s3$bucket, file = tmp, base_url = cfg$s3$base_url,
                      use_https = TRUE, region = "",
                      key = Sys.getenv("OSN_KEY"), secret = Sys.getenv("OSN_SECRET"))
  d <- readr::read_csv(tmp, show_col_types = FALSE)
  unlink(tmp)
  d
}

write_s3_csv <- function(cfg, path, fname, data) {
  tmp <- tempfile(fileext = ".csv")
  readr::write_csv(data, tmp)
  key <- s3_key(cfg, path, fname)
  aws.s3::put_object(file = tmp, object = key, bucket = cfg$s3$bucket, base_url = cfg$s3$base_url,
                     use_https = TRUE, region = "",
                     key = Sys.getenv("OSN_KEY"), secret = Sys.getenv("OSN_SECRET"))
  unlink(tmp)
  message("Uploaded s3://", cfg$s3$bucket, "/", key)
}

# Dates to run for based on run_mode in the config
resolve_run_dates <- function(cfg) {
  if (cfg$run_mode == "daily") {
    return(Sys.Date())
  }
  seq(as.Date(cfg$date_range$start), as.Date(cfg$date_range$end), by = "day")
}
