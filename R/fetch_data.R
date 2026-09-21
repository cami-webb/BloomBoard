# Pulls chlorophyll and covariate data for one buoy and joins it into one daily table

source("R/utils.R")

fetch_buoy_data <- function(cfg, buoy_id) {
  buoy <- cfg$buoys[[buoy_id]]

  message("Loading chlorophyll for ", buoy_id)
  # date_range (backtest) mode uses the calibrated chlorophyll file; daily
  # mode uses the live feed. Both sources' value_sd/value_range/value_trend
  # columns are renamed the same way below so downstream code doesn't care
  # which source produced them.
  if (cfg$run_mode == "date_range" && !is.null(buoy$chlorophyll$calibrated_url)) {
    chl <- load_calibrated_chlorophyll(cfg, buoy$chlorophyll$calibrated_url)
  } else {
    chl <- load_sensor_daily(buoy$chlorophyll, include_hourly_stats = TRUE)
  }
  chl <- chl %>%
    dplyr::rename(chlorophyll = value, chlora_hourly_sd = value_sd,
                  chlora_hourly_range = value_range, chlora_hourly_trend = value_trend)

  covariate_names <- names(buoy$covariates)
  data <- chl
  for (cov_name in covariate_names) {
    message("Loading covariate ", cov_name, " for ", buoy_id)
    cov_data <- load_sensor_daily(buoy$covariates[[cov_name]]) %>%
      dplyr::rename(!!cov_name := value)
    data <- dplyr::full_join(data, cov_data, by = "date")
  }

  data %>% dplyr::arrange(date)
}

# Saves the full daily table for one buoy to S3 as the training dataset
save_training_data <- function(cfg, buoy_id, data) {
  fname <- paste0(buoy_id, "_training_data.csv")
  write_s3_csv(cfg, cfg$s3$write_path, fname, data)
}

# Most recent real chlorophyll observation at or before ref_date (NA if none)
find_seed_date <- function(data, ref_date) {
  obs_dates <- data$date[!is.na(data$chlorophyll) & data$date <= ref_date]
  if (length(obs_dates) == 0) return(as.Date(NA))
  max(obs_dates)
}
