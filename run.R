# BloomBoard main workflow. Run this for a daily update or a backtest via run_mode in config.yml

# For each buoy and each date to run:
#   1. fetch new chlorophyll and covariate data
#   2. save the training data to S3
#   3. forecast each covariate forward
#   4. forecast chlorophyll forward using those covariate forecasts
#   5. save the forecast output to S3

source("R/utils.R")
source("R/fetch_data.R")
source("R/forecast_covariates.R")
source("R/forecast_chlorophyll.R")

library(foreach)
library(doParallel)

n_threads <- as.integer(Sys.getenv("NSLOTS", unset = "1"))
message("Using ", n_threads, " threads")
doParallel::registerDoParallel(cores = n_threads)

cfg <- load_config("config.yml")
run_dates <- resolve_run_dates(cfg)

for (buoy_id in names(cfg$buoys)) {
  message("Buoy: ", buoy_id)

  data <- fetch_buoy_data(cfg, buoy_id)
  save_training_data(cfg, buoy_id, data)

  foreach(ref_date = as.list(run_dates), .packages = c("dplyr", "ranger")) %dopar% {
    ref_date <- as.Date(ref_date)

    seed_date <- find_seed_date(data, ref_date)
    if (is.na(seed_date)) {
      message("No chlorophyll observation at or before ", ref_date, ", skipping")
      return(NULL)
    }
    message("Forecasting from reference date ", ref_date, " (seeded from ", seed_date, ")")

    cov_fc <- forecast_covariates(cfg, buoy_id, data, seed_date, cfg$max_horizon)
    chl_fc <- forecast_chlorophyll(cfg, buoy_id, data, cov_fc, seed_date, cfg$max_horizon)

    if (is.null(chl_fc)) {
      message("Skipping ", ref_date, ", no forecast produced")
      return(NULL)
    }

    fname <- paste0(buoy_id, "-", seed_date, "-forecast.csv")
    write_s3_csv(cfg, cfg$s3$write_path, fname, chl_fc)
    NULL
  }
}

message("Done")
