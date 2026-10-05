# BloomBoard ensemble workflow. Run this for a daily update or a backtest via run_mode in config.yml
# To run a backtest over a specific date range without editing config.yml:
#   Rscript run_ensemble.R 2021-01-01 2021-12-31

# For each buoy (just A01 right now) and each date:
# 1. fetch new chlorophyll and covariate data
# 2. save the training data to S3 (daily mode only)
# 3. forecast an ensemble of chlorophyll trajectories (covariates are forecast inside the calm path)
# 4. save the ensemble in EFI standard format to S3

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

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) {
  cfg$run_mode <- "date_range"
  cfg$date_range <- list(start = args[1], end = args[2])
}

run_dates <- resolve_run_dates(cfg)
n_members <- cfg$ensemble$n_members
message("Run mode: ", cfg$run_mode, ", ", n_members, " ensemble members")

for (buoy_id in names(cfg$buoys)) {
  message("Buoy: ", buoy_id)
  
  data <- fetch_buoy_data(cfg, buoy_id)
  if (cfg$run_mode == "daily") save_training_data(cfg, buoy_id, data)
  
  # each forecast is seeded from the most recent date with a chlorophyll observation
  seed_dates <- as.Date(vapply(run_dates, function(d) as.numeric(find_seed_date(data, d)), numeric(1)))
  seed_dates <- unique(seed_dates[!is.na(seed_dates)])
  
  # backtests skip dates that already have a forecast on S3
  if (cfg$run_mode == "date_range") {
    existing <- aws.s3::get_bucket_df(
      bucket = cfg$s3$bucket, prefix = paste0(cfg$s3$write_path, "/", buoy_id, "-"),
      base_url = cfg$s3$base_url, use_https = TRUE, region = "",
      key = Sys.getenv("OSN_KEY"), secret = Sys.getenv("OSN_SECRET"), max = Inf
    )$Key
    existing <- existing[grepl(paste0("-", cfg$efi$model_id, "\\.csv$"), existing)]
    done_dates <- as.Date(sub(paste0(buoy_id, "-(.*)-", cfg$efi$model_id, "\\.csv"), "\\1", basename(existing)))
    seed_dates <- seed_dates[!(seed_dates %in% done_dates)]
    message(length(done_dates), " forecasts already on S3, ", length(seed_dates), " left to run")
  }
  
  foreach(seed_date = as.list(seed_dates), .packages = c("dplyr", "ranger")) %dopar% {
    seed_date <- as.Date(seed_date)
    start_time <- Sys.time()
    set.seed(as.integer(seed_date))
    
    members <- forecast_chlorophyll_ens(cfg, buoy_id, data, seed_date, cfg$max_horizon, n_members)
    
    if (is.null(members) || all(is.na(members))) {
      message("Skipping ", seed_date, ", no forecast produced")
      return(NULL)
    }
    
    rows <- efi_ensemble_rows(members, seed_date, cfg, buoy_id)
    rows <- rows[!is.na(rows$prediction), ]
    fname <- paste0(buoy_id, "-", seed_date, "-", cfg$efi$model_id, ".csv")
    write_s3_csv(cfg, cfg$s3$write_path, fname, rows)
    
    message("Forecast for ", seed_date, " took ", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1), " s")
    NULL
  }
}

message("Done")