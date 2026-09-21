# One-time prep script combines the two calibrated chlorophyll spreadsheets into one daily CSV and uploads it to S3
# Not part of the daily run.R workflow. rerun this manually only if new calibrated sheets show up

# run from the repo root: cd BloomBoard && Rscript scripts/combine_calibrated_chla.R

library(readxl)
library(dplyr)

source("R/utils.R")

f1 <- "/projectnb/dietzelab/cwebb16/FRP/Coastal/Calibrated_CHLA_A01_2005-2024.xlsx"
f2 <- "/projectnb/dietzelab/cwebb16/FRP/Coastal/Calibrated_CHLA_A01_2024-2025.xlsx"

d1 <- read_excel(f1, sheet = "DATA") %>%
  select(time = MOORING_TIME_EST, value = CHLA_UG_L)

d2 <- read_excel(f2, sheet = "Export Worksheet") %>%
  mutate(time = as.POSIXct(MOORING_TIME_EST, format = "%m/%d/%Y %I:%M %p", tz = "EST")) %>%
  select(time, value = CHLA_UG_L)

combined <- bind_rows(d1, d2) %>%
  arrange(time) %>%
  distinct(time, .keep_all = TRUE)

# Daily mean plus within-day stats (sd, range, hour-of-day trend) from the hourly data
daily <- combined %>%
  filter(!is.na(value)) %>%
  mutate(date = as.Date(time), hour = lubridate::hour(time)) %>%
  group_by(date) %>%
  # NOTE: must compute value_sd/value_range/value_trend BEFORE value is reassigned to the daily mean below 
  summarise(
    value_sd = if (dplyr::n() >= 2) stats::sd(value) else 0,
    value_range = max(value) - min(value),
    value_trend = if (dplyr::n() >= 2) stats::coef(stats::lm(value ~ hour, data = data.frame(value, hour)))[2] else 0,
    value = mean(value),
    .groups = "drop"
  ) %>%
  filter(!is.nan(value))

message("Combined calibrated chlorophyll: ", nrow(daily), " days, ",
        min(daily$date), " to ", max(daily$date))

cfg <- load_config("config.yml")
write_s3_csv(cfg, cfg$s3$read_path, "a01_calibrated_chlorophyll.csv", daily)
