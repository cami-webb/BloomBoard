# BloomBoard

Forecasts chlorophyll up to `max_horizon` days ahead at a buoy, seeded from
a chlorophyll observation and the buoy's water condition sensors. 
No weather forecast covariates.

## How the forecast works

Each reference date is routed to one of two branches depending on whether
today's real chlorophyll is above or below a fixed bloom threshold
(5 µg/L, `compute_bloom_threshold()` in `R/forecast_chlorophyll.R`):

1. **Bloom days** (`forecast_bloom()`): three random forests, each trained
   once on all historical bloom-period data, predict how high the
   remaining bloom will peak, how many more days it will last, and its
   post-peak decline rate. All three are single, non-recursive
   predictions.The day-by-day trajectory is then built geometrically from those 
   three numbers (rise linearly to the predicted peak, then decline linearly at
   the predicted rate), not by stepping the model forward one day at a
   time.

2. **Calm days** (`forecast_calm()`): a single random forest trained only
   on calm-period data, predicting tomorrow's change in chlorophyll from
   today's state and water covariates. This one is recursive,
   each day's prediction becomes the next day's starting state, using
   `R/forecast_covariates.R`'s own recursively-forecasted covariate values
   for days beyond the first.

Both branches use a shared set of buoy water covariates (`covariate_combo`
in `config.yml`), plus, for the bloom branch only, features
like days since onset, rise rate since onset, day of year, each
covariate's 7-day trend, and within-day chlorophyll variability computed
from the raw hourly readings.

Chlorophyll isn't read every single day, so the reference date used to
seed a forecast is the most recent date at or before the run date with a
real chlorophyll reading (`find_seed_date()` in `R/fetch_data.R`), not
necessarily the run date itself.

## Design notes

**Bloom threshold is fixed at 5 µg/L, not dynamic.** A dynamic full-record
mean was tried and reverted as it came out way too low and calm days were being
classified as blooms

**The bloom sub-models predict an upper quantile (0.7), not the mean, for
peak magnitude and duration.** Both populations are right-skewed w/ most
events being short/mild and a few long/large), so the mean prediction
underestimates actual large events. `decline_rate` also uses the 0.7 quantile, 
but since it's negative, the upper quantile means a milder
decline rather than a more extreme one; this was used to counteract
declines that were going down too fast.

**The rise phase collapses to 0 days if there's little predicted further
rise.** Without this, a bloom already near its predicted peak would still
get 1-3 nominal "rise days" (sized off the predicted duration, not the
actual remaining gap), plateauing before the decline
instead of declining when it should decline.

**`rise_fraction`** is the empirical median fraction of an event's
duration that elapses before its peak, computed across every historical
bloom event. Replaces a flat "peak lands halfway through the predicted
duration" assumption with the shape of historical blooms.

**`recent_forecast_error`** (in `compute_bloom_training_panel()`) is a
naive-persistence "surprise" signal: today's actual chlorophyll minus a
naive extrapolation of yesterday's trend (yesterday's value + yesterday's
own day-over-day change). A positive value means it kept rising when a simple 
extrapolation would have expected a slowdown. A version using what the model actually
forecast for today, one day ago, was tried instead and reverted bc
no parallelism was possible (not worth the backtest slowdown). 

**Bloom training panel columns** (`compute_bloom_training_panel()`, one
row per day inside a detected bloom event):
- `days_since_onset` - 1 = onset day
- `rise_rate_since_onset` - average daily change since the bloom started
- `remaining_peak` / `remaining_duration` - the "how high from here" /
  "how much longer" training targets for `peak_model`/`duration_model`
- `doy` - day of year, so the model can learn that spring/fall-onset
  blooms behave differently (often longer, larger) than short summer
  pulses, instead of predicting the same "typical" duration regardless of
  season
- `<cov>_trend7` - each covariate's change over the last 7 days
- `decline_rate` - this event's own peak-to-end decline rate (µg/L/day,
  negative), assigned to every row in the event, so `decline_rate_model`
  can learn that some blooms crash faster or slower than others instead of
  using one fixed historical average rate for every forecast; NA for
  events where the peak was the last observed day (no decline phase was
  ever seen)

**`target_mean` attribute on fitted models.** `train_bloom_models()` and
`train_calm_model()` attach each fitted model's own training target mean
as an R attribute (not an extra return value), so it rides along with the
model object itself and can't be mismatched against the wrong one. Used
by `validation/validate_backtest.Rmd` to report relative RMSE.

**`lookup_with_fallback()`** (`R/utils.R`) looks up a covariate on a
target date, falling back to the most recent value within `max_lookback`
days, then that day-of-year's historical climatological mean. Buoy
sensors have gaps (equipment downtime, comms outages, sometimes a
week-plus), so a strict same-day lookup on a single covariate was
producing entirely NA forecasts on otherwise good days

**Within-day chlorophyll stats** (`chlora_hourly_sd/range/trend`,
computed in `load_sensor_daily()` and `scripts/combine_calibrated_chla.R`)
let the bloom model see whether a day's mean chlorophyll came from a
stable day or one with a lot of inter-daily variability, not just the
collapsed daily average.

## Folders and files

- `config.yml` - buoys, their coordinates, their chlorophyll and covariate
  sensors, which covariates to use (`covariate_combo`), and run
  settings (`max_horizon`, `run_mode`). Add a buoy by adding a key under
  `buoys`. Add a covariate by adding a key under that buoy's `covariates`.
  Each sensor entry has a `historical_url` and `realtime_url`
  (GoMOOS/NERACOOS data portal file paths) and a `use_realtime` toggle for
  whether to include the daily-updating feed or just the historical
  archive. The chlorophyll entry also has a `calibrated_url`, an S3 path
  to QC'd calibrated chlorophyll data (see `scripts/combine_calibrated_chla.R`
  below) used instead of the live portal feed whenever `run_mode` is
  `date_range`.
- `R/utils.R` - config loading, S3 read/write, sensor data fetching,
  gap-tolerant covariate lookups.
- `R/fetch_data.R` - pulls chlorophyll and covariate data for a buoy and
  joins it into one daily table. Chlorophyll comes from the calibrated S3
  file in `date_range` mode, or the live portal feed in `daily` mode.
- `R/forecast_covariates.R` - trains and recursively forecasts each water
  covariate forward (used by the calm branch).
- `R/forecast_chlorophyll.R` - the bloom/calm dispatch and both
  forecasting branches described above.
- `run.R` - runs the whole workflow for every configured buoy: fetch data,
  save training data, forecast covariates, forecast chlorophyll, save
  forecast output. Controlled by `run_mode` in the config: `daily` runs
  for today, `date_range` runs a backtest over a date range instead. Same
  code path either way. Backtests run in parallel across dates via
  `foreach`/`doParallel` when `NSLOTS` (set by `qsub -pe omp N`) is greater
  than 1.
- `run.qsub` - SCC batch job for `run.R` (16 cores, for backtests).
- `.github/workflows/daily.yaml` - runs `run.R` every day and can also be
  triggered manually. Pings a healthchecks.io check-in after a successful
  run so a missed or failed daily run gets flagged.
- `scripts/combine_calibrated_chla.R` - one-time prep script, not part of
  the daily workflow. Combines the two calibrated chlorophyll spreadsheets
  (QC'd directly by the people who collected the data, covering deployments
  the live GoMOOS portal splits across separate per-deployment archives)
  into one daily CSV plus within-day stats, and uploads it to S3 at
  `calibrated_url`. Rerun manually only if new calibrated sheets show up.
  `scripts/combine_calib.qsub` is its SCC batch job.
- `validation/` - SCC-only, gitignored (not pushed to GitHub). Holds
  `validate_backtest.Rmd`, which knits a backtest validation report
  (accuracy by year/lead time/branch, bloom event accuracy, variable
  importance, out-of-bag model fit quality) from whatever forecast output
  is currently on S3 for the buoy/date range in `config.yml`, plus its SCC
  batch job `knit_validation.qsub`, and the report's own output (HTML,
  CSVs, PNGs) once run. Not part of the daily workflow -- run manually
  after a backtest.
  
*Readme generated/organized with help from Claude :)
