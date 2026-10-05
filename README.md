# BloomBoard

BloomBoard forecasts chlorophyll up to `max_horizon` days out at a buoy,
using a chlorophyll reading and the buoy's own water sensors as the
starting point. No weather forecast data goes into it.

Check out the live site: [github.com/cami-webb/BloomBoard](https://github.com/cami-webb/BloomBoard).

## How the forecast works

Every time it runs, it checks whether today's chlorophyll is above or
below a bloom threshold (5 µg/L, set in `compute_bloom_threshold()` in
`R/forecast_chlorophyll.R`) and picks one of two paths:

1. **Bloom days** (`forecast_bloom_ens()`): four random forests, trained on
   every past bloom, predict how high the current bloom will still peak,
   how many days it has left, how much of that time is still spent
   rising, and how fast it'll decline once it peaks. None of these are
   recursive, they each just make one prediction per ensemble member.
   Each member's four numbers get built into its own day-by-day
   trajectory, which rises to the predicted peak, then declines at the
   predicted rate (so each path is triangular), and then gets daily noise.

2. **Calm days** (`forecast_calm_ens()`): one random forest trained only on
   calm periods, predicting tomorrow's change in chlorophyll from today's
   value and the water covariates. This one is recursive, so each day's
   prediction feeds into the next day's, using each member's own
   covariate forecasts from `R/forecast_covariates.R` past day one.

Both paths use the same buoy covariates (`covariate_combo` in
`config.yml`). The bloom path also gets a few extra features: days since
the bloom started, how fast it's been rising, day of year, recent
forecast error, each covariate's 7-day trend, and how much chlorophyll
varies within a day.

The forecast is an ensemble of `ensemble: n_members` members (100 in
`config.yml`). Instead of using the forest's average prediction, each
member draws a random quantile level and takes the forest's prediction
at that level (quantile regression forests, `quantreg = TRUE` in
`ranger`), so every member is a different plausible outcome. The spread
across members is the uncertainty and their average is the central
forecast.

Chlorophyll readings aren't available every single day, so whenever a
forecast runs, it seeds from the most recent day that actually has a
reading, not necessarily today (`find_seed_date()` in `R/fetch_data.R`).

## Output format

Each forecast is saved to S3 (`write_path` in `config.yml`) as
`<buoy_id>-<reference date>-<model_id>.csv`, for example
`a01-2026-10-05-bloomboard.csv`, in the EFI standard ensemble format with
one row per member per day (100 members x 7 days = 700 rows):

| Column | Value |
|---|---|
| `project_id` | `bu4cast` |
| `model_id` | `bloomboard` |
| `datetime` | the date being forecast |
| `reference_datetime` | the date the forecast was seeded from |
| `duration` | `P1D` (daily) |
| `site_id` | the buoy's site id (`2` for A01) |
| `family` | `ensemble` |
| `parameter` | ensemble member number (1 to 100) |
| `variable` | `chlorophyll` |
| `prediction` | chlorophyll forecast in µg/L |

The identifiers come from the `efi` block and each buoy's `site_id` in
`config.yml`.

## Covariates

The six water/air conditions the models use, all from A01's sensors
(set in `config.yml`'s `covariate_combo`):

| Covariate | Description | Sensor | Depth | Units |
|---|---|---|---|---|
| `temp_1m` | Water temperature | sbe37 | 1m | °C |
| `temp_20m` | Water temperature | sbe37 | 20m | °C |
| `cond_1m` | Conductivity | sbe37 | 1m | mS/cm |
| `cond_20m` | Conductivity | sbe37 | 20m | mS/cm |
| `sal_20m` | Salinity | sbe37 | 20m | psu |
| `air_temp` | Air temperature | met | -3m (above surface) | °C |

## Model development

Random Forest wasn't the first pick, a handful of model types got tried
before settling on it: ARIMA, XGBoost, a dynamic GAM, and RF. Every
combination of model type and covariate subset was cross-validated
against each other, with about 237,000 combinations total. RF performed the best
with RMSE around 0.27 vs about 0.30 for XGBoost, 0.38
for GAM, and 0.77 for ARIMA. The same was done to choose the most useful covariates: 
water temp, conductivity, salinity, and
air temp kept showing up in the best performing combinations, while
things like mixed-layer depth and current speed did not perform as well.

## Design notes

**Bloom threshold is a fixed 5 µg/L, not dynamic.** A dynamic mean was
tried at one point and came out way too low so calm days kept getting
flagged as blooms.

**Ensemble members sample each forest's predicted distribution.** 
With an ensemble, each member draws a random quantile level (`draw_quantile()` 
in `R/utils.R`), so the full predicted distribution is sampled. The
draws for the different models are independent of each other.

**Bloom ensemble members are shifted upward and given daily noise.**
Sampling each forest's distribution straight underestimates big blooms,
so for the four bloom models the draws are centered on the 0.7 quantile
instead of the median (`bloom_shift_quantile` in `config.yml`). This
keeps the full spread but moves it up. Each member's trajectory then gets
multiplicative daily noise (`bloom_noise_sd`) that grows from about 20%
on day 1 to about 40% on day 7, since real bloom days swing around more
than a smooth rise and decline can. The two noise values were tuned in
the backtest so the 80% interval covers about 80% of observations.

**The rise phase drops to 0 days if there's barely any rise left.**
Otherwise a bloom that's already at or near its peak would still get a
couple of "rise days" tacked on and plateau before declining, instead of
just declining right away.

**`recent_forecast_error`** is basically a "how surprised should the
model be right now" signal. It's basically today's observed chlorophyll minus what a naive
guess (yesterday's value plus yesterday's trend) would've predicted. A
version using the model's actual forecast from a day ago was tried
instead, but it broke the parallel backtesting and wasn't worth the
slowdown, so it was reverted to the naive guess.

**Bloom training panel columns** (`compute_bloom_training_panel()`, one
row per day inside a bloom):
- `days_since_onset` - day 1 is the onset day
- `rise_rate_since_onset` - average daily change since the bloom started
- `remaining_peak` / `remaining_duration` - what `peak_model` and
  `duration_model` are trying to predict
- `doy` - day of year, so spring/fall blooms (usually bigger, longer)
  aren't treated the same as short summer ones
- `<cov>_trend7` - each covariate's change over the last week
- `decline_rate` - this bloom's own peak-to-end decline rate, so the
  model can learn some blooms crash faster than others instead of using
  one average rate for everything; NA if the peak was the last day
  observed (no decline ever seen)
- `rise_fraction` - this bloom's own fraction of its duration spent
  rising, the target for `rise_fraction_model`

**`target_mean`/`oob_bias` attributes.** Each model's training target
mean and average out-of-bag error get attached directly onto the fitted
model as attributes, so they travel with the model instead of
being a separate thing that could get mismatched. `dashboard/backtest.qmd`
uses these for the rBias/rRMSE numbers in the model fit quality table.

**`lookup_with_fallback()`** (`R/utils.R`) looks for a covariate on the
exact date first, then falls back to the closest value within a few days,
then to that day-of-year's historical average. Buoy sensors go down
sometimes, and a strict same day lookup was turning otherwise fine days
into all-NA forecasts.

**Hourly chlorophyll stats** (`chlora_hourly_sd/range/trend`, from
`load_sensor_daily()` and `scripts/combine_calibrated_chla.R`) let the
bloom model tell a stable day apart from a volatile one, instead of only
seeing the flattened daily average.

## Dashboard site

`dashboard/` is a Quarto site published on GitHub Pages. Pages:
- `index.qmd` - "Today": current bloom status, 7-day outlook, monthly and
  annual trend charts, water column outlook and stratification charts, and a
  currents map, per buoy.
- `about.qmd` - "About": the A01 buoy, who uses it, and why.
- `backtest.qmd` - "Backtest Validations": the 2021-2024 offline backtest,
  accuracy tables, out-of-bag model fit quality.
- `documentation.qmd` - "Documentation": covariates, trajectory features,
  the forecast pipeline diagram, and the model development story (most of
  what's summarized above).
- `acknowledgements.qmd` - funding and data source credits.
- `_quarto.yml` - shared site config and CSS (navbar, fonts, layout).
- `img/` - navbar logos.

Three GitHub Actions workflows run it:
- `daily.yaml` - runs `run_ensemble.R` every day at 17:00 UTC, plus a manual trigger. Pings a
  healthchecks.io check-in so a failed or missed run gets flagged.
- `dashboard.yaml` - renders and publishes `index.qmd`, `documentation.qmd`,
  and `acknowledgements.qmd`, daily at 17:30 UTC and on manual trigger.
- `backtest-page.yaml` - renders and publishes `backtest.qmd` on its own,
  manual trigger only, since it pulls ~1500 historical files and doesn't
  actually change day to day.

Both publishing workflows push to the `gh-pages` branch, and each one only
touches the files it owns so they don't overwrite each other.

## Folders and files

- `config.yml` - buoys, their coordinates, sensors, which covariates to
  actually use (`covariate_combo`), and run settings (`max_horizon`,
  `run_mode`). Add a new buoy by adding a key under `buoys`, add a
  covariate by adding a key under that buoy's `covariates`. Each sensor
  has a `historical_url` and `realtime_url` (GoMOOS/NERACOOS file paths)
  and a `use_realtime` toggle. The chlorophyll entry also has a
  `calibrated_url` pointing at the QC'd data (see
  `scripts/combine_calibrated_chla.R`), used instead of the live feed
  whenever `run_mode` is `date_range`. Also holds the ensemble settings 
  (`ensemble: n_members`, `bloom_shift_quantile`, `bloom_noise_sd`), the `efi`
  block (project_id, model_id, duration, variable), and a `site_id` for each buoy.
- `R/utils.R` - config loading, S3 read/write, sensor fetching, the
  gap-tolerant covariate lookup. Also has `draw_quantile()` (random draws from 
  a quantile forest) and `efi_ensemble_rows()` (formats a forecast in EFI format).
- `R/fetch_data.R` - pulls chlorophyll and covariates for a buoy and joins
  them into one daily table.
- `R/forecast_covariates.R` - trains and forecasts each water covariate forward 
as an ensemble (used by the calm path).
- `R/forecast_chlorophyll.R` - the bloom/calm dispatch and both
  forecasting paths described above.
- `run_ensemble.R` - runs everything for every configured buoy: fetch data,
  save it, forecast the ensemble, save it in EFI format. `run_mode` in the
  config decides if it's a `daily` run or a `date_range` backtest, same
  code path either way. A backtest range can also be passed directly,
  e.g. `Rscript run_ensemble.R 2021-01-01 2021-12-31`. Dates that already
  have a forecast on S3 are skipped, and dates run in parallel via
  `foreach`/`doParallel` when there's more than one core available.
- `run_ensemble.qsub` - SCC array job that runs the backtest with one task
  per year (2021 to 2024), in parallel.
- `dashboard/` - the Quarto site, see above.
- `.github/workflows/` - the three workflows, see above.
- `scripts/combine_calibrated_chla.R` - one-time script, not part of the
  daily run. Combines the two calibrated chlorophyll spreadsheets (collected and QC'd by
  the MWRA) into one daily CSV with
  within-day stats, and uploads it to S3. Only needs rerunning if new
  calibrated sheets show up. `scripts/combine_calib.qsub` is its batch job.
- `validation/` - SCC-only, not pushed to GitHub. Has
  `validate_backtest.Rmd`, an older backtest report (accuracy by
  year/lead time/branch, variable importance, OOB fit) plus its batch job
  and output. Mostly replaced now by `dashboard/backtest.qmd` for anything
  public-facing.
- `logs/` - SCC-only, batch job logs.
