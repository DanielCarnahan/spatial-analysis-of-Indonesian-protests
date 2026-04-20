# Protest Contagion and the Martyrdom Effect

Daniel Carnahan & Risa Toha

Do protests in one location trigger protests elsewhere, or do protest waves reflect independent responses to common shocks? This project adapts self-exciting point-process logic to discrete district-day panel data to separate background protest risk from cross-district spillovers, and applies it to Indonesian protest activity, 2015–2024.

## The paper

`reproducible_analysis.Rmd` is the single source of truth. Knitting it regenerates the full analysis, all figures, and the PDF/HTML output. It contains:

- **Data construction.** ACLED protest events collapsed into unique protest-location-days, joined to GADM district boundaries, merged with district population, poverty rate, and national CPI.
- **Baseline model (M0).** Quasi-Poisson regression of district-day protest counts on district random intercepts, province-year fixed effects, a national daily protest count (excluding the focal district), and structural covariates, with `log(pop)` as offset.
- **Contagion model (M1).** Adds a normalized cross-district exposure term with spatial (θ) and temporal (κ) decay kernels, estimated via 2-D profile likelihood on a grid.
- **Magnitude.** Counterfactual decomposition giving the endogenous share of predicted protest intensity attributable to spillovers.
- **Martyrdom test.** Separate exposure terms for fatal and non-fatal protests, each with their own optimal decay parameters, plus a Wald test for coefficient equality.
- **Robustness.** Geographic spline, exclusion of the national-shock control, VIF checks, Moran's I on residuals, a logistic specification on a binary outcome, and collapse-window sensitivity (0 / 2 / 3 / 5 / 7 days).

## Main findings

- Cross-district exposure significantly improves fit over the baseline (F-test, p < 0.001).
- Roughly 7% of predicted protest intensity is attributable to diffusion rather than background factors.
- Fatal and non-fatal protests diffuse through different channels: fatal events show national reach consistent with media-driven spread, while non-fatal events diffuse more locally.

## Reproducing

```r
rmarkdown::render("reproducible_analysis.Rmd")
```

The Rmd reads, at the repo root:
- `ACLED Data_2025-10-27_Indonesia20142024_OFFICIAL.xlsx` — protest events
- `district_level_population_2014_2020.csv` — population and poverty
- `indonesia_cpi_fred.csv` — consumer price index
- `Indonesian district-level shapefile copy 2/gadm41_IDN_2.shp` — GADM level-2 districts

The 2-D grid search runs in parallel via `future_lapply`; expect a multi-minute knit on a modern laptop.

**Not in git:** shapefile levels 3 and 4 (`gadm41_IDN_3.*` at 90 MB and `gadm41_IDN_4.*` at 208 MB) exceed GitHub's size limits and are not used by the analysis. Download separately from [GADM](https://gadm.org/download_country.html) if needed.

## Dependencies

R packages: `tidyverse`, `sf`, `readxl`, `lubridate`, `MASS`, `geosphere`, `scales`, `spdep`, `sandwich`, `lmtest`, `mgcv`, `zoo`, `numDeriv`, `knitr`, `kableExtra`, `future`, `future.apply`.

PDF output uses XeLaTeX.

## Repo layout

```
reproducible_analysis.Rmd              # the paper
reproducible_analysis.{pdf,html,tex}   # rendered outputs
presentation.{tex,pdf}                 # conference presentation
presentation_script.md                 # speaker notes
references.bib

ACLED Data_*.xlsx                      # raw protest data
district_level_population_*.csv        # raw population + poverty
indonesia_cpi_fred.csv                 # raw CPI
Indonesian district-level shapefile copy 2/   # GADM

figures/, plots/, report_figures/, visualizations/   # figure outputs
checkpoints_*/                         # model fit state (gitignored)
archive/                               # earlier iteration scripts, drafts, and their outputs
```

## Data sources

- **ACLED** — Armed Conflict Location & Event Data Project. Protest and riot events with geocoordinates and fatality counts. [acleddata.com](https://acleddata.com)
- **World Bank / Indonesia DAPOER** — subnational population and poverty rate.
- **FRED, St. Louis Fed** — Indonesian consumer price index, monthly.
- **GADM 4.1** — district administrative boundaries. [gadm.org](https://gadm.org)
