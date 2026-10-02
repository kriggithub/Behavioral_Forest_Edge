# Behavioral Forest Edge

Analysis code for **a mathematical approach for determining behavioral forest edge effects on mantled howler monkeys (*Alouatta palliata*) in Costa Rica.**

## Questions

1. Does howler monkey behavior (% time resting, moving and feeding) change with distance from an anthropogenic edge or a riparian (river) edge?
2. Does social spacing (number of nearest neighbors and distance to nearest neighbors) change with distance from each edge type?
3. Which functional form (null, linear, power, exponential, logistic, segmented, step/changepoint or unimodal) best describes each response?
4. How far into the forest does each edge effect extend, i.e. what is the depth of edge influence (DEI), and how sensitive is that estimate to the chosen threshold?

## Methods

- **Data preparation** (`dplyr`, `stringr`, `Hmisc`): focal observations are standardized, summarized per monkey, then grouped into 15 m distance bins from each edge, with observation-weighted means and standard errors per bin (`weighted.mean`, `wtd.var`)
- **Model comparison**: bin-level models weighted by the number of monkeys per bin, fitted with `lm` (null, linear), `minpack.lm::nlsLM` (power, exponential, logistic, unimodal), `segmented::segmented` (breakpoint) and `chngpt::chngptm` (step), compared by AIC; Nagelkerke pseudo-R² for the best model via `rcompanion::nagelkerke`
- **Depth of edge influence**: for linear and power models, DEI is the distance at which the fitted response reaches two-thirds of its change across the observed distance range, estimated by inverse prediction (`investr::invest`) with 10,000-replicate nonparametric bootstrap percentile 95% CIs; for logistic models, DEI is the inflection point with a delta-method CI (`msm::deltamethod`); for segmented models, it is the estimated breakpoint with its `confint` interval
- **Threshold sensitivity**: DEI re-estimated at thresholds of 0.10, 0.25, 0.33, 0.50, 0.66, 0.75 and 0.90
- **Figures and tables**: `ggplot2` / `ggpubr` panels; AIC and DEI summary table built with `gt`

## Repository contents

| File | Description |
|---|---|
| `data/BehavioralData.csv` | Raw focal-observation behavioral data |
| `data/monkeyIdData.csv` | Per-monkey summary of activity budgets, nearest neighbors and edge distances |
| `data/anthBinData.csv`, `data/rivBinData.csv` | Binned (15 m) datasets for anthropogenic and riparian edges |
| `R/DataCleaning.R` | Cleans raw data and builds the per-monkey and binned datasets |
| `R/initialPlotting.R` | Exploratory plots of per-monkey and binned data |
| `R/anthModelFitting.R`, `R/rivModelFitting.R` | Fit and compare candidate models for each response and edge type |
| `R/anthDEImodels.R`, `R/rivDEImodels.R` | Best-AIC models with DEI point estimates and 95% CIs |
| `R/anthDEImodelsAdj.R`, `R/rivDEImodelsAdj.R` | Formatted versions of the DEI figures |
| `R/anthDEIvaluePlots.R`, `R/rivDEIvaluePlots.R` | DEI estimates across a range of thresholds |
| `R/tableCreation.R` | AIC, pseudo-R² and DEI summary table |
| `output/initialPlots/`, `output/anthModels/`, `output/rivModels/`, `output/DEIplots/` | Exploratory, model-comparison and DEI figures (PDF) |
| `output/allDEIplotsAnthAdj.pdf`, `output/allDEIplotsRivAdj.pdf`, `output/AICtable.pdf` | Formatted DEI figures and summary table |
| `output/anthDEImodels.RData`, `output/rivDEImodels.RData` | Workspaces saved by the DEI model scripts (generated locally, not tracked) |

## Reproducing the analysis

Open `BehavioralForestEdge.Rproj` in RStudio, install the packages below and run the scripts in this order. All paths are relative to the project root, which RStudio sets as the working directory when the project is opened.

1. `R/DataCleaning.R` (then `R/initialPlotting.R`, which uses the data frames it creates)
2. `R/anthModelFitting.R`, `R/rivModelFitting.R`
3. `R/anthDEImodels.R`, `R/rivDEImodels.R` (or the `*Adj.R` versions for formatted figures)
4. `R/anthDEIvaluePlots.R`, `R/rivDEIvaluePlots.R`
5. `R/tableCreation.R`

The binned CSVs are already included, so steps 2–5 can be run without step 1.

```r
install.packages(c("tidyverse", "Hmisc", "ggpubr", "segmented", "strucchange", "chngpt",
                   "minpack.lm", "rcompanion", "investr", "msm", "knitr", "gt", "webshot2"))
```

## Author

**Kurt Riggin**: [GitHub](https://github.com/kriggithub) · [ORCID](https://orcid.org/0009-0004-4700-1251)
