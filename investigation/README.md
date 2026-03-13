# Investigation: Post-1990 Cyclone-GDP Divergence

This folder contains all analysis scripts, figures, and slides investigating why post-1990 cyclone shocks produce larger GDP losses than pre-1990 shocks in the local projection framework.

## Key Finding

The post-1990 divergence is largely explained by **heterogeneous regional GDP trends correlated with treatment status**. 142 of 182 countries never experience a p95 cyclone shock. These never-treated countries (concentrated in Europe and Africa) dominate the year fixed effects, defining the implicit counterfactual. When region-specific year FE are used instead of global year FE, the pre/post divergence disappears.

## Slides

- `lp_gdp_shocks_slides.tex` -- Beamer presentation summarizing all findings. Compile with `pdflatex`.

## Analysis Scripts

All scripts are run from the **repo root** with `Rscript investigation/<script>.R`. They depend on `emileRegs.R` and data in `raw_data/`.

| Script | Description |
|--------|-------------|
| `investigation_composition.R` | Panel balance and sample composition checks |
| `investigation_measurement.R` | Cyclone measurement artifacts (satellite era, wind distributions) |
| `investigation_shocks.R` | Shock frequency, clustering, and timeline |
| `investigation_breakpoint.R` | Breakpoint robustness (rolling split, decade interactions) |
| `investigation_vulnerability.R` | Economic vulnerability, LOO by top-exposed countries |
| `investigation_influential.R` | Leave-one-year-out and leave-one-country-out influence |
| `investigation_income_het.R` | Income quartile/quintile heterogeneity, exposure composition |
| `investigation_structural.R` | Capital intensity, time trend, and human capital interactions |
| `investigation_globalization.R` | Trade openness and investment share (full PWT 10.01) |
| `investigation_never_treated.R` | Drop never-treated countries robustness check |
| `investigation_never_treated_deep.R` | Characterize never-treated vs treated: regions, GDP trends, regional exclusion |
| `investigation_counterfactual.R` | Year FE decomposition, placebo test, continent/region x year FE |

## Figures

All figures are saved to `figures/` and referenced in the slides. See the slides for descriptions of each figure.
