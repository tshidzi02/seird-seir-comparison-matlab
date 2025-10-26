# SEIR vs SEIRD: Does Explicit Mortality Improve Fit?

This repo contains MATLAB code and results for comparing SEIR and SEIRD models on identical synthetic COVID-like data.

## What’s inside
- code/run_results.m — master script (simulate → build Poisson-noisy observations → fit → export figures)
- code/graphs.m — plots incidence & residuals; returns RMSE/AIC
- code/mypoissrnd.m — Poisson sampler (no Stats Toolbox required)
- code/plot_cumulative.m — cumulative cases comparison
- code/sensitivity_alpha.m — sensitivity sweeps for α (mortality)
- figures/ — exported PNGs for the paper
- results/ — metric tables (RMSE, AIC) and CSVs

## How to run
1. Open MATLAB in the repo root.
2. cd code
3. run_results  
This generates Figures 4.1–4.5 and prints RMSE/AIC to the Command Window.

## Reproducibility
- Synthetic data generation uses a fixed random seed (rng(1)).
- All file paths are relative.
- Large outputs tracked via Git LFS.

## License
MIT — see LICENSE.

## Cite this work
See CITATION.cff. A short form:
> Mutshidzi Madzivhandila 24999709 (2025). SEIR vs SEIRD: MATLAB comparison on synthetic COVID-like data. GitHub repository. (https://github.com/tshidzi02/seird-seir-comparison-matlab.git)
