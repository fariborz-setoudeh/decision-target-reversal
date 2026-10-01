# Decision-Target Reversal under Informative Label Missingness

Code and analysis inputs for the manuscript **“Decision-Target Reversal under Informative Label Missingness.”**

The study examines how changing deployment misclassification costs can reverse the risk ordering between complete and informative partial classification.

## Files

| File | Purpose |
|---|---|
| `01_population_theory.R` | Population information matrices, phase map and theoretical critical costs |
| `02_symmetric_simulation.R` | Main simulation: 5,000 paired replications at each sample size |
| `04_unequal_prior_simulation.R` | Unequal-prior simulation: 3,000 paired replications |
| `07_gastro_common_covariance_FINAL.R` | Gastrointestinal-lesion analysis with training-only leave-one-out fitting |
| `03_reproduce_paper_figures.py` | Paper figures, Monte Carlo uncertainty, bootstrap crossing summaries and second-moment checks |
| `gastro_decision_target_analysis.rds` | Four selected features, histological reference classes and agreement-defined label mask |

## Requirements

R packages:

```r
install.packages(c("ggplot2", "dplyr", "readr", "scales"))
```

Python packages:

```bash
python -m pip install numpy pandas scipy matplotlib rdata
```

The real-data script uses base R only.

## Reproduction

Run these commands from the repository folder:

```bash
Rscript 01_population_theory.R .
Rscript 02_symmetric_simulation.R .
Rscript 04_unequal_prior_simulation.R .
Rscript 07_gastro_common_covariance_FINAL.R . gastro_decision_target_analysis.rds
python 03_reproduce_paper_figures.py .
```

The simulations can take substantial time. Script 02 resumes compatible saved replications. Script 07 rechecks compatible saved fits; to test fitting from scratch, run it in a new folder containing only Script 07 and the input RDS.

Outputs are written to:

- `decision_target_reversal_results/`
- `review_gastro_common_covariance/`
- `paper_figures/`

The figure script reads generated results and does not fit models or generate new Monte Carlo samples.

## Real-data source and interpretation

The original gastrointestinal-lesion dataset is available from the [UCI Machine Learning Repository](https://doi.org/10.24432/C5V02D).

The supplied input retains the four previously selected features and agreement-defined label mask. The analysis is conditional on that historical feature selection.

The real-data comparison is exploratory. Numerical checks establish reproducible local optima, not guaranteed global maxima. Observed risk differences do not establish a population predictive advantage, and deployment cost ratios are sensitivity parameters rather than clinical recommendations.
