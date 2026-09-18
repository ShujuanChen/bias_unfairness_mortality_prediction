# Selective health data drive bias and unfairness in mortality prediction

This repository contains the code for the primary analysis and tutorial of implementing it for the study *Selective health data drive bias and unfairness in mortality prediction*. The overall study design is shown in Figure 1. Please refer to the [Overview](#overview) section below for the layout of this repository. The two step-by-step tutorials, for reproducing the primary study results and for generating the participation weights, are published at [shujuanchen.github.io/bias_unfairness_mortality_prediction](https://shujuanchen.github.io/bias_unfairness_mortality_prediction/) and run on synthetic data without any restricted-data access.

![Overall study design](Schematic_study_design.png)

<sub>To dissect prediction bias arising from selective training data, we designed a controlled benchmark framework in which the training dataset varied while the modelling framework, predictors, outcomes, and evaluation population were held constant. Mortality models were trained on selective health data from UK Biobank (UKB) and, separately, on the nationally representative benchmark dataset, the personal-level mortality registrations (PMR) linked with census. Both models were then deployed to the same target population, PMR, to predict mortality risk. Prediction bias was estimated by comparing predicted risk with observed mortality. PMR-trained model served as a representative benchmark to verify whether the same modelling framework produced well-calibrated predictions when trained on representative data.  To assess bias correction, UKB was reweighted to better resemble the wider target population using two additional representative reference datasets. Super Learner was used to estimate participation weights.  Sensitivity analyses assessed whether the main findings were robust to alternative prediction models, reference datasets, and weighting methods.</sub>

## Overview

The [two step-by-step tutorials](https://shujuanchen.github.io/bias_unfairness_mortality_prediction/) use synthetic data, so neither requires access to restricted datasets. Tutorial 1 shows how to reproduce the primary study analysis. Tutorial 2 focuses on generating the participation weights for UK Biobank sample. This participation weighting framework can potentially be applied to other similar settings of selective biobank databases. These weights can be used to handle the covariate shift between the selective health data and their target population, conditioned on their shared variables. Both tutorials identify the scripts used at each step. This README summarises the pipeline, configuration, outputs, runtime, software requirements and data access.

- [Pipeline](#pipeline)
- [Configuration](#configuration)
- [Outputs](#outputs)
- [Run time](#run-time)
- [System requirements](#system-requirements)
- [Installation](#installation)
- [Data and access](#data-and-access)
- [Repository layout](#repository-layout)
- [Licence](#licence)

## Pipeline

The analysis has five phases: harmonisation, prediction-bias assessment, participation weighting, bias correction and uncertainty estimation. [run_pipeline.R](run_pipeline.R) is the single entry point and calls the script responsible for each step. The [tutorial](https://shujuanchen.github.io/bias_unfairness_mortality_prediction/) provides the commands and explains the scripts used in each phase.

| Phase             | What it does                                                                                                                   |
| ----------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| 1 Harmonisation   | Harmonises all four sources onto shared coding schemes, applies the row rules in one place and defines the analytic population |
| 2 Prediction bias | Trains models on UKB and PMR, evaluates both in PMR and reports the prediction gap; this phase runs without phase 3            |
| 3 Weighting       | Estimates participation weights using the method specified by`phase3.weight_method` and reports diagnostics                  |
| 4 Bias correction | Retrains the outcome models on reweighted UKB and reports how much of the gap is reduced                                       |
| 5 Uncertainty     | Runs multiplier-bootstrap replicates of the pipeline and calculates uncertainty of the full pipeline modelling                 |

`--pipelines` specifies whether the analytic population is built for `outcome`, `hse`, `census` or a combination of them. The default includes both reference datasets, as in the study. Users can select only one reference sample to generate weights without the possession of the other source; stages requiring the other omitted source are skipped. A single-reference population may be larger than the intersection used in the study because it does not apply inclusion rules for the other pipelines. 

A phase 5 replicate is generated by [multiplier_bootstrap.py](code/phase5_uncertainty/multiplier_bootstrap.py). It applies strictly positive, mean-one multipliers to cohort participants, target-population records within design strata and reference-sample sampling units, without dropping rows. The participation and outcome models are refitted in each replicate, so the intervals include the uncertainty from full-pipeline modelling including the weight-estimation uncertainty. Each replicate has its own intermediate directory at `temp/<run>/replicates/<method>_NNNN/`. It links to the point-estimate harmonisation, population and folds, and stores its own refitted outputs. Replicates do not write directly to the reported results.

To run the pipeline, generate the synthetic inputs, update any required paths in [framework_config.json](framework_config.json), and run the five phases in order.

## Configuration

[framework_config.json](framework_config.json) specifies the analysis settings, grouped by phase.

| Block      | What it holds                                                                                                                                                                                                       |
| ---------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `phase1` | The study window, the age range, the two horizons and the three predictor sets                                                                                                                                      |
| `phase2` | The device every fit runs on, the outcome definitions, the five weight-source keys from different weighting approaches, the whole network specification, the splits and their seed                                  |
| `phase3` | The weighting method this run uses, the five it could use, the cross-fitting folds and seeds, the winsorisation percentiles, the participation learners and their hyperparameters, and the moment-matching settings |
| `phase5` | The replicate count, the interval level and the draw's seed                                                                                                                                                         |

Phase 4 has no separate configuration block. It uses the weights produced in phase 3; figure settings are defined in the reporting stage.

`phase3.weight_method` selects the weighting method. The other settings describe the reported analysis, and the code stops if it encounters a value it cannot implement. Five weighting methods (including alternative reference dataset) are available:

| `phase3.weight_method`  | What it does                                                                                                                        |
| ------------------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| `hse_superlearner`      | Cross-fitted Super Learner ensemble of five learners, against HSE, to generate inverse-odds weight.                                |
| `hse_lassologit`        | Lasso-penalised logistic regression with second-order interactions, against HSE, on the same outer folds                            |
| `hse_raking`            | Iterative proportional fitting onto the HSE margins, continuous auxiliaries binned at the reference sample's own weighted quintiles |
| `hse_entropy_balancing` | Minimum-entropy weights matching the first two moments of every continuous auxiliary and every categorical share, crossed with sex  |
| `census_superlearner`   | The Super Learner again, against the 2011  Census microdata as reference dataset                                                  |

Changing this key selects the method used for weight estimation, diagnostics, weighted outcome fits, correction tables and the bootstrap. Model-based methods estimate the probability of cohort membership in the pooled samples and use the winsorisation percentiles specified by `phase3.winsorise_percentiles` to handle the extreme weights. Raking and entropy balancing do not fit membership probabilities and are not winsorised, because clipping their weights would violate the balance constraints used to derive them.

Raking and entropy balancing iterate until each margin meets its configured tolerance or `phase3.moment_matching.max_iterations` is reached. On the real cohort, they met their tolerances after 115 and 7 iterations, respectively. If a solver stops short, the stage warns with the remaining gap and tolerance but still writes the weights. Such a warning matters because weights that do not meet the constraints may not represent the intended estimator. On the synthetic inputs, raking meets its tolerance while entropy balancing stops at the iteration limit, 5,000 iterations with a largest remaining gap of 0.00111 against a tolerance of 1e-08. Raking also warns that one binned auxiliary has no cohort record in its top bin, so no weighting can move that share onto its target. Both messages are consequences of the synthetic columns being drawn independently, and neither occurs on the real cohort.

## Outputs

`results/<run>/` contains one directory per phase, with its figures and tables together. Intermediate files, including tables used only to draw figures, are stored under `temp/<run>/`.

Each new run receives a timestamped identifier, keeping its files separate from previous runs. `--run-id` continues an existing run so later phases can use earlier outputs. Each invocation also writes a manifest under `temp/<run>/manifest/` that records its commands, exit statuses, runtimes, R version, platform, Python interpreter and computing device.

The row counts below are expected for the synthetic inputs. Most are set by the analysis structure rather than the generated data, so differences may indicate that a run used different settings or inputs.

### `results/phase1_harmonisation/`

| Output              | Rows | What it is                                 |
| ------------------- | ---- | ------------------------------------------ |
| `flow_table.xlsx` | 39   | Population-flow steps for the four sources |

### `results/phase2_prediction_bias/`

| Output                                                       | What it is                                                                                                   |
| ------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------ |
| `prediction_bias_by_cause_{5,10}y.png`                     | Prediction bias in mortality risk for UKB- and PMR-trained models, by cause                                 |
| `prediction_bias_by_strata_{5,10}y.png`                    | All-cause prediction bias by subgroups                                                                       |
| `additional_unflagged_death_unweighted_{5,10}y.png`        | Additional observed deaths left unflagged by the UKB-trained model relative to PMR, across risk thresholds |
| `difference_in_false_negative_rate_unweighted_{5,10}y.png` | Difference in false-negative rate between models across risk thresholds                                      |
| `calibration_predictive_parity.xlsx`                       | Comparison of between-subgroup variation in relative bias and positive predictive value                      |
| `model_performance.xlsx`                                   | Weighted AUC at both horizons, concordance index and Brier score for the UKB- and PMR-trained models         |

### `results/phase3_weighting/`

| Output                      | What it is                                                                |
| --------------------------- | ------------------------------------------------------------------------- |
| `weight_overlap.png`      | Fitted participation-probability distributions with common-support bounds |
| `weight_distribution.png` | Weight distribution, with percentile cut points                           |
| `weight_summary.xlsx`     | Weight distribution and effective sample size                             |
| `weight_positivity.xlsx`  | Overlap bounds, proportions outside common support and inverse-odds range |

### `results/phase4_bias_correction/`

| Output                                     | What it is                                                                    |
| ------------------------------------------ | ----------------------------------------------------------------------------- |
| `bias_correction_by_cause_{5,10}y.png`   | Bias before and after weighting by cause                                      |
| `bias_correction_by_strata_{5,10}y.png`  | All-cause bias correction by subgroup                                         |
| `additional_unflagged_death_{5,10}y.png` | Additional unflagged deaths across risk thresholds before and after weighting |
| `model_performance.xlsx`                 | Weighted AUC at both horizons, concordance index and Brier score before and after weighting |

### `results/phase5_uncertainty/`

| Output                                   | Rows | What it is                                                        |
| ---------------------------------------- | ---- | ----------------------------------------------------------------- |
| `bootstrap_intervals_by_cause.xlsx`    | 10   | Percentile intervals for bias and correction by cause and horizon |
| `bootstrap_intervals_by_subgroup.xlsx` | 68   | All-cause intervals by horizon and subgroup level                 |

These two workbooks contain all reported intervals, at the level specified by `phase5.interval_level`. The replicate count has to support the level: at 0.95 each tail needs 80 replicates, so a demonstration run of 10 leaves 0.3 replicates in each tail and returns the extreme replicates instead of estimated quantiles. The stage warns when this happens. Each row gives the point estimate, bootstrap mean and unadjusted percentile bounds. Once phase 5 is complete, it redraws the phase 4 correction figures with intervals. 

## Run time

The following runtimes were measured on synthetic inputs using the reported method, with `--jobs 6` for phases 2 and 4 and five Super Learner fold workers, each using four threads. “Command seconds” sums the time spent on individual commands; “wall clock” is the elapsed time for each phase.

| Phase                        | Commands     | Command seconds | Wall clock         |
| ---------------------------- | ------------ | --------------- | ------------------ |
| 1 Harmonisation              | 6            | 4               | 4 s                |
| 2 Prediction bias            | 50           | 153             | 31 s               |
| 3 Weighting                  | 3            | 70              | 1.2 min            |
| 4 Bias correction            | 14           | 37              | 11 s               |
| **Phases 1 to 4**      | **73** | **264**   | **1.9 min**  |
| 5 Uncertainty, 10 replicates | 2            | 855             | 14.2 min           |
| **Phases 1 to 5**      | **75** | **1,119** | **16.2 min** |

Phase 5 dominates the total and scales with the replicate count, so it is also given per replicate. Each replicate refits the participation model and every outcome model, one participation model and 42 outcome fits here, and `--all` runs those fits two at a time. On these synthetic inputs a replicate takes 85 seconds with the reported method.

The replicate cost is what separates the five weighting methods, because they share the same outcome fits and differ only in the participation model each replicate re-estimates. Run end to end at ten replicates on the same inputs and machine, `hse_raking` takes 4.4 minutes and 21 seconds per replicate, `hse_entropy_balancing` 5.5 minutes and 27 seconds, `census_superlearner` 15.2 minutes and 80 seconds, `hse_superlearner` 16.2 minutes and 85 seconds, and `hse_lassologit` 49.6 minutes and 277 seconds, the second-order interaction design making it the most expensive of the five.

## System requirements

The analysis code is written in R and Python, with no compiled components in the repository. It was tested using R 4.5.3 and Python 3.10.20. Package requirements and versions are listed in [R_requirements.txt](R_requirements.txt) and [requirements.txt](requirements.txt).

The pipeline reads three environment variables:

| Variable                  | What it is for                                                                                                |
| ------------------------- | ------------------------------------------------------------------------------------------------------------- |
| `REWEIGHTING_PYTHON`    | Python interpreter containing the packages in`requirements.txt`; phases 2, 4 and 5 check it before starting |
| `OPENBLAS_NUM_THREADS`  | Set to 1 to avoid competing threads in the linear-algebra library                                             |
| `REWEIGHTING_SL_MIN_GB` | Lowers the Super Learner memory requirement for synthetic demonstrations                                      |

The outcome-model device is set by `phase2.device` in [framework_config.json](framework_config.json). It accepts `gpu` or `cpu` and defaults to `gpu`. On a machine without a CUDA device, set it to `cpu` before running the outcome-model phases.

CPU training is single-threaded with a fixed summation order, supporting repeatable fits on the same machine. Fits may not be bit-identical across platforms or devices. Participation models run on the CPU regardless of the outcome-model device.

## Installation

```
# R, exact versions for reproduction, for example
#   remotes::install_version("survival", "3.8.3")
# or the current versions
install.packages(c("survival","SuperLearner","glmnet","ranger","xgboost","dbarts","nnet",
  "fastDummies","Matrix","dplyr","stringr","lubridate","tibble","forcats","plyr",
  "data.table","haven","readxl","writexl","jsonlite"))

# Python
python -m pip install -r requirements.txt
```

Installation typically takes 20–40 minutes on a desktop computer, mostly for PyTorch and the R machine-learning packages.

## Data and access

Restricted individual-level study data are not distributed in the repository. The `data/` folders show the expected layout and can hold the synthetic inputs generated by the tutorials.

## Repository layout

```
code_publish/
├── README.md                          # this file
├── LICENSE                            # MIT licence
├── framework_config.json              # every setting the analysis runs on
├── run_pipeline.R                     # the only entry point
├── R_requirements.txt                 # R dependencies, with versions
├── requirements.txt                   # Python dependencies, with versions
├── Schematic_study_design.png         # Figure 1 study-design schematic
├── docs/
│   ├── index.html                     # rendered online tutorial (GitHub Pages, /docs)
│   └── .nojekyll                      # serve the HTML as-is (skip Jekyll)
│
├── code/
│   ├── common/                        # shared library, read by every phase
│   ├── phase1_data/                   # harmonisation, row rules, population, record flow
│   ├── phase2_prediction_bias/        # folds, training, assembly, cumulative incidence, results
│   ├── phase3_weighting/              # the five weighting methods, diagnostics and results
│   ├── phase4_correction/             # correction results
│   └── phase5_uncertainty/            # multiplier bootstrap and intervals
│
├── data/                              # place each restricted input here, or generate a stand-in
│   ├── regenerate_synthetic_inputs.sh # write every synthetic input
│   ├── UKB/                           # the UKB extract, and its generator
│   ├── PMR/                           # the PMR extract, and its generator
│   ├── HSE/                           # the HSE survey waves, and their generator
│   ├── census/                        # the Census microdata, and its generator
│   ├── imd/                           # public deprivation scores, and their generator
│   └── lookup/                        # public geography lookups, and their generator
│
├── temp/<run>/                        # intermediates, manifests and replicate trees (git-ignored)
└── results/<run>/                     # one directory per phase (git-ignored)
```

## Licence

This repository is released under the MIT licence. The full text is in [LICENSE](LICENSE).
