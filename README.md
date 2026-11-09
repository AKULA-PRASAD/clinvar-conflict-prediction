# Predicting conflicting clinical interpretations of genetic variants (ClinVar)

A machine-learning pipeline in R that predicts whether a ClinVar variant has **conflicting clinical interpretations** between submitting laboratories, from 65,188 variants in the public [Kaggle ClinVar-conflicting dataset](https://www.kaggle.com/datasets/kevinarvai/clinvar-conflicting) (25.2% conflicting). It began as my DA5030 (Introduction to Data Mining and Machine Learning, Northeastern University) signature project in Spring 2026. It was revised in November 2026 to fix leakage in cross-validation, implement genuine out-of-fold stacking, and test how much performance depends on recognising genes seen in training.

**Report:** https://akula-prasad.github.io/clinvar-conflict-ml/ (rendered from [`analysis.Rmd`](analysis.Rmd))

## Results (held-out 30% test set)

| Model | 5-fold CV AUC | Test AUC | Test AUPRC | F1 |
|---|---|---|---|---|
| Logistic regression | 0.676 | 0.669 | 0.379 | 45.8% |
| Decision tree (CART) | 0.670 | 0.699 | 0.399 | 50.0% |
| **Random forest** (best by CV AUC) | **0.713** | **0.734** | **0.454** | **51.0%** |
| Bagged trees | — | 0.701 | 0.402 | 42.3% |
| Majority-vote ensemble | — | — | — | 50.2% |
| Stacked ensemble (out-of-fold logistic meta-learner) | — | 0.717 | 0.430 | 49.5% |

- **Reference points:** a random ranking gives AUPRC 0.252 (the prevalence). Always predicting "Concordant" gives 74.8% accuracy. Every model is below that accuracy, because they are trained on balanced classes and trade accuracy for recall.
- **Gene-grouped check:** random-fold CV AUC is 0.715, but **0.684 when every gene's variants are held out together** (15,001 variants, 1,703 genes). About 0.03 of AUC comes from recognising genes rather than variant properties.
- **Interpretation:** discrimination is modest. This is a methods demonstration, not a clinical tool.

## Pipeline

- **Split:** stratified 70/30, made before any preprocessing.
- **Training-set-only fits, applied to the test set:**
  - dropping near-empty and high-cardinality columns;
  - random-forest imputation models;
  - IQR capping;
  - `log1p` on allele frequencies;
  - engineered features (`avg_AF`, `is_rare_variant`);
  - scaling;
  - the correlation filter.
- **Exploration:** chi-square association tests, PCA, and k-means (a cluster-membership feature, where test rows are assigned to the nearest training centroid).
- **Models:** logistic regression, CART and random forest, tuned by 5-fold CV with **upsampling inside each training fold** (`caret::trainControl(sampling = "up")`). All three share the same folds.
- **Ensembles:**
  - bagged trees, validated on a split carved from the original (not upsampled) data;
  - a majority vote;
  - a stacked logistic meta-learner trained on out-of-fold probabilities, with its threshold chosen on out-of-fold data.

## What the November 2026 revision changed

| Issue in the original | Fix | Effect |
|---|---|---|
| Cross-validation ran on data upsampled **before** the folds were cut, so copies of minority rows sat on both sides of a fold | Upsampling inside each training fold, on the original training set | RF CV AUC 0.762 → **0.713**, now consistent with its test AUC |
| Bagging validation split drawn from upsampled data | Split the original data first, then upsample the training part | Validation score no longer inflated |
| The "stacked ensemble" was a majority vote | Kept the vote under its correct name; added out-of-fold stacking | A true stacked model (test AUC 0.717) |
| The chi-square step selected `is.integer()` columns after scaling, so it tested 0 features | Select by number of distinct values; test on unduplicated data | 9 features tested, 8 associated (p < 0.05) |
| "Best model" chosen by test-set F1 | Choose by cross-validated AUC | Random forest |
| Accuracy and F1 only | Added AUC and AUPRC with baselines | Threshold-free comparison |
| Text described median/mode imputation and "new observations" | Text matches the code (RF imputation; examples are held-out test rows) | — |
| Random split only | Added a gene-grouped cross-validation check | Quantifies gene memorisation |
| Data fetched from a URL with no check | SHA-256 verification; `sessionInfo()` recorded | Reproducible inputs |

The original submission is kept unchanged in [`original/`](original/), and re-knitting it reproduces its reported numbers exactly (RF test AUC 0.7227, CV AUC 0.762).

## Reproduce

Requires R ≥ 4.5 and pandoc. Install the packages, then knit:

```r
install.packages(c("tidyverse", "caret", "corrplot", "rpart", "rpart.plot", "randomForest",
                   "e1071", "ipred", "digest", "rmarkdown"))
rmarkdown::render("analysis.Rmd", output_dir = "docs", output_file = "index.html")   # about 4 minutes
```

The data (`data/clinvar_conflicting.csv.zip`, 3.7 MB) is downloaded on first run from [`AKULA-PRASAD/da5030-data`](https://github.com/AKULA-PRASAD/da5030-data), and its SHA-256 is verified (`a3a988d5…6aaf`). Package versions used are printed by `sessionInfo()` at the end of the report.

## Limitations

- **Labels:** "Conflicting" partly reflects how many laboratories submitted a variant, and ClinVar interpretations change over time. Neither is a biological property of the variant.
- **Missing annotations:** PolyPhen and SIFT are more than 60% missing and were imputed, not measured.
- **Single split:** one random train/test split, so the gene-grouped check is the better guide to performance on unseen genes.
- **Not for clinical use.**

## Layout

| Path | Contents |
|---|---|
| `analysis.Rmd` | The full analysis (revised) |
| `docs/index.html` | Rendered report (GitHub Pages) |
| `original/` | Original Spring 2026 submission, unchanged |
| `data/` | Downloaded data (not committed) |
