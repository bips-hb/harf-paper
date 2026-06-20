## High-Dimensional Adversarial Random Forests

This repository contains the code accompanying the paper *"High-Dimensional Adversarial Random Forests"*. Before running any `R` script, set the working directory to the root of this repository.

The code is organized in the `R-code` directory, where scripts are numbered according to the analysis workflow.

### Downstream Clustering Analysis

- `01-prepare-hember-data.R`: Prepares the Hember dataset.
- `02-preprocess-hember-data.R`: Preprocesses the Hember dataset.
- `03-hember-problem.R`, `04-hember-algorithm.R`, `05-hember-submission.R`, `06-hember-results.R`: Implement and evaluate the downstream clustering experiments using `batchtools`.
- `07-hember-pairwise-cor.R`: Computes pairwise Pearson correlations between features.

### Downstream Prediction Analysis

- `08-prepare-curated-tcga-data.R`: Downloads and preprocesses the TCGA datasets.
- `09-tcga-problem.R`, `10-tcga-algorithm.R`, `11-tcga-submission.R`, `12-tcga-results.R`: Implement and evaluate the downstream prediction experiments.

### Case Study: Alzheimer's Disease Brain Metabolomics

- `15-AD/`: Contains all scripts for data preparation and analysis of the Alzheimer's disease brain metabolomics case study.
- `15-AD/01-ad-problem.R`, `15-AD/02-ad-algorithm.R`, `15-AD/03-ad-submission.R`, `15-AD/04-ad-results.R`: Implement and evaluate the downstream prediction analysis for the Alzheimer's disease dataset.
