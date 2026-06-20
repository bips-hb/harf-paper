## High-dimensional adversarial random forests

This repository contains the code for the paper "High-dimensional adversarial random forests". Before running an `R` script, please make sure to set the working directory to the root of this repository.
The code is organized as follows:
- `R-code`: contains all `R` scripts, chronologically numerated.

### Downstream clustering analysis 
- `01-prepare-hember-data`: contains the code to prepare the Hember data for the analysis.
- `02-preproceed-hember-data.R`: contains the code to preprocess the Hember data for the analysis.
- `03-hember-problem.R`, `04-hember-algorithm.R`, `05-hember-submission.R`, `06-hember-results`: contain the batchtools code to conduct the downstream clustering analysis.
- `07-hember-pairwise-cor.R`: contains the code to conduct the empirical pairwise Pearson correlation between features.

### Downstream prediction analysis
- `08-prepare-currated-tcga-data`: contains the code to download and preprocess the TCGA data for the analysis.
- `09-tcga-problem.R`, `10-tcga-algorithm.R`, `11-tcga-submission.R`, `12-tcga-results.R`: contain the code to conduct the downstream prediction analysis.

### Case Study: Alzheimer's Disease Brain Metabolomics
- `15-AD`: contains the code to prepare and analyze the Alzheimer's disease brain metabolomics data for the analysis.
- `15-AD/01-ad-problem.R`, `15-AD/02-ad-algorithm.R`, `15-AD/01-ad-submission.R`, `15-AD/03-ad-results.R`: contain the code to conduct the downstream prediction analysis for the Alzheimer's disease brain metabolomics data.