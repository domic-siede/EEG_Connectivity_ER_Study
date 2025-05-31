# EEG_Connectivity_ER_Study
MATLAB scripts for computing wPLI-based EEG connectivity during emotion regulation

# EEG Connectivity Analysis for Emotion Regulation

This repository contains MATLAB code used to compute and analyze oscillatory brain connectivity across different emotion regulation conditions, focusing on theta, alpha, beta, and delta bands. The scripts are part of a study exploring how attachment orientations modulate brain connectivity during cognitive reappraisal and expressive suppression.

## 🧠 Overview

The pipeline includes three scripts:

### 1. `connectivity_per_condition_neutral_BANDS_1.m`

- Computes condition-specific connectivity matrices using FieldTrip's `wpli_debiased` method for each EEG frequency band.
- Compares each experimental condition (Reappraise, Suppress, Negative) against the Neutral baseline.
- Applies a data-driven threshold (95th percentile) per band to retain salient connectivity increases.
- Saves group-level connectivity and masked matrices for follow-up analysis.

### 2. `extract_subject_connectivity_diff_vs_neutral_2.m`

- Uses the group-derived masks to extract individual-level connectivity differences versus Neutral for each band and condition.
- Outputs a CSV file with subject-wise connectivity differences, which serves as input for the statistical modeling.

### 3. `run_LME_models_connectivity_3.m`

- Loads the individual-level connectivity CSVs for each band.
- Fits Linear Mixed-Effects Models (LMMs) using condition, attachment anxiety (ANX), and attachment avoidance (AVD) as predictors.
- Outputs model results for interpretation.

## ⚙️ Dependencies

- MATLAB R2022b (or later)
- [EEGLAB](https://sccn.ucsd.edu/eeglab/index.php) (tested with version 2023.0)
- [FieldTrip toolbox](https://www.fieldtriptoolbox.org/) (tested with release 20230118)

## 📁 Directory Structure

├── connectivity_per_condition_neutral_BANDS_1.m
├── extract_subject_connectivity_diff_vs_neutral_2.m
├── run_LME_models_connectivity_3.m
├── avgMatricesByBand.mat
├── electrodePairsByBand.mat
└── subject_connectivity_differences.csv

## 📬 Contact

For questions or support, please contact:  
**Marcos Domic-Siede**  
[mdomic@ucn.cl](mailto:mdomic@ucn.cl)
