# Tau-Astrocyte Analysis

This repository contains custom MATLAB and R scripts used for image analysis, quantification, statistical analysis, and visualization in a study investigating tau uptake, localization, and processing in human iPSC-derived astrocytes.

## Repository Contents

### MATLAB

The MATLAB scripts are used for image-based analysis of tau and astrocyte/lysosomal localization.

#### `Tau_LysoTracker_Colocalization_2D.m`
Calculates Pearson correlation coefficients between tau and LysoTracker fluorescence signals in 2D images. The script processes paired fluorescence channels and summarizes correlation values across experimental groups.

#### `Tau_LysoTracker_Colocalization_ZStack.m`
Calculates Pearson correlation coefficients between tau and LysoTracker signals using confocal z-stack images.

#### `Tau_Inside_Astrocytes_ZStack_Quantification.m`
Identifies astrocyte regions and tau-positive pixels in z-stack images and calculates the number of tau-positive pixels located within astrocyte regions, the total number of astrocyte pixels, and the percentage of tau-positive pixels within astrocytes.

### R

The R scripts are used for statistical analysis and visualization of live-cell imaging data.

#### `ART_ANOVA_2Condition_LiveCell.R`
Performs aligned rank transform (ART) analysis for experiments containing two experimental conditions measured across multiple time points. The script includes a mixed-effects model with well-level random effects, post-hoc comparisons, multiple-comparison corrections, summary statistics, and visualization.

#### `ART_ANOVA_4_Conditions.R`
Performs aligned rank transform (ART) analysis for experiments containing four experimental conditions measured across multiple time points. The script includes a mixed-effects model with well-level random effects, same-timepoint post-hoc comparisons, multiple-comparison corrections, summary statistics, and visualization.

## Requirements

### MATLAB

The MATLAB scripts require MATLAB with the Image Processing Toolbox.

### R

The R scripts require R and the following packages:

- `readxl`
- `dplyr`
- `ARTool`
- `ggplot2`
- `writexl`
- `stringr`
- `wesanderson`

## Data

Raw experimental imaging data and experimental spreadsheets are not included in this repository.

The MATLAB scripts expect image data to be organized according to the folder structure specified within each script. The R scripts require appropriately formatted Excel input files containing the variables referenced in each analysis script.

## Reproducibility

The scripts were developed and modified for the analyses performed in this study. Parameters such as image-processing thresholds, background correction settings, and statistical analysis settings are defined within the respective scripts.

Users should review the parameter sections of each script before applying the code to new datasets.

## Statistical Analysis

Statistical analyses of live-cell imaging data were performed using aligned rank transform (ART) models with experimental condition and time as fixed effects and well-level random effects where applicable. Post-hoc comparisons and multiple-comparison corrections are specified within the individual R scripts.

## Citation

The associated manuscript is currently in preparation. Citation information will be added upon publication.

Repository DOI:

The code associated with this study is archived on Zenodo:
XXXXX
