# Nutritional Impact of India’s First COVID-19 Lockdown

This repository contains R codes used for data management, statistical analysis, and figure generation for the manuscript:

Marri AR, Dauphinais MR, et al. Nutritional Impact of India’s First COVID-19 Lockdown: Was It Equitable?

## Data Source

Data come from the National Family Health Survey (NFHS-5), part of the Demographic and Health Surveys (DHS) Program.

NFHS-5 data are publicly available upon registration at: https://dhsprogram.com

Due to data use agreements, raw data are not included in this repository.

## Analytical Overview

Analyses include:

- Survey-weighted prevalence estimation
- Survey-weighted logistic regression
- Concentration index estimation
- Forest plots for weighted odds ratios
- Concentration curve visualization
- Sensitivity analyses restricted to overlapping states and union territories

All analyses were conducted in R (version 4.5.1).

## Repository Structure
- men_nfhs5_may2026.R – Data cleaning, analysis, and figure generation for men aged 15–54 years
- women_nfhs5_may2026.R – Data cleaning, analysis, and figure generation for women aged 15–49 years
- children_nfhs5_may2026.R – Data cleaning, analysis, and figure generation for children under five years of age
- nfhs5_phase2_sensitivity_combined.R – Sensitivity analyses restricted to the 13 states and union territories with observations available in both pre- and post-lockdown periods

## Required R Packages
Key packages used include:
- survey
- dplyr
- tidyr
- ggplot2
- convey
- Hmisc
- readr

Users should install all required packages before running the scripts.

## Reproducibility
To reproduce the analyses:
1. Register for access to NFHS-5 data through the DHS Program.
2. Download the required NFHS-5 datasets.
3. Update file paths within the scripts to match local directories.
4. Run the analysis scripts for men, women, and children.

Because DHS data are not redistributed, complete replication requires independent access to the original NFHS-5 datasets.

## Citation
If you use code from this repository, please cite:
Marri, A. R., Dauphinais, M. R., Martinez, L., Karoly, M., McQuaid, F., & Sinha, P. (2025, October 8). Nutritional impact of India’s first COVID-19 lockdown: Was it equitable? [Preprint]. SSRN. https://doi.org/10.2139/ssrn.5569486
