# Nutritional Impact of India's First COVID-19 Lockdown

This repository contains the R code used for data management, statistical analysis, and figure generation for the manuscript:

Marri AR, Dauphinais MR, et al. Nutritional Impact of India's First COVID-19 Lockdown: Was It Equitable?

## Data Source

Data come from the National Family Health Survey (NFHS-5), part of the Demographic and Health Surveys (DHS) Program. NFHS-5 data are publicly available upon registration at https://dhsprogram.com. Due to DHS data use agreements, raw data are not included in this repository.

## Analytical Overview

Analyses include:

- Survey-weighted (unadjusted) prevalence estimation
- Covariate-adjusted prevalence estimation by marginal standardization (adjusting for age, caste, and religion)
- Survey-weighted logistic regression (adjusted odds ratios)
- Concentration index estimation and concentration curve visualization
- Forest plots of adjusted odds ratios
- Sensitivity analyses restricted to the 13 states and union territories with observations in both pre- and post-lockdown periods

All analyses were conducted in R (version 4.5.1).

## Repository Structure

Primary analyses:
- men_nfhs5_aug2026.R – Data cleaning, analysis, and figure generation for men aged 15–54 years
- women_nfhs5_aug2026.R – Data cleaning, analysis, and figure generation for women aged 15–49 years
- child_nfhs5_aug2026.R – Data cleaning, analysis, and figure generation for children under five years of age

Sensitivity analyses (restricted to the 13 overlapping states and union territories):
- 99a_Sensitivity_analyses_men.R
- 99b_Sensitivity_analyses_women.R
- 99c_Sensitivity_analyses_children.R

## Required R Packages

Key packages used include: survey, dplyr, tidyr, ggplot2, convey, Hmisc, readr. Users should install all required packages before running the scripts.

## Reproducibility

To reproduce the analyses:

1. Register for access to NFHS-5 data through the DHS Program.
2. Download the required NFHS-5 datasets.
3. Update file paths within the scripts to match local directories.
4. Run the primary analysis scripts (men, women, children), then the sensitivity scripts.

Because DHS data are not redistributed, complete replication requires independent access to the original NFHS-5 datasets.

## License

This code is released under the MIT License (see LICENSE file).

## Citation

If you use code from this repository, please cite:
Marri, A. R., Dauphinais, M. R., Martinez, L., Karoly, M., McQuaid, F., & Sinha, P. Nutritional impact of India's first COVID-19 lockdown: Was it equitable?
