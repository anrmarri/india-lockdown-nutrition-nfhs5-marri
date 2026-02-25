# Nutritional Impact of India’s First COVID-19 Lockdown

This repository contains R code used for data management, statistical analysis, and figure generation for the manuscript:

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
- Concentration curve visualization

All analyses were conducted in R (version 4.5.1).

## Repository Structure

- `scripts/` – R scripts for cleaning, analysis, and figure generation for men, women, and children under five

## Reproducibility

Users must download NFHS-5 datasets directly from DHS and update file paths accordingly before running scripts.
