# Evaluating User Satisfaction with NTHU's Academic Information System (AIS)

**Author:** Dale John Baltazar  
**Program:** International MBA (IMBA), National Tsing Hua University  
**Method:** Partial Least Squares Structural Equation Modeling (PLS-SEM)

## Overview

This repository contains the data and analysis code for an IMBA thesis examining the determinants of user satisfaction with NTHU's Academic Information System. The study uses PLS-SEM to test five hypothesized paths from information quality (IQ), system quality (SQ), service quality (SERVQ), mobile user experience (MUX), and complementary assets (CA) to overall user satisfaction (US).

## Repository Structure

```
├── data/
│   └──  Raw survey data are not publicly shared to protect participant 
confidentiality. De-identified summary statistics are available upon 
reasonable request to the corresponding author.
├── scripts/
│   └── analysis.R            # Complete
├── figures/                   # Generated figures
├── NTHU_AIS_Analysis.Rproj   # RStudio project
├── .gitignore
└── README.md
```

## Data

The dataset contains 81 survey responses (77 usable after cleaning) with 33 columns:

| Columns | Content |
|---------|---------|
| 1–6 | Administrative (ID, timestamps, email, name, international status) |
| 7–10 | Demographics (age group, gender, level of study, AIS use frequency) |
| 11–14 | Information Quality (IQ1–IQ4) |
| 15–18 | System Quality (SQ1–SQ4) |
| 19–21 | Service Quality (SERVQ1–SERVQ3) |
| 22–25 | Mobile User Experience (MUX1–MUX4) |
| 26–29 | Complementary Assets (CA1–CA4) |
| 30–33 | User Satisfaction (US1–US4) |

All survey items use a 5-point Likert scale (1 = Strongly Disagree to 5 = Strongly Agree).

## How to Reproduce

1. Clone this repository.
2. Open `NTHU_AIS_Analysis.Rproj` in RStudio.
3. Install required packages:
   ```r
   install.packages(c("psych", "dplyr", "tidyr", "ggplot2", "seminr", "gridExtra"))
   ```
4. Source the analysis script:
   ```r
   source("scripts/analysis.R")
   ```

The script will:
- Load and recode the survey data
- Print respondent demographics
- Compute item-level and construct-level descriptive statistics
- Assess reliability (Cronbach's α, ρ_A, CR) and convergent validity (AVE)
- Evaluate discriminant validity (Fornell–Larcker, HTMT)
- Estimate the PLS-SEM model and run bootstrap inference (5,000 resamples, seed = 2024)
- Report hypothesis test results, effect sizes (f²), and VIF
- Generate and save Figures 2–8 to `figures/`

## R Packages

| Package | Purpose |
|---------|---------|
| `psych` | Descriptive statistics, reliability analysis |
| `dplyr` | Data manipulation |
| `tidyr` | Data reshaping |
| `ggplot2` | Visualization |
| `seminr` | PLS-SEM estimation and bootstrapping |
| `gridExtra` | Multi-panel figure layout |

## License

This project is shared for academic purposes. Please cite appropriately if using any part of this work.
