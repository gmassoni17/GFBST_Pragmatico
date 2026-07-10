# GFBST_Pragmatico

Implementation of the **Generalized Full Bayesian Significance Test (GFBST)** for pragmatic hypotheses.

This repository contains the R code developed as part of my PhD research on Bayesian hypothesis testing with pragmatic hypotheses. The goal is to provide implementations and examples used in the numerical studies presented in the thesis.

## Abstract Tese
Agnostic tests are an approach used in statistical hypothesis testing, with three decisions: accept the hypothesis of interest, reject the hypothesis of interest or remain agnostic. The Generalized Full Bayesian Significance Test (GFBST) is considered a good test, as it maintains desirable logical properties and statistical characteristics. However, under certain conditions, GFBST never accepts a precise null hypothesis, which limits its use. To overcome this limitation, pragmatic hypotheses are considered, which are non-precise hypotheses constructed in a way that is sufficiently close to the original precise hypotheses and useful to overcome this problem. Pragmatic hypothesis tests allow the acceptance of the null hypothesis in appropriate cases, making them a viable alternative for dealing with precise hypotheses in scientific research. In this work, we develop solutions for pragmatic hypothesis tests in problems involving categorical data in homogeneity testing, and we present simulations and real-data applications that illustrate the methodology we propose.

https://teses.usp.br/teses/disponiveis/45/45133/tde-26062026-181614/

## Repository structure

```
.
├── R/            # R functions
├── data/         # Example datasets (if applicable)
├── scripts/      # Simulation and analysis scripts
├── figures/      # Generated figures
└── README.md
```

## Requirements

* R (version 4.2 or later recommended)
* Required R packages listed in the scripts

## Usage

Clone the repository:

```bash
git clone https://github.com/gmassoni17/GFBST_Pragmatico.git
```

Open the project in RStudio and run the desired scripts.

## Citation

If you use this code in academic work, please cite the corresponding thesis or publication.

## Author

Gabriela Massoni

PhD  in Statistics.

This project is licensed under the MIT License.
