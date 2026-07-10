# GFBST_Pragmatico

Implementation of the **Generalized Full Bayesian Significance Test (GFBST)** for pragmatic hypotheses.

This repository contains the R code developed as part of my PhD research on Bayesian hypothesis testing with pragmatic hypotheses. The goal is to provide implementations and examples used in the numerical studies presented in the thesis.

## Abstract Tese
Agnostic tests are an approach used in statistical hypothesis testing, with three decisions: accept the hypothesis of interest, reject the hypothesis of interest or remain agnostic. The Generalized Full Bayesian Significance Test (GFBST) is considered a good test, as it maintains desirable logical properties and statistical characteristics. However, under certain conditions, GFBST never accepts a precise null hypothesis, which limits its use. To overcome this limitation, pragmatic hypotheses are considered, which are non-precise hypotheses constructed in a way that is sufficiently close to the original precise hypotheses and useful to overcome this problem. Pragmatic hypothesis tests allow the acceptance of the null hypothesis in appropriate cases, making them a viable alternative for dealing with precise hypotheses in scientific research. In this work, we develop solutions for pragmatic hypothesis tests in problems involving categorical data in homogeneity testing, and we present simulations and real-data applications that illustrate the methodology we propose.

https://teses.usp.br/teses/disponiveis/45/45133/tde-26062026-181614/

## Repository structure

```
├── Códigos - Simulações/
├── Códigos - Aplicações Reais/
├── 00_TH_kx2_funcoes_dist.R
├── 00_TH_kx2_funcoes_hpd.R
├── 01_TH_kx2_gerar_dados.R
├── 02_TH_kx2_resultados_dist.R
├── 03_TH_kx2_resultados_hpd.R
├── 04_TH_kx2_resultados_estimaprobs.R
├── 05.1_TH_kx2_resultados_funcoes_aplicacao.R
└── GFBST_Pragmatico.Rproj
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

## Reference

Massoni, G. (2026). *Generalized Full Bayesian Significance Test for Pragmatic Hypotheses*. PhD Thesis, University of São Paulo.
