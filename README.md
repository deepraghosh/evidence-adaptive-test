# Evidence-adaptive global test (EA)

The **evidence-adaptive test (EA)** combines marginal and precision-adjusted evidence to test whether a multivariate mean vector is zero. Correlation adjustment can reveal a signal by removing predictable noise, but it can also attenuate changes shared by correlated coordinates. EA retains both representations and allows either, or their agreement at a coordinate, to support the global decision.

The canonical model is 

$$\boldsymbol{X} \sim N_d(\boldsymbol{\mu} ,\boldsymbol{\Sigma}),\qquad H_0:\boldsymbol \mu = \boldsymbol 0,\qquad H_1:\boldsymbol \mu \ne \boldsymbol 0,$$

with positive-definite covariance $\boldsymbol \Sigma$. The repository considers **known covariance** (Setting K) and **estimated covariance** (Setting W), in which an independent scatter matrix satisfies $\boldsymbol A\sim W_d(n,\boldsymbol \Sigma)$ with $n\ge d$.

EA combines the marginal and adjusted two-sided p-values $P_j,Q_j$ at each coordinate through a local Simes rule, then applies a Bonferroni correction across coordinates:

$$S_j=\min \lbrace 2\min(P_j,Q_j),\max(P_j,Q_j) \rbrace,\qquad P_{\mathrm{EA}}=\min\lbrace 1,d\min_j S_j\rbrace.$$

The analytic rule rejects when $P_{\mathrm{EA}}\le\alpha$. Under the stated Gaussian model and the appropriate known-covariance or independent-Wishart construction of the inputs, this rule controls the global Type I error in finite samples at level $\alpha$. The guarantee uses the model's dependence structure within each pair.

The real-data examples test $H_0: \boldsymbol \mu_1 = \boldsymbol \mu_2$ for two groups of protein profiles. Under independent Gaussian sampling with a common covariance, their scaled mean difference and pooled within-group scatter give the same canonical model.

| Folder | Contents |
|---|---|
| [`simulations/`](simulations/README.md) | Conditional Monte Carlo simulations, analytic EA comparisons, and instructions for reproducing tables and figures. |
| [`real data analysis/`](real%20data%20analysis/README.md) | Frozen TCGA measurements, p-value calculations, reference verification, and model/calibration notes. |

For each simulated dataset, the code generates a fresh set of null reference draws conditional on the observed augmented scatter matrix. All methods share these reference draws within the replicate, and estimated quantities are refitted on each draw. The reference sampler requires only the observed pair and the Wishart degrees of freedom. In Setting K, the test scores can use the known covariance. The real-data analyses use the analytic and asymptotic calibrations described in the real-data guide.

## Start here

RStudio is convenient for selecting individual run sections. Standard R is sufficient. Simulation acceleration through Rcpp and a C++ compiler is optional; the scripts never install packages automatically.

- **Observed data:** from the repository root, run `Rscript "real data analysis/RUN_data_analysis.R"`. This reads the frozen TCGA data, computes analytic $p$-values, and checks the 14 displayed values against the bundled references. 
- **Simulation tables and figures.** From the repository root, run `source("simulations/EA_run.R")` in R. The [simulation guide](simulations/README.md) and its [plain-text copy](simulations/README.txt) show how to load the saved outputs or rerun each design.
- **Analytic EA comparisons.** Table S7 compares analytic EA with conditionally calibrated tests on the same simulated datasets used for Tables 4 and 5. The script `simulations/EA_component_benchmarks.R` generates the component comparisons in Tables S8 and S9.

```text
evidence-adaptive-test/
├── README.md
├── simulations/
│   ├── README.md
│   ├── README.txt
│   ├── EA_run.R
│   ├── EA_methods.R
│   ├── EA_simulations.R
│   ├── EA_fast.cpp
│   ├── EA_component_benchmarks.R
│   └── results/
└── real data analysis/
    ├── README.md
    ├── RUN_data_analysis.R
    ├── EA_unknown_covariance_realdata.R
    ├── EA_manuscript_simulations.R
    ├── TCGA_Datahub_navigation_instructions.txt
    ├── data/
    └── [reference values, checksums, and source attribution]
```

## Reproducibility

The conditional runs use master seed `20261118` and `B = 99` reference draws per observed dataset. The table runs use `R = 1000` evaluation replicates per setting, and the figure runs use `R = 200`. The analytic component study uses seed `20261119` and `R = 1000`. All studies use level `.05`.

The `simulations/results/` folder contains the saved simulation results and the PDFs for Figures 1 and S1. The conditional-run results record the simulation settings and random-number streams.

The real-data folder includes frozen TCGA measurements with checksum verification and upstream attribution in [SOURCE_DATA.txt](real%20data%20analysis/SOURCE_DATA.txt). Data terms are recorded there separately from the analysis code.

