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
| [`simulations/`](simulations/README.md) | Reproducible simulation settings, table and figure instructions, analytic EA comparisons, and the paired two-stage bootstrap. |
| [`real data analysis/`](real%20data%20analysis/README.md) | Frozen TCGA measurements, p-value calculations, reference verification, and model/calibration notes. |

The main power comparisons use a common Monte Carlo null-rank calibration for all methods. This is distinct from EA's analytic finite-sample rule, which is evaluated separately in Table S7. The estimated-covariance simulation panel uses the true simulation covariance to generate its calibration bank. The real-data analyses report the analytic and asymptotic calibrations documented in their documentation.

## Start here

RStudio is convenient for selecting individual run sections. Standard R is sufficient. Simulation acceleration through Rcpp and a C++ compiler is optional; the scripts never install packages automatically.

- **Observed data:** from the repository root, run `Rscript "real data analysis/RUN_data_analysis.R"`. This reads the frozen TCGA data, computes analytic $p$-values, and checks the 14 displayed values against the bundled references. 
- **Simulation tables and figures:** open the [simulation reproduction guide](simulations/README.md). It maps every requested table and figure to its run sections, fixed settings, and saved outputs. Full runs use 19,999 calibration replicates and 10,000 evaluation replicates per covariance group.
- **Bootstrap and analytic EA:** `simulations/EA_additional_analyses.R` contains the two-stage paired bootstrap and Table S7. Its A0–A3 sections can be run separately; sourcing the whole file runs all of them.

```text
evidence-adaptive-test/
├── README.md
├── simulations/
│   ├── README.md
│   ├── README.txt
│   ├── EA_manuscript_simulations.R
│   ├── EA_explore_settings.R
│   ├── EA_additional_analyses.R
│   └── Run_confirmation.R
└── real data analysis/
    ├── README.md
    ├── RUN_data_analysis.R
    ├── EA_unknown_covariance_realdata.R
    ├── EA_manuscript_simulations.R
    ├── TCGA_DataHub_navigation_instructions.txt
    ├── data/
    └── [reference values, checksums, and source attribution]
```

The two copies of `EA_manuscript_simulations.R` are intentional: the data-analysis folder retains its frozen shared library, while the simulation folder includes the additional design functions. Keep each driver with the library beside it.

## Reproducibility

Simulation seeds are fixed integers. The five-variable and moving-average master seeds are 20261026 and 20261028; the two-stage bootstrap seed is 20261030. Full-precision outputs and RNG/configuration records are saved when the corresponding run sections execute. Large generated score banks, caches, and figure outputs are not included in this source repository.

The existing seeded bootstrap record has 90 positive simultaneous lower bounds. Its smallest unrounded bound is 0.0188522653536516, so all bounds exceed 0.0188 (standard four-decimal rounding of the minimum gives 0.0189). The code recalculates this result rather than assuming it. See the reproduction guide for the comparison family and exact commands.

The real-data folder includes frozen TCGA measurements with checksum verification and upstream attribution in [SOURCE_DATA.txt](real%20data%20analysis/SOURCE_DATA.txt). Data terms are recorded there separately from the analysis code.

