# Real data analysis

This folder reproduces the TCGA protein-profile comparisons from the original patient-level observations. The data originate from [NCI’s TCGA program](https://www.cancer.gov/ccg/research/genome-sequencing/tcga). Run `RUN_data_analysis.R` to calculate the $p$-values for the colorectal and stomach-cancer examples and compare them with the included reference results.

The analysis tests equality of two population mean vectors,

$$H_0:\boldsymbol\mu_1=\boldsymbol\mu_2,$$

against a difference in at least one coordinate. These are global tests: a small $p$-value indicates evidence of a difference somewhere in the protein profiles; it does not identify the individual proteins responsible.

## Requirements and files

Only standard R is required. The four frozen input files are included, so the analysis runs offline and does not download data or install packages.

| File | Purpose |
|---|---|
| [`RUN_data_analysis.R`](RUN_data_analysis.R) | Driver for loading, analysis, combined results, and reference verification. |
| [`EA_unknown_covariance_realdata.R`](EA_unknown_covariance_realdata.R) | Data readers, group definitions, two-sample summaries, and p-value calculations. |
| [`EA_manuscript_simulations.R`](EA_manuscript_simulations.R) | Shared score and precision-estimation helpers, loaded as a function library. |
| [`data/`](data/) | The four frozen RPPA and clinical annotation files. |
| [`EXPECTED_VALUES.csv`](EXPECTED_VALUES.csv) | Previously calculated p-values and logarithms, used only for verification. |
| [`INPUT_MD5.csv`](INPUT_MD5.csv) | Checksums of the four input files, checked automatically by the driver. |
| [`SOURCE_DATA.txt`](SOURCE_DATA.txt) | Data attribution, upstream links, and redistribution terms. |
| [`REFERENCE_SESSION.txt`](REFERENCE_SESSION.txt) | Reference R version and platform. |
| [`MANIFEST_SHA256.txt`](MANIFEST_SHA256.txt) | Integrity checksums for the distributed files. |

Keep this folder's structure intact. In particular, the data must remain under `data/`.

## Run in RStudio

1. Open `RUN_data_analysis.R` in this folder.
2. Choose **Session → Set Working Directory → To Source File Location**.
3. Run **1. SETUP**, which loads both supporting scripts and verifies the input checksums.
4. Run **2. COLORECTAL**, followed by **3. STAD**.
5. Run **4. COMBINED TABLE**, then **5. VERIFY**.

Alternatively, click **Source** to run the complete driver. Do not run the other two scripts separately; the driver sources them in the required order. Loading the shared simulation engine does not run simulations.

From the repository root, the command-line equivalent is:

```sh
Rscript "real data analysis/RUN_data_analysis.R"
```

From another R working directory, source the driver by its path:

```r
source("/path/to/evidence-adaptive-test/real data analysis/RUN_data_analysis.R")
```

The driver locates the other files relative to its own saved location.

## Read and optionally save the results

After a successful run, inspect:

```r
View(comparison_table)
View(verification)
View(colorectal_result$table)
View(stad_result$table)
```

`comparison_table` contains the seven displayed methods for both datasets, their p-values, and the calibration descriptions. `verification` compares the calculated base-10 logarithms with `EXPECTED_VALUES.csv`, using a tolerance of `1e-6` to allow small numerical differences between platforms. The reference values do not enter the p-value calculation or determine significance; matching them verifies reproducibility of the implemented calculations.

The driver prints results and retains them in R memory. It does not automatically save result files. To save the combined table explicitly:

```r
write.csv(comparison_table, "TCGA_p_values.csv", row.names = FALSE)
```

More detailed components remain available as `colorectal_result$qcq_details`, `colorectal_result$zws_details`, and `colorectal_result$pooled_summaries`, with corresponding components in `stad_result`.

Some tail probabilities are smaller than ordinary R numbers can represent. In those cases, `p_value` is `NA`, `p_value_display` preserves a scientific-notation string, and `log10_p_value` remains finite. Use the display or logarithm columns rather than interpreting `NA` as a failed logarithmic calculation or replacing it with zero.

## Data and preprocessing

| Dataset | Groups | Group sizes | Proteins, $d$ | Pooled-scatter degrees of freedom, $n$ |
|---|---|---:|---:|---:|
| Colorectal | `COAD_CIN` versus `COAD_GS` | 175 and 45 | 189 | 218 |
| Stomach adenocarcinoma (STAD) | `N0` versus stages beginning `N1`, `N2`, or `N3` | 108 and 230 | 193 | 336 |

The colorectal groups contain colon adenocarcinomas. Rows of the analysis matrices are patients, and columns are proteins. The loader reads RPPA values as supplied without an additional log transformation, retains proteins with finite measurements across the whole cohort before selecting groups, matches samples to patient annotations, and checks for repeated patients.

The frozen measurements and annotations are from **The Cancer Genome Atlas (TCGA), PanCancer Atlas studies**, distributed through [cBioPortal DataHub](https://github.com/cBioPortal/datahub): [colorectal](https://github.com/cBioPortal/datahub/tree/master/public/coadread_tcga_pan_can_atlas_2018) and [stomach adenocarcinoma](https://github.com/cBioPortal/datahub/tree/master/public/stad_tcga_pan_can_atlas_2018). The bundled bytes, identified by `INPUT_MD5.csv`, define the reproduction inputs; a particular upstream commit has not been established for these copies. See [`SOURCE_DATA.txt`](SOURCE_DATA.txt) for the upstream database terms and study publications. The analysis code does not claim ownership of the TCGA measurements.

## What is calculated

Under independent Gaussian samples with a common positive-definite covariance matrix, let $m_1,m_2$ be the group sizes and define

$$X=\sqrt{\frac{m_1m_2}{m_1+m_2}}(\overline Y_1-\overline Y_2),\qquad n=m_1+m_2-2.$$

Let $A$ be the sum of the two within-group scatter matrices. Then $X$ and $A$ are independent, $A\sim W_d(n,\Sigma)$, and the null gives $X\sim N_d(0,\Sigma)$. This connects the two-sample analysis to the canonical Gaussian–Wishart model used by EA; $A/n$ is the pooled covariance estimate. The inverse-based construction requires $n\ge d$.

EA, Z-Max, Feng-COM, MinP, and Quadratic use this scaled mean difference and within-group scatter. QCQ and ZWS operate on the original observation matrices `Y1` and `Y2`. There are **no reconstructed observations, Gaussian draws, permutations, bootstrap samples, or seeds** in this analysis. The original coordinate scales are retained, except that the Z-Max nodewise Lasso uses pooled coordinate standardization to put its penalty on a comparable scale across proteins.

## Methods and interpretation

| Method | Calibration used in this data analysis |
|---|---|
| EA | Analytic finite-sample Gaussian–Wishart calibration under the common-covariance model. |
| MinP | Bonferroni combination of marginal Student p-values under that model. |
| Quadratic | Exact Gaussian-model F calibration. |
| Z-Max | Asymptotic Gumbel approximation, with conditions on dependence and precision estimation. |
| Feng-COM | Retained summary-based statistic with asymptotic calibration. |
| QCQ | Documented two-sample adaptation: a Chen–Qin quadratic base with the retained maximum/higher-criticism enhancement and provisional normal calibration. |
| ZWS | Full order-2 and order-6 two-sample mean U-statistics with a centered analytic variance adaptation and normal/minimum-p calibration. |

Z-TopK is unavailable because the retained implementation requires simulation-based calibration. IHC is unavailable because the retained implementation requires known covariance. Both remain `NA` with explanatory labels in the full result objects.

The QCQ and ZWS two-sample versions differ from the canonical simulation scores. In particular, the complete primary-paper QCQ two-sample specification has not been verified; its implemented adaptation is documented in the source comments. The extreme colorectal QCQ/ZWS probabilities are formal tail extrapolations whose finite-sample tail accuracy has not been established.

EA's finite-sample Type I error statement applies under the Gaussian/common-covariance assumptions, rather than to arbitrary observational distributions. STAD was identified in an exploratory screen; these p-values are not adjusted for that selection. Matching the supplied reference results establishes computational reproducibility and does not validate those model assumptions or asymptotic approximations.
