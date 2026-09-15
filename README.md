# TMTL

Code for **Joint-Sparse Transfer Learning for High-Dimensional Multi-Output Regression**.

## Repository structure

- `Codes/`: estimation functions and C++ helpers.
- `Simulations/`: simulation experiments with $p = 200$ and $K$ = 10, 20, and 30.
- `Realdata/Site1_p2500/2500HVG/`: Site 1 real-data analysis starting from 2,500 highly variable genes before additional filtering. For the final $p$ and $K$, refer to the main text.

## Estimation functions

- `Codes/functions_joint.R`: joint-sparse multitask estimation functions used for **TMTL(Fused)** and **MTL(Target)**. The multitask estimator is also used for **MTL(Full)** with pooled target and source data.
- `Codes/functions_stl.R`: taskwise single-response transfer estimation functions used for **TSTL(Fused)**.

The corresponding method-specific scripts implement the target-based
debiasing steps for **TMTL(Debiased)** and **TSTL(Debiased)**.


## Requirements

The estimation code uses the following R packages:

```r
install.packages(c(
  "glmnet", "Matrix", "dplyr",
  "Rcpp", "RcppArmadillo", "RcppEigen"
))
```

The C++ helpers require a working C++ compiler compatible with R.
The optional preprocessing script additionally requires `zellkonverter`,
`SingleCellExperiment`, and `Seurat`. 

Before running the scripts, adjust `source()`, `sourceCpp()`, and data paths
to match the location of the downloaded repository. Shared functions and
C++ helpers are provided in `Codes/`. Relative paths are resolved from
the current R working directory, which can be checked with `getwd()`.

## Methods

The directory names correspond to the following methods in the manuscript:

| Method in the manuscript | Simulation directory | Real-data directory |
|---|---|---|
| **TMTL(Fused)** | `MTL_transfer` | `MTL_transfer_cluster` |
| **TMTL(Debiased)** | `MTL_transfer_debiased` | `MTL_transfer_debiased` |
| **MTL(Target)** | `MTL_notransfer` | `MTL_notransfer` |
| **MTL(Full)** | `MTL_notransfer_merged` | `MTL_notransfer_merged` |
| **TSTL(Fused)** | `STL_transfer_cluster_sparse` | `STL_transfer_cluster` |
| **TSTL(Debiased)** | `STL_transfer_debiased` | `STL_transfer_debiased` |

- **TMTL(Fused):** joint-sparse multitask transfer estimator before target-based debiasing.
- **TMTL(Debiased):** joint-sparse multitask transfer estimator with target-based debiasing.
- **MTL(Target):** multitask estimator fitted using target data only.
- **MTL(Full):** multitask estimator fitted by pooling the target and source data under a common coefficient matrix.
- **TSTL(Fused):** taskwise single-response transfer estimator before debiasing.
- **TSTL(Debiased):** taskwise single-response transfer estimator with target-based debiasing.

The real-data directory also contains `nullmodel/`, which fits the
intercept-only reference model using target-training means for each protein.

## Simulation experiments

The simulation directories are:

- `Simulations/p200K10/`
- `Simulations/p200K20/`
- `Simulations/p200K30/`

Each contains four source-shift settings:

| Setting directory | Scenario |
|---|---|
| `threesource_noskew` | Balanced shifts |
| `threesource_aligned_shift0.25` | Aligned shifts with alpha = 0.25 |
| `threesource_aligned_shift0.5` | Aligned shifts with alpha = 0.5 |
| `threesource_aligned_shift1` | Aligned shifts with alpha = 1 |

Each setting uses **100 replications**.

### Execution order

1. Choose a `p200K` directory and a setting.
2. Load the included `threesource.Rdata`, which contains the simulation
   settings and coefficient matrices. Alternatively, regenerate it by
   running `datacreation.R` from the setting directory, then load the
   generated file.
3. Load the required functions from `Codes/`.
4. Set the working directory to the relevant method directory and run `sim.R`.
5. Repeat for all six methods and all four settings.

Run `MTL_transfer` before `MTL_transfer_debiased`, and
`STL_transfer_cluster_sparse` before `STL_transfer_debiased`.
The debiased methods load the corresponding fused estimates.

When starting a fresh R session, load the settings and required functions
again before running a method.

Each method generates observations using `set.seed(j)` for replication
`j` and saves its results to `summary.Rdata`. These generated result files
are omitted from the repository. Rerunning a method overwrites its result file.

### Summaries and figures

After completing all methods and settings, run
`mseplotting_threesource_frobnorm_negtransfer_predictionerror.R`
from the selected `p200K` directory.

This script summarizes coefficient error and prediction RMSE, computes
negative-transfer frequencies relative to **MTL(Target)**, and generates:

- `logfrobnorm_threesource.pdf`
- `predrmse_threesource.pdf`






## Real-data analysis

Considering the repository size, we provide **Site 1 as a representative example**
of the four-site analysis. The same analysis workflow applies to the
other sites, using their corresponding preprocessed data, domain indices,
and source-domain configurations.

The uploaded example is located in:

`Realdata/Site1_p2500/2500HVG/`

It contains 12 target cell types and uses **20 random-split replications**. 
Number of source domains $L$, $p$ (dimensionality), and $K$ (the number of tasks) may differ by Sites, and the details are in the main text.

### Files

| File | Purpose |
|---|---|
| `domain_index.csv` | Mapping between domain indices and cell types |
| `clr_data_donor_nofiltered_2500hvg*.Rdata` | Preprocessed domain-specific input data |
| `simulationsetting.R` | Tuning parameters, split proportions, and replication count |
| `dataprocess_site1.R` | Preprocessing code for Site 1 |
| `target1.R` through `target12.R` | Analysis scripts for the respective target domains |
| `rmse_domainwise_comparedtonull.r` | Domain-wise relative-MSE summaries |

The preprocessed inputs are included. To rerun preprocessing, obtain the
[GSE194122 dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122)
and update the `h5ad_file` path in `dataprocess_site1.R` to point to the
BMMC H5AD file.

### Execution order

1. Load the required functions from `Codes/` and the settings in
   `simulationsetting.R`.
2. Set the working directory to a method directory and run the appropriate
   `target*.R` script.
3. Repeat for all 12 targets and all six methods.
4. Run the corresponding target scripts in `nullmodel/`.
5. Run `rmse_domainwise_comparedtonull.r` from the `2500HVG` directory.

For each target, run `MTL_transfer_cluster` before
`MTL_transfer_debiased`, and `STL_transfer_cluster` before
`STL_transfer_debiased`.

When starting a fresh R session, load the required functions and settings
again, including for the null-model scripts.

Each target script saves its results as
`summary_logtransform_target*.Rdata` in the method directory.
Rerunning a target script overwrites its corresponding result file.

### Evaluation metric

The real-data summary script reports **relative MSE**, despite the
`rmse` prefix in its filename:

```math
Q_{m,r}^{(\ell)} = \frac{\mathrm{MSE}_{m,r}^{(\ell)}}{\mathrm{MSE}_{\mathrm{null},r}^{(\ell)}}.
```

For each method and target domain, the script reports the mean and sample
standard deviation of this ratio across 20 random-split replications.
Smaller values indicate better predictive performance.
