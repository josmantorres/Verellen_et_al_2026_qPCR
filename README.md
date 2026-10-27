# Statistical Analyses of PER and qPCR Data

This repository contains the R scripts and input datasets used for the statistical analyses of the proboscis extension response (PER) and quantitative PCR (qPCR) experiments reported in Verellen et al. (2026).

The analyses cover sucrose-evoked PER in different behavioral categories of honey bees and gene expression in the antennal lobe, brain, and antennae.

## Repository Contents

```text
.
├── README.md
├── data/
│   ├── PER in situ.xlsx
│   └── qPCR data.xlsx
├── scripts/
│   ├── PER_insitu_script_Verellenetal2026.R
│   └── qPCR_script_Verellenetal2026.R
└── LICENSE
```

Adjust the file paths in the scripts if the input workbooks are stored in a different location. The file and directory names shown above should match those used in the repository.

## Requirements

The analyses were performed in R and use the following packages.

**PER analysis**
- `readxl`
- `pastecs`
- `ggplot2`
- `lme4`
- `ResourceSelection`
- `sjPlot`
- `car`
- `emmeans`
- `dplyr`

**qPCR analysis**
- `readxl`
- `ggplot2`
- `ggpubr`
- `glmmTMB`
- `car`
- `DHARMa`
- `emmeans`
- `dplyr`
- `statmod`
- `tweedie`

Install any missing packages before running the scripts.

## 1. Proboscis Extension Response (PER)

### Data

The PER dataset is provided in `PER in situ.xlsx`. The script uses the following variables:

- `Season`: experimental season.
- `date`: experimental date.
- `hive`: colony identifier (renamed from `#hive` in the input file).
- `category`: behavioral category.
- `Sucrose concentration (%w/v)`: sucrose concentration used in the assay (renamed to `sucrose` in the script).
- `PER`: binary response, coded as **1 when the proboscis extension response occurs and 0 when it does not occur**.

Behavioral categories are:
- `NF`: nectar followers.
- `PF`: pollen followers.
- `C`: non-follower controls.

The sucrose concentrations included in the script are 0, 0.1, 0.3, 1, 3, 10, and 30 (% w/v).

### Statistical analysis

The script summarizes the data and fits candidate binomial models with a logit link. The candidate models include behavioral category and sucrose concentration as additive predictors, a model including their interaction, a model treating sucrose concentration as numeric, and a generalized linear mixed model (GLMM) with colony as a random intercept.

The candidate models are compared using Akaike's Information Criterion (AIC). The script selects the mixed model:

```r
PER ~ category + sucrose + (1 | hive)
```

Model evaluation includes a Hosmer–Lemeshow goodness-of-fit test. The script reports model coefficients and odds ratios, performs tests of model effects, and uses estimated marginal means (`emmeans`) for pairwise comparisons among behavioral categories. It also calculates estimated PER probabilities across sucrose concentrations and generates a publication-style plot with standard-error bars.

## 2. Quantitative PCR (qPCR)

### Data and behavioral categories

The qPCR dataset is provided in `qPCR data.xlsx`, with separate worksheets for:

- `Lobulo antenal`: antennal lobe.
- `Cerebro`: brain.
- `Antena`: antennae.

The script uses the following variables:

- `categoria`: behavioral category (`SN`, `SP`, or `NS` in the input data).
- `gen`: target gene.
- `muestra`: sample identifier.
- `ct`: expression value used in the statistical models.

The behavioral-category labels are interpreted in the figures as:
- `SN` = nectar followers (`NF`).
- `SP` = pollen followers (`PF`).
- `NS` = non-follower controls (`C`).

### Raw Ct values and relative expression

**Important:** The workbook `qPCR data.xlsx` contains raw Ct values extracted directly from the thermocycler. However, the variable named `ct` in the statistical analysis script represents relative gene expression (*fold change*) calculated using the \(2^{-\Delta\Delta Ct}\) method, as described in the manuscript's Materials and Methods.

The qPCR analysis script does **not** calculate relative expression from raw Ct values. Before running the statistical models and generating the figures, the raw Ct measurements must therefore be transformed into relative expression values using the reference genes and normalization procedure described in the manuscript.

The manuscript reports the use of *GAPDH* and *RPS18* as reference genes. Relative expression quantification follows the \(2^{-\Delta\Delta Ct}\) approach adapted for multiple reference genes and gene-specific amplification efficiencies. Inter-plate normalization was performed using the global average of calibrator samples across plates. Consult the manuscript's Materials and Methods for the experimental and normalization details required to reproduce these values.

The manuscript also describes quality-control procedures for technical triplicates, including exclusion of individual Ct values that differed by more than 0.6 cycles within a triplicate and values identified as outliers due to technical failures in reference-gene amplification.

In addition, the qPCR R script excludes six rows from the antenna dataset (rows 24–26 and 78–80 after import), described in the script as samples with abnormally high Ct values for the Rp18s internal control in the control–NPF group.

### Statistical analysis

The qPCR data are analyzed separately for the antennal lobe, brain, and antennae. Models include behavioral category, gene, and their interaction as fixed effects, with sample identity (`muestra`) as a random effect.

The script evaluates candidate distributions and model diagnostics. The selected model specifications are:

- **Antennal lobe:** Gamma GLMM with a log link.
- **Brain:** Gamma GLMM with a log link.
- **Antennae:** Tweedie GLMM with a log link and a dispersion model that includes behavioral category, gene, and their interaction.

Type III analyses of deviance are used to assess model effects. Pairwise comparisons between behavioral categories are performed within each gene using estimated marginal means (`emmeans`). For the brain and antennae analyses, pairwise p-values are adjusted using the Benjamini–Hochberg procedure.

The scripts generate plots of expression values by behavioral category and gene. The publication-style plots show boxplots and individual observations, with significance brackets for selected pairwise comparisons.

## Reproducibility

1. Clone or download this repository.
2. Install the required R packages listed above.
3. Place the Excel workbooks in the locations expected by the scripts, or update the paths in the scripts.
4. For qPCR, transform the raw Ct values into relative expression values using the reference genes and normalization procedure described in the manuscript.
5. Run the R scripts to fit the statistical models and generate the figures.

The qPCR script assumes that the relative expression values have already been calculated and supplied in the format expected by the models.

## Data and Code Availability

The R scripts and input datasets used for the PER and qPCR analyses are available in this repository.

**Manuscript:** Verellen et al. (2026). *[Differential sucrose sensitivity and neuropeptide expression in peripheral and central sensory tissues among dance followers of the honey bee, Apis
mellifera]*.
