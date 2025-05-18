---
editor_options: 
  markdown: 
    wrap: 72
---

# mzQuality

[![R-CMD-check](https://github.com/hankemeierlab/mzQuality/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hankemeierlab/mzQuality/actions/workflows/R-CMD-check.yaml)
[![BiocCheck](https://github.com/hankemeierlab/mzQuality/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/hankemeierlab/mzQuality/actions/workflows/bioc-check.yml)
[![Codecov test
coverage](https://codecov.io/gh/hankemeierlab/mzQuality/graph/badge.svg)](https://app.codecov.io/gh/hankemeierlab/mzQuality)

MzQuality is a user-friendly R package for quality control of
metabolomics studies. It features outlier detection, batch-correction
using pooled study quality control samples (SQC), filters for removing
unreliable compounds, various plots for inspecting and generating
reports for further processing. See our
[preprint](https://www.biorxiv.org/content/10.1101/2025.01.22.633547v1)
for more information.

This R package forms the backbone of our shiny dashboard application
[mzQualityDashboard](https://github.com/hankemeierlab/mzQualityDashboard),
which is recommended for interactive use. To install and use the
dashboard, see the [mzQualityDashboard
repository](https://github.com/hankemeierlab/mzQualityDashboard).

## Installing mzQuality

To install mzQuality, you can run the following script:

``` r
pkgs <- installed.packages()
if (!"remotes" %in% pkgs) {
    # Needed to retrieve development packages from github
    install.packages("remotes", type = "binary")
}

if (!"BiocManager" %in% pkgs) {
    # Needed for some mzQuality dependencies
    install.packages("BiocManager", type = "binary")
}

if (!"GenomeInfoDbData" %in% pkgs) {
    # Required by some dependencies, but not automatically installed 
    BiocManager::install("GenomeInfoDbData")
}

if (!"mzQuality" %in% pkgs) {
    # Install mzQuality
    remotes::install_github("hankemeierlab/mzQuality",type = "binary")
}
```

## Getting Started

To get an idea of the capabilities of mzQuality, an example dataset
containing x compounds and y samples has been added. These can be loaded
by running the following code. The result is a SummarizedExperiment,
which is the core object that mzQuality uses. See
[here](https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html)
for an overview on a SummarizedExperiment.

``` r
exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
```

To use your own data, either a SummarizedExperiment with can be used, or
a tab-delimited text file. See the vignette [Data
input](https://github.com/hankemeierlab/mzQuality/vignettes/Data_Input.html)
for an explanation for the format and mandatory columns to be present.
An example tab-delimited file can be seen
[here](https://github.com/hankemeierlab/mzQuality/inst/example.tsv)

Once your files are ready, you can use the `readData` function to read
in your data. It will check if all mandatory columns are present and if
the data is in the correct format. Finally, the function
`buildExperiment` will convert the data into a SummarizedExperiment
object. This is the core object that mzQuality uses to perform all
calculations and analyses and is used throughout the package.

``` r
# Path to the file
path <- system.file("example.tsv", package = "mzQuality")

# Read the file and check if all mandatory columns are present:
combined <- readData(path)

# Build the SummarizedExperiment object used throughout mzQuality:
exp <- buildExperiment(combined)
```

Once the SummarizedExperiment has been build, various filters and
calculations can be applied. For convenience, a wrapper function has
been designed that performs all analyses included in this package:

``` r
exp <- doAnalysis(exp = exp, doAll = TRUE)
```

The `doAnalysis` function will perform the following steps:

1.  Calculate the ratio between the compounds and assigned internal
    standards
2.  Perform batch correction using the pooled study quality control
    samples (SQC),
3.  Calculate the percentage of background signal compared to the study
    samples,
4.  Calculate the matrix effect, calculate the ratio of the QC sample,
5.  Calculate the presence of the compounds in the samples,
6.  Calculate the median area of the compounds in the samples,
7.  Suggest Internal Standards based on the calculated values.

If known concentrations for calibration lines have been supplied, the
`doAnalysis` function will also calculate the calibration line and the
corresponding R2 value. See the vignette [Data
input](https://github.com/hankemeierlab/mzQuality/vignettes/Data_Input.html)
for more information on how to supply concentrations.

All calculations will be added to the `assay`, `rowData` and `colData`
slots of the experiment, or overwrite the values if they are already
present. Calling `doAnalysis` again will overwrite previous values.

## Selecting compounds & samples

It is likely that mzQuality deems some compounds and/or samples
unreliable for reporting. It bases the decision for compounds on the
combination of RSDQC, the background signal percentage, and the presence
of the compounds in QC samples. The thresholds for these values can be
set in `doAnalysis`. For samples, this is based on the outcome of the
Rosner Test, which tests for statistical outliers in QC samples.

mzQuality adds a column called `use` in both the `rowData` and `colData`
slots of the SummarizedExperiment. These contain either a `TRUE` or
`FALSE` value, indicating if the compound or sample is reliable for
reporting based on the set thresholds. These values can be set manually
to `TRUE` or `FALSE` to override the automatic selection.

## Plots

Once a selection of the desired compounds and samples has been made,
various plots can be made to visualize the data. The following plot
functions are present:

``` r
# Boxplot of all values in aliquots/samples 
aliquotPlot(exp)

# Scatterplot of the compound / IS ratio for the first compound
compoundPlot(exp, assay = "ratio")

# Principal Component Analysis plot showing the aliquot distribution per
# type and batch for batch-corrected ratios
pcaPlot(exp, assay = "ratio_corrected")

# Violinplot showing the distribution for a given sample type
violinPlot(exp)
```

## Create an export

Once analysis and inspection has been completed, a report can be created
containing various plots and tables. The function `createReports` will
create a folder containing HTML files with plots, tab-delimited files
containing the **colData**, **rowData**, and the various **assays**, and
an Excel file that contains all the data in a single file.

Based on the set thresholds, mzQuality distinguishes between
`High Confidence`, `Caution`, `Low SNR` and `Non-Reportable` compounds
in the compoundSummary tab.

``` r
# Create an export in the current folder containing the reports
createReports(
    folder = getwd(),
    project = "mzQuality",
    exp = exp,
    backgroundPercent = 40, 
    cautionRSD = 15, 
    nonReportableRSD = 30
)
```

# Citing mzQuality

To cite mzQuality in publications, please use:

van der Peet M, Maas P, Wegrzyn A, Lamont L, Fleming R, Harms A,
Hankemeier T, Kindt A (2025). “mzQuality: A tool for quality monitoring
and reporting of targeted mass spectrometry measurements.” *bioRxiv*.
<doi:10.1101/2025.01.22.633547>
<https://doi.org/10.1101/2025.01.22.633547>,
<https://www.biorxiv.org/content/early/2025/01/24/2025.01.22.633547>.
