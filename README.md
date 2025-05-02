# mzQuality

MzQuality is a user-friendly R package for quality control of metabolomics 
studies. It features outlier detection, batch-correction using pooled study 
quality control samples (SQC), filters for removing unreliable compounds, 
various plots for inspecting and generating reports for further processing. 
See our [preprint](https://www.biorxiv.org/content/10.1101/2025.01.22.633547v1) 
for more information.

This R package forms the backbone of our shiny dashboard application 
[mzQualityDashboard](https://github.com/hankemeierlab/mzQualityDashboard), 
which is recommended for interactive use. To install and use the dashboard,
see the [mzQualityDashboard repository](https://github.com/hankemeierlab/mzQualityDashboard).

## Installing mzQuality

To install mzQuality, you can run the following script:

```{r install}
if (!"remotes" %in% installed.packages()) {
    install.packages("remotes", type = "binary)
}

if (!"mzQuality" %in% installed.packages()) {
    remotes::install_github("hankemeierlab/mzQuality", type = "binary")
}
```

## Quick Start guide

To get an idea of the capabilities of mzQuality, an example dataset containing x compounds and y samples has been added. These can be loaded by running the following code. The result is a SummarizedExperiment, which is the core object that mzQuality uses. See here for an overview on a SummarizedExperiment.

```{r}
exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
```

To use your own data, either a SummarizedExperiment with can be used, or a tab-delimited text file. A convenience function has been added to convert a tab-delimited file into a SummarizedExperiment:

```{r}
combined <- buildCombined(system.file("example.tsv", package = "mzQuality"))
exp <- buildExperiment(combined)
```

Once the SummarizedExperiment has been build, various filters and calculations can be applied. For convenience, a wrapper function has been designed that performs all analyses included in this package:

```{r}
exp <- doAnalysis(
    exp = exp
)
```

Selecting compounds & samples

Plots

Create an export


