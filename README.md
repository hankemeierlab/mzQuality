# mzQuality

MzQuality is a user-friendly R package for quality control of metabolomics studies. It features outlier detection, batch-correction using pooled study quality control samples (SQC), filters for removing unreliable compounds, various plots for inspecting and generating reports for further processing. See our preprint for more information.

This R package forms the backbone of our shiny dashboard application `mzQualityDashboard`, which is recommended for interactive use. To install and use the dashboard, see the `mzQualityDashboard` repository.

## Installing mzQuality

To install mzQuality, you can run the following script:

```{r install}
if (!"remotes" %in% installed.packages()) {
    install.packages("remotes")
}

remotes::install_github(".../mzQuality")
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

Once the SummarizedExperiment has been build, various filters and calculations can be applied. 
