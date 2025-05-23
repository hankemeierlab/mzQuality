# External Data Descriptions

This document describes the contents of the `extdata` folder.
The files `data.RDS`, `sciex.tsv` and `example.tsv` contained in the folder 
are part of the same data set, while `concentrations.tsv` contains a table
of known concentrations of spiked compounds in the calibration lines. 

The data set comprises of 6 random yet consecutively measured batches 
containing 419 EDTA plasma study samples originally part of the 
[CoSTREAM consortium](https://costream.eu/index.html). 
Also included were different types of QC samples, blank samples, 
and calibration curves. The samples were prepared and measured according 
to the method described by [Yang. et al (2024)](https://www.sciencedirect.com/science/article/pii/S0003267024001491).

## sciex.tsv
The data has been integrated using SCIEX OS (version 2.1.6.59781) and 
exported using the tab-delimited text format. It contains a subset of 
the formentioned study by only containing two study QC samples. 
This file is used to test if the export format of Sciex OS is parsed 
correctly and thus can be used as an alternative input format for mzQuality. 

## example.tsv
This file contains a small subset of the data set and is meant for 
instructing the user on how to use the recommended format as input 
for mzQuality. The format is explained in the `Data-Input`
vignette. 

## data.RDS
This file is a saved _SummarizedExperiment_ object containing data 
used as both example and test data. This object contains the full data set
as described in the first section and was also integrated and exported 
using Sciex OS to a text file. This file was read using the 
`readData` function and converted to a SummarizedExperiment object 
using the `buildExperiment` function. Lastly, this file was saved 
using `saveRDS` and included in this package as a showcase of the 
capabilities of mzQuality for a multi-batch data set.

## concentrations.tsv
The concentrations.tsv file is a tab-delimited file containing a 
template for an alternative method for supplying known concentrations
of spiked compounds. This can be useful when the supplied data file
is in sciex-os text format. The format is explained in the `Data-Input`
vignette. 
