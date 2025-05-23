# mzQuality 0.99.2
- Moved data files from inst/* to inst/extdata/* .
- Added description of data files in inst/script/README.md.
- Changed extensions of .txt files to .tsv to reflect the tab-separated format.
- Moved the files in inst/rmd/* to inst/templates/* and renamed them to
  reflect their purpose.
- Added `mzQuality` vignette, showcasing a typical mzQuality workflow.
- Removed `doAll` and `aliquots` arguments from `doAnalysis` as they were
  no longer used

# mzQuality 0.99.1

- Renamed `Getting_Started` vignette to `Data_Input` and added information on 
how to read and load data in mzQuality.
- Simplified adding concentrations by supplying the `concentration` column
- Changed assay name that calculates batch-corrected concentrations with
the prefix of the sample type of the calibration line. 
- Changed example data to match example experiment used throughout testing.

# mzQuality 0.99.0

- Initial version
- Reformatted code for Bioconductor submission.
- Added README, LICENSE, and NEWS files.
