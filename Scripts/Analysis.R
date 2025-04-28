# # Perform analysis
path <- r"(D:\2023\ADNI\adni_low_ph_cleaned_3.tsv)"
#
# # Build combined should be simplified
df <- buildCombined(path)

conc_path <- r"(D:\2023\ADNI\Concentration list 50uL_Sample_mzQuality_Signaling 20221020_JG_ACAL.txt)"
concs <- read.delim(conc_path)

#combined <- arrow::read_delim_arrow(path, delim = "\t")
exp <- buildExperiment(buildCombined(path)) %>%
    addConcentrations(df = concs) %>%
    doAnalysis(doAll = TRUE)


exp[rowData(exp)$Use, ]
assay(exp)

df2 <- mzQuality:::buildCombined(path)


head(df2)
head(df)
exp2 <- mzQuality::buildExperiment(df2)
rowData(exp2)
rowData(exp)

exp3 <- mzQuality::buildExperiment(as.data.frame(df))
rowData(exp3)


profvis::profvis(doAnalysis(buildExperiment(buildCombined(path)), doAll = TRUE))

#
#
exp <- buildExperiment(df)
exp
library(mzQuality2)
rowData(exp)
#
# exp <- doAnalysis(exp, aliquots = which(exp$Use))
#
# rowData(exp)
#
# bench::mark(
#   doAnalysis(exp, doAll = TRUE)
# )
#
# rowData(exp)
#
# bench::mark(
#   doAnalysis(exp, compounds = which(rowData(exp)$Use)) #aliquots = which(colData(exp)$Use)
# )
#
# doAnalysis(exp, compounds = which(rowData(exp)$Use), aliquots = which(colData(exp)$Use))
#
# exp
# doAnalysis(mzQualityExp, compounds = 1:10)
#
# mzQualityExp

library(mzQuality2)
comb <- buildCombined(system.file(package = "mzQuality2", "example.tsv"))

exp <- buildExperiment(comb)


doAnalysis(exp[1, ])
