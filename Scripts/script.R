library(ggplot2)
library(mzQuality2)


dataPaths <- c(
    #High = r"(D:\ConcentrationsCode\NMC-22-002_ADNI_hph.tsv)",
    Low = r"(D:\ConcentrationsCode\adni_low_ph_cleaned.tsv)"
)

concPath <- c(
    ACAL = r"(D:\ConcentrationsCode\concentrations_wei_ACAL.txt)",
    CAL = r"(D:\ConcentrationsCode\concentrations_wei.txt)"
)

outputs <- expand.grid(
    names(concPath), names(dataPaths)
)

to_remove <- c(
    "2202RAB_0001_A1",
    "2202RAB_0044_A1",
    "2202RAB_0079_A1",
    "2202RAB_0449_A1",
    "2202RAB_0620_A1",
    "2202RAB_0624_A1",
    "2202RAB_0644_A1",
    "2202RAB_0739_A1",
    "2202RAB_0767_A1",
    "2202RAB_0859_A1",
    "2202RAB_0006SQC_D1",
    "2202RAB_0005SQC_H1",
    "2202RAB_0006SQC_B1",
    "2202RAB_0013SQC_E1"
)



for (i in seq_len(nrow(outputs))) {

    calType <- outputs$Var1[i]
    ph <- outputs$Var2[i]

    cPath <- concPath[calType]
    path <- dataPaths[ph]


    exp <- buildExperiment(buildCombined(path))

    exp <- exp[,-which(colnames(exp) %in% to_remove)] %>%
        filterSST() %>%
        filterISTD() %>%
        doAnalysis(removeBadCompounds = FALSE) %>%
        addConcentrations(df = read.delim(cPath), filterComps = FALSE) %>%
        calculateConcentrations(type = names(cPath)[1]) %>%
        addBatchCorrectionAssay(assay = "Concentration")


    file <- glue::glue("ADNI_{type}_{ph}.xlsx", type = calType, ph = ph)

    exp <- exp[, exp$type %in% c("SAMPLE", "NIST", "SQC", "LQC")]
    writeConcentrations(file, exp)

#     file <- glue::glue("ADNI_{type}_{ph}_concentrations.xlsx", type = calType, ph = ph)
#     reportConcentrationsPerBatch(exp, file)
}



path <- r"(D:\ConcentrationsCode\adni_low_ph_cleaned.tsv)"
path <- r"(D:\ConcentrationsCode\NMC-22-002_ADNI_hph.tsv)"



calType <- "ACAL"
ph <- "High"

cPath <- concPath[calType]
path <- dataPaths[ph]


exp <- buildExperiment(buildCombined(path))

exp <- exp[,-which(colnames(exp) %in% to_remove)] %>%
    filterSST() %>%
    filterISTD() %>%
    filterOutliers() %>%
    doAnalysis(removeBadCompounds = TRUE) %>%
    addConcentrations(df = read.delim(cPath), filterComps = TRUE) %>%
    calculateConcentrations(assayX = "Ratio",type = calType, subtractCal = TRUE,
                            checkForOutliers = TRUE,
                            useWeights = FALSE, forceOrigin = TRUE
    )


assay(exp["BSH_LPI20.4", exp$batch == 2], "Concentration")



suggestedExp <- replaceInternalStandards(exp) %>%
    doAnalysis(doAll = TRUE) %>%
    addConcentrations(df = concs, filterComps = TRUE) %>%
    calculateConcentrations(assayX = "Ratio",type = "ACAL",
                            checkForOutliers = TRUE, useWeights = FALSE)



types <- c("SQC", "NIST", "ACAL")
compound <- 44

compoundPlotNew(exp, compound = 1, batches = 1, types = c("ACAL", "SAMPLE"), assay = "Ratio", withinTrend = TRUE, trendTypes = types) %>%
    facetPlot(ncol = 3)

compoundPlotNew(exp, compound = compound, types = types, assay = "Ratio_Corrected", withinTrend = TRUE, trendTypes = types)
compoundPlotNew(exp, compound = compound, types = types, assay = "WithinRatio", withinTrend = TRUE, trendTypes = types)
compoundPlotNew(exp, compound = compound, types = types, assay = "Concentration", withinTrend = TRUE, trendTypes = types)
compoundPlotNew(exp, compound = compound, types = types, assay = "Concentration_Corrected", withinTrend = TRUE, trendTypes = types)

compoundPlotNew(suggestedExp, compound = compound, types = types, assay = "Ratio", withinTrend = TRUE, trendTypes = types)
compoundPlotNew(suggestedExp, compound = compound, types = types, assay = "Ratio_Corrected", withinTrend = TRUE, trendTypes = types)
compoundPlotNew(suggestedExp, compound = compound, types = types, assay = "WithinRatio", withinTrend = TRUE, trendTypes = types)
compoundPlotNew(suggestedExp, compound = compound, types = types, assay = "Concentration", withinTrend = TRUE, trendTypes = types)
compoundPlotNew(suggestedExp, compound = compound, types = types, assay = "Concentration_Corrected", withinTrend = TRUE, trendTypes = types)



cowplot::plot_grid(
    compoundPlotNew(exp, compound = compound, types = types, assay = "WithinRatio", withinTrend = TRUE, trendTypes = types),
    compoundPlotNew(suggestedExp, compound = compound, types = types, assay = "WithinRatio", withinTrend = TRUE, trendTypes = types),
    compoundPlotNew(exp, compound = compound, types = types, assay = "Concentration_Corrected", withinTrend = TRUE, trendTypes = types),
    compoundPlotNew(suggestedExp, compound = compound, types = types, assay = "Concentration_Corrected", withinTrend = TRUE, trendTypes = types)
)




x <- biological[, biological$Batch == batch]
concentrations <- t(assay(x, "Concentration"))
r2 <- t(assay(x, "R2"))
calRange <- t(assay(x, "CalRange"))
calRange[is.na(calRange)] <- FALSE
highConfLinearRange <- colSums(calRange) / nrow(calRange) >= 0.7
highR2 <- r2[1, ] >= 0.95

codes <- do.call(c, lapply(1:length(highConfLinearRange), function(i){
    if (highConfLinearRange[i] & highR2[i]) {
        code <- 1
    } else if (highConfLinearRange[i]) {
        code <- 2
    } else if (highR2[i]) {
        code <- 3
    } else {
        code <- 4
    }
    code
}))

rsdqc <- rowData(x)$RSDQC.Corrected
background <- rowData(x)$backgroundSignal

header <- t(data.frame(
    `Sample Name` = rownames(x),
    RSDQC = as.double(rsdqc),
    Background = as.double(background),
    Concentration = codes,
    LOD = LoD
))



rbind(
    header,
    round(concentrations, digits)
)
})
