library(dplyr)
library(ggplot2)
library(mzQuality2)

path <- r"(H:\ConcentrationsCode\adni_low_ph_cleaned.tsv)"


concsACAL <- read.delim(r"(H:\ConcentrationsCode\concentrations_wei_ACAL.txt)")
comb <- buildCombined(path)
 # comb <- comb %>%
 #     select(-c(area_is, rt_is, compound_is)) %>%
 #     as.data.frame()


# write.table(comb, file = r"(D:\ConcentrationsCode\adni_low_ph_noIS.tsv)", sep = "\t", row.names = F)

 exp <- buildExperiment(comb) %>%
        filterSST() %>%
        filterISTD() %>%
        identifyOutliers() %>%
        filterOutliers() %>%
        doAnalysis(removeBadCompounds = FALSE, removeOutliers = F) %>%
        addConcentrations(concsACAL, filterComps = FALSE) %>%
        calculateConcentrations(type = "ACAL") %>%
        addBatchCorrectionAssay(assay = "concentration", removeOutliers = F) %>%
        carryOverEffect(assay = "concentration_corrected")

rowData(exp)$carryOver

rowData(exp)$backgroundSignal
rowData(exp["BSL_5_iPF2a_VI", ])$carryOver


compoundPlotNew(exp, compound = "BSL_TUDCA", batches = 1, assay = "area") # good


compoundPlotNew(exp, compound = "BSL_5_iPF2a_VI", batches = 5, assay = "area") # too low

compoundPlotNew(exp, compound = "BSL_4_HDoHE", batches = 5, assay = "area") # too low


calculateEffect(exp = exp["BSL_5_iPF2a_VI", exp$batch == 4], type = "BLANK")
calculateEffect(exp = exp["BSL_10_NO2_LA", exp$batch == 1], type = "BLANK")


assay(exp["BSL_16_HDoHE", c("2202RAB_0001CAL8_A1", "2202RAB_0001BLANK_A2")], "concentration_corrected")


sum(exp$type == "SQC")
nrow(exp)


rowData(exp_outliers)

rowData(exp)

t.test(
    as.vector(assay(exp_outliers["BSL_thromboxane_B2", exp$type == "SQC"], "ratio_corrected")),
    as.vector(assay(exp["BSL_thromboxane_B2", exp$type == "SQC"], "ratio_corrected")),
    paired = T
)






cowplot::plot_grid(ncol = 1,
    compoundPlotNew(exp, compound = "BSL_thromboxane_B2", batches = 5, assay = "concentration_corrected"),
    compoundPlotNew(exp_outliers, compound = "BSL_thromboxane_B2", batches = 5, assay = "concentration_corrected")
)

compoundPlotNew(exp_outliers, compound = "BSL_thromboxane_B2", batches = 5, assay = "area")


boxplot(
    data.frame(
        RSDQC_corrected = rowData(exp)$rsdqcCorrected,
        RSDQC_outlier_corrected = rowData(exp_outliers)$rsdqcCorrected
    )

)





m <- t(mzQuality2:::calculateCorrectedRSDQCs2(exp))
rownames(m)[apply(m, 1, which.min)]

apply(m, 1, min)


concentrationPlotNew(exp[which(rowData(exp)$hasKnownConcentrations), ], batch = c(1,2)) %>% facetPlot()


violinPlotNew(exp)


unlist(pcaPlotNew(exp)$text[1])


metadata <- as.data.frame(colData(exp))

res <- lapply(1:nrow(metadata), function(i){
    names <- colnames(metadata)
    values <- as.vector(metadata[i, ])
    paste0(sprintf("%s: %s", names, values), collapse = "\n")
})
names(res) <- rownames(metadata)


res[[1]]

paste(c(1,2,3), c(4,5,6))



pcaPlotNew(exp, assay = "concentration")

aliquotPlot(exp, batch = 1)

assay(exp, "Concentration")[1:5, 1:5]


mzQuality2:::downloadZip("test", exp, "test.zip")

exp
rowData(exp)



compoundReports(exp, folder = getwd())
summaryReport(exp, folder = getwd(), assays = c("Ratio_Corrected", "Concentration_Corrected"),
              plots = c("aliquot", "pca", "violin", "heatmap"))

library(foreach)

mzQuality2:::create_report(
    path = getwd(),
    parameters = list(exp = exp, plots = "", assays = c("Ratio", "Ratio_Corrected", "Concentration_Corrected")),
    template = "mzquality-compounds-se.Rmd",
    output = "test.html"
)


assayNames(exp)

assay(exp, "WithinRatio")[1:5, 1:5]
assay(exp, "Ratio_Corrected")[1:5, 1:5]


exp2 <- replaceInternalStandards(exp) %>%
    doAnalysis(doAll = TRUE) %>%
    calculateConcentrations(type = "ACAL") %>%
    addBatchCorrectionAssay(assay = "Concentration")

summaryReport(exp2, getwd(), assays = c("Ratio_Corrected", "Concentration_Corrected"))

library(ggplot2)
concentrationPlotNew(exp, compound = 3)



library(ggplot2)
mzQualityDashboard::toPlotly(compoundPlotNew(exp, compound = rownames(exp)[1], batches = c("1", "2")))


exp <- exp[intersect(rownames(exp), rownames(concsCAL)), ]


concs <- list(
    CAL = concsCAL,
    ACAL = concsACAL
)
assay <- 'WithinRatio'
for (calType in names(concs)) {
    exp <- exp %>%
        addConcentrations(concs[[calType]])

    cols <- exp$type == "CAL" & exp$calno == min(exp$calno, na.rm = T)
    assay(exp, assay)[cols] <- 0

    m <- do.call(cbind, lapply(unique(exp$batch), function(i){
        batch <- exp[, exp$type == calType & exp$batch == i]

        yVals <- assay(batch, "Concentration")
        yVals <- yVals[complete.cases(yVals), ]

        xVals <- assay(batch[rownames(yVals), ], assay)


        slope <- rowWiseSlope(xVals, yVals)
        int <- rowWiseIntercept(xVals, y = yVals, slope = slope)
        subset <- exp[rownames(yVals), exp$batch == i]
        xVals <- assay(subset, assay)

        xVals * slope + int
    }))
    #m[m < 0] <- 0
    assay(exp, glue::glue("predicted_{type}", type = calType)) <- m
}


rowData(exp)


# Good!
compound <- 3
cowplot::plot_grid(
    ncol = 2,
    compoundPlotNew(exp[compound, ], batches = c(1,2, 3, 4), assay = "Ratio_Corrected", shortName = T),
    compoundPlotNew(exp[compound, ], batches = c(1,2, 3, 4), assay = "WithinRatio", shortName = T),
    compoundPlotNew(exp[compound, ], batches = c(1,2, 3, 4), assay = "predicted_ACAL", shortName = T),
    compoundPlotNew(exp[compound, ], batches = c(1,2, 3, 4), assay = "predicted_CAL", shortName = T)
)

## BAD due to batch differences in (A)CALs
compound <- 18
cowplot::plot_grid(
    ncol = 2,
    compoundPlotNew(exp[compound, ], batches = c(1,2, 3, 4), assay = "Ratio_Corrected", shortName = T),
    compoundPlotNew(exp[compound, ], batches = c(1,2, 3, 4), assay = "WithinRatio", shortName = T),
    compoundPlotNew(exp[compound, ], batches = c(1,2, 3, 4), assay = "predicted_ACAL", shortName = T),
    compoundPlotNew(exp[compound, ], batches = c(1,2, 3, 4), assay = "predicted_CAL", shortName = T)
)



compound <- 2
cowplot::plot_grid(
    ncol = 2,
    compoundPlotNew(exp[compound, ], batches = c(1,2, 3, 4), assay = "Ratio_Corrected", shortName = T),
    compoundPlotNew(exp[compound, ], batches = c(1,2, 3, 4), assay = "WithinRatio", shortName = T),
    compoundPlotNew(exp[compound, ], batches = c(1,2, 3, 4), assay = "predicted_ACAL", shortName = T),
    compoundPlotNew(exp[compound, ], batches = c(1,2, 3, 4), assay = "predicted_CAL", shortName = T)
)




withinBatchCorrectionVectorize <- function(exp, assay, samples, qcType){
    x <- exp[, samples]
    qcs <- x[, x$type == qcType]
    vals <- assay(qcs, assay)

    orderMatrix <- matrix(rep(qcs$order, nrow(vals)), ncol = ncol(qcs), byrow = T)
    slope <- rowWiseSlope(orderMatrix, vals)
    int <- rowWiseIntercept(orderMatrix, y = vals, slope = slope)
    orderMatrix <- matrix(rep(1:ncol(x), nrow(vals)), ncol = ncol(x), byrow = T)
    return(orderMatrix * slope + int)
}

rowData(exp)

colData(exp) %>%
    as.data.frame() %>%
    filter(.data$type == "CAL") %>%
    filter(.data$batch == 3)

colData(exp[, exp$batch == 3])

x <- replaceInternalStandards(exp)

df <- data.frame(
    rsdqc = rowData(exp)$rsdqc,
    rsdqcBetween = rowData(exp)$rsdqcCorrected,
    rsdqcSuggested = rowData(exp)$suggestedRSDQC,
    rsdqcWithin = rowData(exp)$rsdqcWithin,
    rsdqcSuggestedWithin = rowData(x)$rsdqcWithin
)




# Rowwise Calibration Line Modelling --------------------------------------


x <- exp[, exp$batch == 1]

cals <- x[, x$type == "ACAL"]
vals <- assay(cals, "Ratio_Corrected")
concMatrix <- assay(cals, "Ratio")

which(x$calno == min(x$calno, na.rm = T))

slope <- rowWiseSlope(vals, concMatrix)
int <- rowWiseIntercept(vals, y = concMatrix, slope = slope)
assay(x, "Ratio_Corrected") * slope + int



df <- stack(df)
ggplot(df %>% filter(values <= 30), aes(y = values, group = ind, fill = ind)) + geom_boxplot()

batches <- length(unique(exp$batch))
facetCols <- round(sqrt(batches))


cowplot::plot_grid(
    compoundPlotNew(exp, compound = 2, batches = 1:3, assay = "Ratio", shortName = FALSE,
                    trendTypes = c("CAL", "ACAL", "SQC"), withinTrend = TRUE),
    compoundPlotNew(exp, compound = 2, batches = 1:3, assay = "Ratio_Corrected", shortName = FALSE,
                    trendTypes = c("CAL", "ACAL", "SQC"), withinTrend = TRUE),
    compoundPlotNew(exp, compound = 2, batches = 1:3, assay = "WithinRatio", shortName = FALSE,
                    trendTypes = c("CAL", "ACAL", "SQC"), withinTrend = TRUE),
    ncol = 3
)



cowplot::plot_grid(
    compoundPlotNew(exp, compound = 2, batches = 1:3, assay = "Ratio", shortName = FALSE,
                    trendTypes = c("SQC"), withinTrend = F),
    compoundPlotNew(exp, compound = 2, batches = 1:3, assay = "Ratio_Corrected", shortName = FALSE,
                    trendTypes = c("SQC"), withinTrend = F),
    compoundPlotNew(exp, compound = 2, batches = 1:3, assay = "WithinRatio", shortName = FALSE,
                    trendTypes = c("SQC"), withinTrend = F),
    ncol = 3
)




compoundPlotNew(exp, compound = 3, types = "SQC", assay = "Ratio_Corrected", shortName = FALSE, withinTrend = TRUE)

type <- "SQC"
assay <- "Ratio"
cowplot::plot_grid(
    compoundPlotNew(exp, compound = 5, types = type, assay = "Ratio_Corrected", shortName = FALSE,
                    trendTypes = c(type), withinTrend = TRUE),

    compoundPlotNew(exp, compound = 5, types = type, assay = "WithinRatio", shortName = FALSE,
                    trendTypes = c(type), withinTrend = TRUE),
    ncol = 1
)


compoundPlotNew(exp, compound = 5, types = type, batches = 1, assay = assay, shortName = FALSE,
                trendTypes = c(type), withinTrend = TRUE)


rsdqc(withinBatchCorrection(exp), assay = "WithinRatio")

calculateCorrectedRSDQCs2(exp)

x <- addBatchCorrectionAssay(exp = exp, assay = "WithinRatio")

sum(rsdqc(exp, assay = "Ratio_Corrected") >= rsdqc(x, "WithinRatio_Corrected"))

rowData(exp)$withinRSDQC <- rsdqc(x, "WithinRatio_Corrected")
rowData(exp)

which(assay(x[, x$type == "SQC"], "WithinRatio_Corrected") < 0, arr.ind = T)

assay(x[, x$type == "SQC"][41,], "Ratio_Corrected")

vals <- as.vector(assay(x[, x$type == "SQC" & x$batch == 1][41, ], "Ratio"))

i <- 41
j <- 1
x <- exp[, exp$batch == j]

vals <- as.vector(assay(x[i, ], "Ratio"))

df <- data.frame(
    order = 1:ncol(x),
    value = vals,
    type = x$type
)

df <- df %>%
    filter(.data$type == "SQC")

df
#df$order <- 1:nrow(df)

m <- lm(value ~ order, df)
m
int <- as.vector(coef(m)[1])
slope <- as.vector(coef(m)[2])
df$value / (1:ncol(x) * slope)[df$order]

df$value

return((1:ncol(x) * slope) + int)
plot(df$value)

assay(x[, x$type == "SQC"][41, 9], "WithinCorrection")


rowData(exp)

rowData(exp[17, ])

compound_is <- unique(rowData(exp)$compound_is)
res <- lapply(compound_is, function(IS){
    comps <- rownames(exp)[which(rowData(exp)$compound_is == IS)]

    res <- lapply(comps, function(compound){
        cowplot::plot_grid(
            ncol = 1,
            compoundPlotNew(exp, compound = compound, types = "SQC", assay = "area", trendTypes = "SQC", shortName = T, withinTrend = T),
            compoundPlotNew(exp, compound = compound, types = "SQC", assay = "area_is", trendTypes = "SQC", shortName = T, withinTrend = T),
            compoundPlotNew(exp, compound = compound, types = "SQC", assay = "Ratio", trendTypes = "SQC", shortName = T, withinTrend = T),
            compoundPlotNew(exp, compound = compound, types = "SQC", assay = "Ratio_Corrected", trendTypes = "SQC",shortName = T, withinTrend = T),
            #compoundPlotNew(exp, compound = compound, types = "SQC", assay = "WithinCorrection", trendTypes = "SQC",shortName = T, withinTrend = T),
            compoundPlotNew(exp, compound = compound, types = "SQC", assay = "WithinRatio", trendTypes = "SQC",shortName = T, withinTrend = T)#,
            #compoundPlotNew(x, compound = compound, types = "SQC", assay = "WithinRatio_Corrected",trendTypes = "SQC", shortName = T, withinTrend = T)
        )
    })
    names(res) <- comps
    res
})
names(res) <- compound_is

# To Show in powerpoint
# Response of area does not match its internal standard..
# whats going on???
res$BSL_d4_6_keto_PGF1a_ISTD$BSL_6_keto_PGF1a
res$BSL_d5_FA22.6_w3_ISTD$BSL_FA22.6_w3
res$BSL_d4_CA_ISTD$BSL_CA
res$BSL_d4_12_13_DiHOME_ISTD$BSL_12_13_DiHOME
res$BSL_d8_12_HETE_ISTD$BSL_12_HETE




type <- "SAMPLE"
compound <- 17
cowplot::plot_grid(
    ncol = 1,
    compoundPlotNew(exp, compound = compound, types = c("SQC", type), assay = "area", trendTypes = c("SQC", type), shortName = T, withinTrend = T),
    compoundPlotNew(exp, compound = compound, types = c("SQC", type), assay = "area_is", trendTypes = c("SQC", type), shortName = T, withinTrend = T),
    compoundPlotNew(exp, compound = compound, types = c("SQC", type), assay = "Ratio", trendTypes = c("SQC", type), shortName = T, withinTrend = T),
    compoundPlotNew(exp, compound = compound, types = c("SQC", type), assay = "Ratio_Corrected", trendTypes = c("SQC", type), shortName = T, withinTrend = T),
    compoundPlotNew(exp, compound = compound, types = c("SQC", type), assay = "WithinRatio", trendTypes = c("SQC", type), shortName = T, withinTrend = T)
)

boxplot(rowData(exp)$rsdqcWithin)

nrow(exp)
compoundPlotNew(exp, compound = 8, trendTypes = c("CAL", "ACAL", "SQC"),
                assay = "Ratio_Corrected", shortName = FALSE, withinTrend = T) %>%
    plotly::ggplotly()

rownames(exp)[8]

aliquotPlotNew(exp, "Ratio_Corrected", batches = c(1:4), types = "ACAL") %>%
    facetPlot(ncol = 4) %>%
    plotly::ggplotly()



violinPlotNew(exp, assay = "area", batches = c(1, 26), types = "SQC") %>%
    facetPlot(ncol = 1, shareX = FALSE)



pcaPlotNew(exp, assay = "WithinRatio", types = c("SQC", "LQC", "SAMPLE"),
           sampleAsBatch = T, addConfidenceInterval = T) #%>%
    facetPlot(ncol = 5, shareY = T, shareX = T)


    plotly::ggplotly()








