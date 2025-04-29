#' @title Studentized residuals
#' @description Studentized residuals represent t-distributed residuals,
#' meaning that each residual is scored according to the t-distribution.
#' This allows for determining outliers in regression lines
#' @details This function effectively calculates the studentized residuals without
#' recalculating the model. The basic idea is to delete the observations one at
#' a time, each time refitting the regression model on the remaining n–1
#' observations. Then, we compare the observed response values to their
#' fitted values based on the models with the ith observation deleted. This
#' produces deleted residuals. Standardizing the deleted residuals produces
#' studentized residuals.
#' @returns Matrix of t-distributed residuals
#' @param residuals Matrix of residuals where each row is one compound and
#' the each column is an observation from the training set. This is usually
#' a CAL or ACAL line
#' @param n Number of observations. Should be equal to the amount of non-NA
#' values in the training set.
#' @param k Number of predictors. Here it is always 1
#' @seealso https://online.stat.psu.edu/stat462/node/247/ for the formula used here
#' @export
studentResiduals <- function(residuals, n, k = 1) {
    rsd <- sqrt(rowSums(residuals^2, na.rm = TRUE) / (n - 2))
    standardResiduals <- residuals / rsd
    studentResiduals <- standardResiduals * ((n - k - 2) / (n - k - 1 - standardResiduals^2))**0.5
    studentResiduals[!is.finite(studentResiduals)] <- NA
    return(studentResiduals)
}

#' @title Calculate concentratios for a full experiment
#' @returns SummarizedExperiment with predicted concentrations in the assay
#' given in `saveAssay`. Defaults to `"concentration"`.
#' @param byBatch Should the concentration be calculated per batch? Defaults to
#' TRUE.
#' @param saveAssay Name of the assay to save the results. Defaults to
#' `"Concentration"`
#' @export
calculateConcentrations <- function(exp, byBatch = TRUE, saveAssay = "concentration", ...) {
    batches <- unique(exp$batch)
    if (byBatch) {
        x <- lapply(batches, function(i) {
            batchConcentrations(exp = exp, batch = i, ...)
        })
    } else {
        x <- list(batchConcentrations(exp, batches, ...))
    }


    r2Matrix <- do.call(cbind, lapply(x, `[[`, "R2"))
    concentrationMatrix <- do.call(cbind, lapply(x, `[[`, "Predicted"))
    residualMatrix <- do.call(cbind, lapply(x, `[[`, "Residuals"))
    outlierMatrix <- do.call(cbind, lapply(x, `[[`, "Outliers"))
    linearRanges <- do.call(cbind, lapply(x, `[[`, "linearRanges"))
    calRatios <- do.call(cbind, lapply(x, `[[`, "calRatios"))

    reportableCompoundsPerBatch <-
        rowData(exp)$rsdqcCorrected <= 30 &
        r2Matrix >= 0.95 &
        linearRanges >= 0.9


    if (byBatch) {
        colnames(r2Matrix) <- batches
        colnames(linearRanges) <- batches
        colnames(reportableCompoundsPerBatch) <- batches
    }


    rowData(exp)$concentrationR2 <- r2Matrix
    rowData(exp)$concentrationResiduals <- residualMatrix
    rowData(exp)$concentrationOutliers <- outlierMatrix
    rowData(exp)$linearRanges <- linearRanges
    rowData(exp)$calRatios <- calRatios

    idx <- is.na(assay(exp, "concentration"))
    rowData(exp)$hasCalculationConcentrations <- rowSums(idx) != ncol(exp)


    concentrationMatrix[concentrationMatrix < 0] <- NA

    dimnames(concentrationMatrix) <- dimnames(exp)
    assay(exp, saveAssay)[idx] <- concentrationMatrix[idx]


    return(exp)
}

#' @title Get weights for weighted regression
#' @description placeholder
#' @details placeholder
#' @returns placeholder
#' @param yValues Values of the response
getWeights <- function(yValues) {
    values <- yValues[, 1]
    if (max(values, na.rm = TRUE) == 0) {
        values <- yValues[, 2]
    }
    return(1 / (values + yValues^2))
}


#' @title Robust Linear Model with Outlier Detection
#' @description Fits a row-wise linear model between x and y, with options
#'   for intercept handling, outlier detection using studentized residuals,
#'   and reweighting.
#' @details This function calculates slope and intercept row-wise and
#'   optionally removes outliers from the first and last samples before
#'   refitting. Studentized residuals greater than 2 are considered outliers.
#' @returns A list with components:
#'   \item{slope}{Numeric vector of slopes for each row.}
#'   \item{intercept}{Numeric vector of intercepts for each row.}
#'   \item{residuals}{Matrix of residuals for each sample.}
#'   \item{outliers}{Logical matrix indicating outliers.}
#' @param xValues A numeric matrix of x-values (e.g., injections or time).
#' @param yValues A numeric matrix of y-values (e.g., intensity values).
#' @param weights A numeric matrix of weights. Used only if \code{useWeights}
#'   is TRUE.
#' @param forceOrigin Logical. If TRUE, fits the model through the origin.
#' @param subtractIntercept Logical. If TRUE, subtracts intercept from x-values.
#' @param checkForOutliers Logical. If TRUE, detects and removes outliers.
#' @param useWeights Logical. If TRUE, computes and applies weights to x-values.
#' @export
modelWithOutliers <- function(xValues, yValues, weights,
                              forceOrigin = FALSE, subtractIntercept = FALSE,
                              checkForOutliers = TRUE, useWeights = FALSE) {

    slope <- rowWiseSlope(xValues, yValues)
    intercept <- rowWiseIntercept(xValues, yValues, slope)

    if (subtractIntercept) xValues <- xValues - intercept

    residuals <- yValues - (intercept * !forceOrigin) - slope * xValues * weights
    studentResiduals <- studentResiduals(residuals, rowSums(!is.na(yValues)))
    outliers <- abs(studentResiduals) > 2

    if (!checkForOutliers) {
        return(list(slope = slope, intercept = intercept,
                    residuals = residuals, outliers = outliers))
    }

    N <- ncol(outliers)
    yValues[, 1][outliers[, 1]] <- NA
    yValues[, N][outliers[, N]] <- NA
    xValues[, 1][outliers[, 1]] <- NA
    xValues[, N][outliers[, N]] <- NA

    if (useWeights) weights <- getWeights(xValues)

    modelWithOutliers(xValues, yValues, weights, forceOrigin,
                      subtractIntercept, FALSE, useWeights)
}

#' @title Calculate concentrations per batch
#' @description placeholder
#' @details placeholder
#' @returns placeholder
#' @param exp SummarizedExperiment object
#' @param batch
#' @param assayX
#' @param assayY
#' @param type
#' @param subtractCal
#' @param checkForOutliers Boolean value, should studentized residuals be
#' used to determine outliers? If TRUE, outliers are detected and only removed
#' at the start or end of a calibration line. If FALSE, outlier detection
#' is ignored.
#' @param useWeights Boolean value, should weighted regression be used? Defaults
#' to FALSE, which indicates no weights are used.
#' @importFrom dplyr pull
#' @importFrom matrixStats rowMins rowMaxs
#' @export
batchConcentrations <- function(exp, batch, assayX = "ratio",
                                assayY = "concentration", type = "ACAL",
                                subtractCal = FALSE, forceOrigin = FALSE,
                                checkForOutliers = TRUE, subtractIntercept = FALSE,
                                useWeights = FALSE) {

    batchExp <- exp[, exp$batch %in% batch]
    training <- batchExp[, batchExp$type %in% type]

    calAliquots <- colnames(training)[order(colData(training)$calno)]
    training <- training[, calAliquots]
    yValues <- assay(training, assayY)
    xValues <- assay(training, assayX)


    cols <- seq_len(ncol(xValues) - 1)
    calRatios <- xValues[, 1 + cols] / xValues[, cols]
    calRatios[!is.finite(calRatios)] <- NA

    xValues <- xValues - xValues[, 1] * subtractCal

    if (useWeights) {
        weights <- getWeights(xValues)
    } else {
        weights <- matrix(1, nrow(xValues), ncol(xValues))
    }

    model <- modelWithOutliers(xValues, yValues, weights, forceOrigin,
                               subtractIntercept, checkForOutliers, useWeights)

    lowestCal <- rowMins(xValues, na.rm = TRUE)
    highestCal <- rowMaxs(xValues, na.rm = TRUE)

    studySamples <- assay(batchExp[, batchExp$type == "SAMPLE"], assayX)
    linearRanges <- rowSums(studySamples > lowestCal & studySamples < highestCal, na.rm = TRUE) /
        rowSums(!is.na(studySamples))

    predicted <- assay(batchExp, assayX) * model$slope + (model$intercept * !forceOrigin)
    r2 <- rowCorrelation(xValues, yValues)

    list(Predicted = predicted, R2 = r2, Residuals = model$residuals,
        Outliers = model$outliers, linearRanges = linearRanges, calRatios = calRatios
    )
}


#' @title Calculate covariance by row
#' @description placeholder
#' @details placeholder
#' @returns placeholder
#' @param x Matrix of X-values
#' @param y Matrix of y-values
rowWiseCov <- function(x, y) {
    valid <- !is.na(x) & !is.na(y) # Only consider non-NA pairs
    n <- rowSums(valid) - 1 # Degrees of freedom
    n[n <= 0] <- NA # Avoid division by zero

    rowSums((x - rowMeans(x, na.rm = TRUE)) * (y - rowMeans(y, na.rm = TRUE)),
        na.rm = TRUE
    ) / n
}

#' @title Calculate the slope of a regression line by row
#' @description
#' @details
#' @param x Matrix of X-values
#' @param y Matrix of y-values
rowWiseSlope <- function(x, y) rowWiseCov(x, y) / rowWiseCov(x, x)

#' @title Calculate the intercept of a regression line by row
#' @description placeholder
#' @details placeholder
#' @returns placeholder
#' @param x Matrix of X-values
#' @param y Matrix of y-values
#' @param slope Vector of named slope values, usually
#' obtained from [rowWiseSlope]
rowWiseIntercept <- function(x, y, slope) {
    rowMeans(y, na.rm = TRUE) - slope * rowMeans(x, na.rm = TRUE)
}


#' @title Calculate the R2 manually
#' @description placeholder
#' @details placeholder
#' @returns placeholder
rowCorrelation <- function(x, y) {
    mx <- x - rowMeans(x, na.rm = TRUE)
    my <- y - rowMeans(y, na.rm = TRUE)

    r2 <- rowSums(mx * my, na.rm = TRUE) / sqrt(rowSums(mx^2, na.rm = TRUE) * rowSums(my^2, na.rm = TRUE))
    r2[!is.finite(r2)] <- NA
    r2
}
