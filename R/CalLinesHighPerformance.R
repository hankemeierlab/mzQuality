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
    residualMatrix <- do.call(cbind, lapply(x, `[[`, "StudentResiduals"))
    outlierMatrix <- do.call(cbind, lapply(x, `[[`, "Outliers"))
    linearRanges <- do.call(cbind, lapply(x, `[[`, "linearRanges"))
    calRatios <- do.call(cbind, lapply(x, `[[`, "calRatios"))

    reportableCompoundsPerBatch <- rowData(exp)$rsdqcCorrected <= 30 &
        r2Matrix >= 0.95 &
        linearRanges >= 0.9


    if (byBatch) {
        colnames(r2Matrix) <- batches
        colnames(linearRanges) <- batches
        colnames(reportableCompoundsPerBatch) <- batches
    }


    rowData(exp)$concentrationR2 <- r2Matrix
    rowData(exp)$studentizedResiduals <- residualMatrix
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
#' @export
batchConcentrations <- function(exp, batch, assayX = "ratio",
                                assayY = "concentration", type = "ACAL",
                                subtractCal = FALSE, forceOrigin = FALSE,
                                checkForOutliers = TRUE, subtractIntercept = FALSE,
                                useWeights = FALSE) {
    # Select the current batch data
    batchExp <- exp[, exp$batch %in% batch]

    # Within the batch, select the type to calulate the slope for
    training <- batchExp[, batchExp$type %in% type]

    # Sort names so that the order will always be from lowest cal to highest
    # training <- training[, sort(colnames(training))]


    metadata <- colData(training) %>%
        as.data.frame()

    metadata$aliquot <- colnames(training)
    calAliquots <- metadata %>%
        arrange(.data$calno) %>%
        pull(.data$aliquot)

    training <- training[, calAliquots]

    # Take the y-values for the given assay of the training data
    yValues <- assay(training, assayY)

    # Similarly, take the x-values of the training-set
    xValues <- assay(training, assayX)

    N <- ncol(xValues)
    calRatios <- xValues[, seq(2, N, 1)] / xValues[, seq(1, N - 1, 1)]
    calRatios[!is.finite(calRatios)] <- NA

    # If subtractCAL, subtract CAL0 from everything
    # because of sorting, we can take the first column
    if (subtractCal) {
        xValues <- xValues - xValues[, 1]
    }


    weights <- matrix(1, nrow = nrow(xValues), ncol = ncol(xValues))

    if (useWeights) {
        weights <- getWeights(xValues)
    }


    # Calculate the slope given the x- and y-values
    slope <- rowWiseSlope(
        x = xValues,
        y = yValues
    )

    # Next, Calculate the intercept given the x- and y-values and slope
    int <- rowWiseIntercept(
        x = xValues,
        y = yValues,
        slope = slope
    )

    if (subtractIntercept) {
        xValues <- xValues - int
    }

    # Calculate the raw residuals by subtracting the intercept and slope * x-vals
    residuals <- yValues - (int * !forceOrigin) - slope * xValues * weights

    # Calculate the studentized-residuals (t-distributed error residuals)
    studentResiduals <- studentResiduals(
        residuals = residuals,
        n = rowSums(!is.na(yValues))
    )

    # If the values are above 2, they are cosidered outliers (according to olsrr)
    outliers <- abs(studentResiduals) > 2

    if (checkForOutliers) {
        # only remove values from the first and last column if they are outliers
        startOutliers <- outliers[, 1]
        N <- ncol(outliers)
        endOutliers <- outliers[, N]

        yValues[, 1][startOutliers] <- NA
        yValues[, N][endOutliers] <- NA

        # Repeat the same for x-values
        xValues[, 1][startOutliers] <- NA
        xValues[, N][endOutliers] <- NA

        if (useWeights) {
            weights <- getWeights(xValues)
        }

        # Recalculate the slope, now with added NA values for outliers
        # These will be ignored due to na.rm = TRUE
        slope <- rowWiseSlope(
            x = xValues,
            y = yValues
        )

        # Similarly for the intercept
        int <- rowWiseIntercept(
            x = xValues,
            y = yValues,
            slope = slope
        )

        # Calculate the raw residuals by subtracting the intercept and slope * x-vals
        # residuals <- yValues - (int * !forceOrigin) - slope * xValues
    }



    # Check if values are within linear range (90% of study stamples between lowest and highest cal)
    lowestCal <- matrixStats::rowMins(xValues, na.rm = TRUE)
    highestCal <- matrixStats::rowMaxs(xValues, na.rm = TRUE)

    studySamples <- assay(batchExp[, batchExp$type == "SAMPLE"], assayX)
    linearRanges <- rowSums(studySamples > lowestCal & studySamples < highestCal, na.rm = TRUE) / rowSums(!is.na(studySamples))

    # Take all the values of the batch-assay and multiply with the slope and
    # add the intercept to obtain the new values
    predicted <- assay(batchExp, assayX) * slope + (int * !forceOrigin)

    r2 <- rowCorrelation(xValues, yValues)

    return(
        list(
            Predicted = predicted,
            R2 = r2,
            StudentResiduals = studentResiduals,
            Outliers = outliers,
            linearRanges = linearRanges,
            calRatios = calRatios
        )
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
