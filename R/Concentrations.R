#' @title Studentized residuals
#' @description Studentized residuals represent t-distributed residuals,
#' meaning that each residual is scored according to the t-distribution.
#' This allows for determining outliers in regression lines
#' @details This function effectively calculates the studentized residuals
#' without recalculating the model. The basic idea is to delete the observations
#' one at a time, each time refitting the regression model on the remaining n–1
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
#' @seealso https://online.stat.psu.edu/stat462/node/247/ for the
#' formula used here
studentResiduals <- function(residuals, n, k = 1) {
    rsd <- sqrt(rowSums(residuals^2, na.rm = TRUE) / (n - 2))
    stdResid <- residuals / rsd
    resid <- stdResid * ((n - k - 2) / (n - k - 1 - stdResid ** 2 )) ** 0.5
    resid[!is.finite(resid)] <- NA
    return(resid)
}

#' @title Calculate Concentrations for a given type in an Experiment
#' @description This function calculates predicted concentrations for all
#'   compounds in a SummarizedExperiment object. It processes each batch
#'   individually, aggregates results, and updates the experiment object
#'   with calculated concentrations and associated metadata.
#' @details The function computes concentrations for each batch in the
#'   experiment using the `batchConcentrations` function. It aggregates
#'   results such as R² values, residuals, outliers, linear ranges, and
#'   calibration ratios. The function also evaluates the reportability of
#'   compounds based on RSD, linearity, and linear range thresholds.
#'   Concentrations below zero are replaced with `NA`.
#' @returns A SummarizedExperiment object with predicted concentrations
#'   stored in the assay specified by `saveAssay` (default: `"concentration"`).
#'   Additional metadata is added to `rowData`, including R² values, residuals,
#'   outliers, linear ranges, and calibration ratios.
#' @param exp A SummarizedExperiment object containing the experimental data.
#' @param type A character vector specifying the calibration
#'   line types to use for concentration calculations.
#' @param ... Additional arguments passed to the `batchConcentrations` function.
calculateConcentrationsPerType <- function(exp, type, ...){

    batches <- unique(exp$batch)
    x <- lapply(batches, function(batch) {
        batchConcentrations(exp = exp, batch = batch, type = type, ...)
    })

    components <- c("R2", "Residuals", "Outliers", "linearRanges", "Predicted")

    matrices <- lapply(components, function(component) {
        do.call(cbind, lapply(x, `[[`, component))
    })
    names(matrices) <- components

    for (component in components[seq_len(3)]) {
        assay(exp, sprintf("%s_%s", type, component)) <- matrices[[component]]
    }

    goodRSD <-  rowData(exp)$rsdqcCorrected <= 30
    goodLinearity <- rowSums(matrices$R2 >= 0.95, na.rm = TRUE) >= 1
    goodLinearRanges <- rowSums(matrices$linearRanges >= 0.9, na.rm = TRUE) >= 1
    reportable <- goodRSD & goodLinearity & goodLinearRanges
    rowData(exp)$ReportableConcentration <- reportable

    concentrationMatrix <- matrices$Predicted
    idx <- is.na(concentrationMatrix)
    rowData(exp)$hasCalculationConcentrations <- rowSums(idx) != ncol(exp)
    concentrationMatrix[concentrationMatrix < 0] <- NA


    cols <- exp$type == type
    saveAssay <- sprintf("%s_%s", type, "concentration")
    concentrationMatrix[, cols] <- assay(exp[, cols], "concentration")
    assay(exp, saveAssay) <- concentrationMatrix

    return(exp)
}

#' @title Calculate Concentrations for a Full Experiment
#' @description This function calculates predicted concentrations for all
#'   compounds in a SummarizedExperiment object. It processes each batch
#'   individually, aggregates results, and updates the experiment object
#'   with calculated concentrations and associated metadata.
#' @details The function computes concentrations for each batch in the
#'   experiment using the `batchConcentrations` function. It aggregates
#'   results such as R² values, residuals, outliers, linear ranges, and
#'   calibration ratios. The function also evaluates the reportability of
#'   compounds based on RSD, linearity, and linear range thresholds.
#'   Concentrations below zero are replaced with `NA`.
#' @returns A SummarizedExperiment object with predicted concentrations
#'   stored in the assay specified by `saveAssay` (default: `"concentration"`).
#'   Additional metadata is added to `rowData`, including R² values, residuals,
#'   outliers, linear ranges, and calibration ratios.
#' @param exp A SummarizedExperiment object containing the experimental data.
#' @param calibrationLineTypes A character vector specifying the calibration
#'   line types to use for concentration calculations. Defaults to the
#'   concentration given in the metadata of the experiment
#' @param ... Additional arguments passed to the `batchConcentrations` function.
#' @export
#' @examples
#' # Example usage:
#' # Load example dataset
#' exp <- readRDS(system.file(package = "mzQuality", "data.RDS"))
#'
#' # Calculate concentrations
#' exp <- calculateConcentrations(exp)
#'
#' # Access the calculated concentrations
#' assay(exp, "concentration")
calculateConcentrations <- function(
        exp, calibrationLineTypes = metadata(exp)$concentration, ...
) {


    # Check if the input is a SummarizedExperiment object
    if (!validateExperiment(exp)) {
        stop("Input must be a SummarizedExperiment object.")
    }

    # Check if the assay is present in the experiment
    if (!"concentration" %in% assayNames(exp)) {
        return(exp)
    }

    # Check if the type is present in the experiment
    if (!any(calibrationLineTypes %in% exp$type)) {
        return(exp)
    }


    for (type in calibrationLineTypes) {
        exp <- calculateConcentrationsPerType(
            exp = exp, type = type, ...
        )
    }

    return(exp)
}

#' @title Get Weights for Weighted Regression
#' @description Computes weights for a weighted regression model based on the
#'   provided response values. The weights are inversely proportional to the
#'   values and their squared magnitude.
#' @details This function determines the appropriate weights for regression
#'   by evaluating the response values. If the first column of the response
#'   matrix contains only zeros, the second column is used instead.
#' @returns A numeric vector of weights for the regression model.
#' @param yValues A numeric matrix of response values.
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
modelWithOutliers <- function(
        xValues, yValues, weights,
        forceOrigin = FALSE, subtractIntercept = FALSE,
        checkForOutliers = TRUE, useWeights = FALSE
) {

    slope <- rowWiseSlope(xValues, yValues)
    int <- rowWiseIntercept(xValues, yValues, slope)

    if (subtractIntercept) xValues <- xValues - int

    residuals <- yValues - (int * !forceOrigin) - slope * xValues * weights
    studentResiduals <- studentResiduals(residuals, rowSums(!is.na(yValues)))
    outliers <- abs(studentResiduals) > 2

    if (!checkForOutliers) {
        return(list(
            slope = slope, intercept = int,
            residuals = residuals, outliers = outliers
        ))
    }

    N <- ncol(outliers)
    yValues[, 1][outliers[, 1]] <- NA
    yValues[, N][outliers[, N]] <- NA
    xValues[, 1][outliers[, 1]] <- NA
    xValues[, N][outliers[, N]] <- NA

    if (useWeights) weights <- getWeights(xValues)

    model <- modelWithOutliers(
        xValues, yValues, weights, forceOrigin,
        subtractIntercept, FALSE, useWeights
    )

    return(model)
}

#' @title Calculate absolute concentrations using calibration curves
#' @description Calculates absolute concentrations of analytes using linear
#'     models based on analyte/internal standard ratios from calibration
#'     samples.
#' @details This function performs calibration and quantification of analytes
#'     in a batch-wise manner. It extracts calibration samples from the
#'     provided SummarizedExperiment object, builds linear models for each
#'     analyte, and applies these models to predict concentrations for all
#'     samples in the specified batch.
#'
#'     The function supports various options for calibration model building,
#'     including:
#'     - Outlier detection using studentized residuals
#'     - Weighted regression for heteroscedastic data
#'     - Forcing the model through the origin
#'     - Subtracting baseline values from calibrators
#'
#'     The function also calculates linear ranges for each analyte by checking
#'     whether sample ratios fall within the calibration range.
#' @param exp SummarizedExperiment object containing the analytical data
#' @param batch Character vector specifying which batch(es) to process
#' @param assayX String specifying the assay name for the ratio values (e.g.,
#'     analyte/IS ratio)
#' @param assayY String specifying the assay name for the concentration values
#' @param type Character vector specifying which sample types to use for
#'     calibration, defaults to "ACAL"
#' @param subtractCal Logical, whether to subtract the first calibrator value
#'     from all calibrators
#' @param forceOrigin Logical, whether to force the regression line through the
#'     origin (0,0)
#' @param checkForOutliers Logical, whether to detect outliers using
#'     studentized residuals. If TRUE, outliers are detected and only removed
#'     at the start or end of a calibration line. If FALSE, outlier detection
#'     is ignored.
#' @param subtractIntercept Logical, whether to subtract the intercept from
#'     predictions
#' @param useWeights Logical, whether to use weighted regression. Defaults to
#'     FALSE, which indicates no weights are used.
#' @return A list containing:
#'     \item{Predicted}{Matrix of predicted concentrations for all samples in
#'         the batch}
#'     \item{R2}{Vector of R-squared values for each analyte's calibration
#'         curve}
#'     \item{Residuals}{Matrix of residuals from the calibration models}
#'     \item{Outliers}{Matrix indicating outliers in the calibration data
#'         (TRUE/FALSE)}
#'     \item{linearRanges}{Vector indicating the proportion of samples within
#'         the calibration range for each analyte}
#'     \item{calRatios}{Matrix of ratios between consecutive calibration
#'         levels}
#' @importFrom dplyr pull
#' @importFrom matrixStats rowMins rowMaxs
batchConcentrations <- function(
        exp, batch, assayX = "ratio",
        assayY = "concentration", type = "ACAL",
        subtractCal = FALSE, forceOrigin = FALSE,
        checkForOutliers = TRUE, subtractIntercept = FALSE,
        useWeights = FALSE
) {

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
    weights <- 1 * (getWeights(xValues) * useWeights)

    model <- modelWithOutliers(
        xValues, yValues, weights, forceOrigin,
        subtractIntercept, checkForOutliers, useWeights
    )

    lowestCal <- rowMins(xValues, na.rm = TRUE)
    highestCal <- rowMaxs(xValues, na.rm = TRUE)

    samples <- assay(batchExp[, batchExp$type == "SAMPLE"], assayX)
    boolMat <- samples > lowestCal & samples < highestCal
    linearRanges <- rowSums(boolMat, na.rm = TRUE) / rowSums(!is.na(samples))

    ax <- assay(batchExp, assayX) * model$slope
    b <- model$intercept * !forceOrigin
    predicted <- ax + b
    r2 <- rowCorrelation(xValues, yValues)

    return(getConcentrationAssays(
        batchExp, colnames(training), model, predicted, r2,
        linearRanges, calRatios
    ))
}

getConcentrationAssays <- function(
        batchExp, trainingAliquots, model, predicted, r2,
        linearRanges, calRatios
) {
    template <- assay(batchExp)
    template[TRUE] <- NA

    r2Matrix <- linearRangeMatrix <- residMatrix <- outlierMatrix <- template

    r2Matrix[, trainingAliquots] <- r2
    linearRangeMatrix[, trainingAliquots] <- linearRanges
    residMatrix[, trainingAliquots] <- model$residuals
    outlierMatrix[,trainingAliquots] <- model$outliers

    list(
        Predicted = predicted, R2 = r2Matrix, Residuals = residMatrix,
        Outliers = outlierMatrix, linearRanges = linearRangeMatrix
    )
}


#' @title Calculate covariance by row
#'
#' @description
#' This function computes the covariance for each row of two matrices,
#' handling missing values appropriately. It calculates the pairwise
#' covariance between corresponding rows of `x` and `y`.
#'
#' @details
#' The function removes `NA` values from each row-wise computation and
#' applies the standard covariance formula with an unbiased estimator
#' (dividing by `n-1`). If a row has less than two non-NA values,
#' the covariance is set to `NA` to prevent division errors.
#'
#' @param x A numeric matrix of X-values.
#' @param y A numeric matrix of Y-values.
#'
#' @returns A numeric vector containing the row-wise covariance values.
rowWiseCov <- function(x, y){

    valid <- !is.na(x) & !is.na(y)  # Only consider non-NA pairs
    n <- rowSums(valid) - 1  # Degrees of freedom
    n[n <= 0] <- NA  # Avoid division by zero

    x <- (x - rowMeans(x, na.rm = TRUE)) * (y - rowMeans(y, na.rm = TRUE))
    return(rowSums(x, na.rm = TRUE) / n)
}


#' @title Calculate the slope of a regression line by row
#'
#' @description
#' Computes the row-wise slope of a simple linear regression model,
#' where each row represents an independent dataset.
#'
#' @details
#' The function calculates the slope as the ratio of the row-wise
#' covariance of `x` and `y` to the row-wise covariance of `x` with itself.
#' This follows the standard linear regression formula:
#' `slope = Cov(x, y) / Cov(x, x)`.
#'
#' @param x A numeric matrix of X-values.
#' @param y A numeric matrix of Y-values.
#'
#' @returns A numeric vector containing the row-wise slope values.
rowWiseSlope <- function(x, y) {
    slope <- rowWiseCov(x, y) / rowWiseCov(x, x)
    return(slope)
}


#' @title Calculate the intercept of a regression line by row
#'
#' @description
#' Computes the row-wise intercept of a simple linear regression model,
#' where each row represents an independent dataset.
#'
#' @details
#' The function calculates the intercept using the standard formula:
#' `intercept = mean(y) - slope * mean(x)`, applied row-wise. The slope
#' values should be computed using [rowWiseSlope].
#'
#' @param x A numeric matrix of X-values.
#' @param y A numeric matrix of Y-values.
#' @param slope A numeric vector of slopes, typically obtained from
#' [rowWiseSlope].
#'
#' @returns A numeric vector containing the row-wise intercept values.
rowWiseIntercept <- function(x, y, slope) {
    int <- rowMeans(y, na.rm = TRUE) - slope * rowMeans(x, na.rm = TRUE)
    return(int)
}


#' @title Calculate row-wise correlation (R-squared)
#' @description Computes the row-wise correlation (R-squared) between two
#'     matrices or data frames.
#' @details This function calculates the row-wise correlation (R-squared)
#'     values between corresponding rows of two matrices or data frames. It
#'     subtracts the row means from each element, computes the dot product
#'     between the rows, and normalizes by the product of the row-wise
#'     standard deviations.
#'
#'     Non-finite values in the result (e.g., due to division by zero) are
#'     replaced with `NA`.
#' @param x A numeric matrix or data frame.
#' @param y A numeric matrix or data frame with the same dimensions as `x`.
#' @return A numeric vector containing the row-wise R-squared values.
#' @examples
#' x <- matrix(rnorm(20), nrow = 5)
#' y <- matrix(rnorm(20), nrow = 5)
#' rowCorrelation(x, y)
#' @export
rowCorrelation <- function(x, y) {
    mx <- x - rowMeans(x, na.rm = TRUE)
    my <- y - rowMeans(y, na.rm = TRUE)

    frac <- rowSums(mx * my, na.rm = TRUE) / sqrt(rowSums(mx^2, na.rm = TRUE))
    r2 <- frac * rowSums(my^2, na.rm = TRUE)
    r2[!is.finite(r2)] <- NA
    return(r2)
}

#' @title Add the ratio of compounds present in the linear range
#' @description Adds a boolean matrix to the experiment object indicating
#'     whether the ratios of compounds fall within the specified calibration
#'     range for each batch.
#' @details This function checks whether the ratios of compounds in the
#'     experiment object fall within the specified calibration range (defined
#'     by `minCalNo` and `maxCalNo`) for each batch. It creates a boolean
#'     matrix indicating whether each compound is within the range and saves
#'     it as a new assay in the experiment object.
#' @param exp A SummarizedExperiment object containing the experimental data.
#' @param assay A string specifying the assay name to use for ratio values.
#'     Defaults to "ratio".
#' @param calType A string specifying the calibration type to use. Defaults to
#'     the concentration metadata of the experiment.
#' @param minCalNo An integer specifying the minimum calibration number to
#'     include in the range. Defaults to 2.
#' @param maxCalNo An integer specifying the maximum calibration number to
#'     include in the range. Defaults to 6.
#' @param saveAssay A string specifying the name of the assay to save the
#'     boolean matrix. Defaults to "CalRange".
#' @return A SummarizedExperiment object with the updated assay containing the
#'     boolean matrix of compounds within the linear range.
#' @importFrom dplyr bind_cols
#' @examples
#' # Example usage:
#' # exp <- SummarizedExperiment(...)
#' # updatedExp <- addLinearRange(exp, assay = "ratio", minCalNo = 2,
#' #     maxCalNo = 6, saveAssay = "CalRange")
#' @export
addLinearRange <- function(
        exp, assay = "ratio", calType = metadata(exp)$concentration,
        minCalNo = 2, maxCalNo = 6, saveAssay = "CalRange"
) {

    if (!calType %in% exp$type) {
        return(exp)
    }

    if (!"calno" %in% colnames(colData(exp))) {
        return(exp)
    }

    # Per batch, check which compounds ratios are between the minimum
    # and maximum
    batches <- unique(exp$batch)

    # Per batch, return a boolean matrix and bind them together
    # Returns a compounds x (nTypes x nBatches) data.frame of ratios
    ratioDataFrame <- bind_cols(lapply(batches, function(batch) {
        # Subset the ratios per batch and calType used
        batchExp <- exp[, exp$batch == batch]

        # Select the middle cal points
        idx <- which(batchExp$calno >= minCalNo & batchExp$calno <= maxCalNo)
        middle <- batchExp[, idx]

        # Should only be 1 column, ensure with [,1]
        lowest <- middle[, which.min(middle$calno)][, 1]
        highest <- middle[, which.max(middle$calno)][, 1]

        # Select the assay in the middle callines as vectors
        minimumRatios <- as.vector(assay(lowest, assay))
        maximumRatios <- as.vector(assay(highest, assay))

        ratios <- assay(batchExp, assay)
        ratios[is.na(ratios)] <- -1

        boolMatrix <- (ratios >= minimumRatios) + (ratios <= maximumRatios) == 2
        return(boolMatrix)
    }))

    ratioDataFrame <- as.data.frame(ratioDataFrame[colnames(exp)])
    dimnames(ratioDataFrame) <- dimnames(exp)

    assay(exp, saveAssay) <- ratioDataFrame
    return(exp)
}
