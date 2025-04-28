#' @title Add the ratio of compounds present in linear range
#' @description placeholder
#' @details placeholder
#' @returns placeholder
#' @param exp
#' @param assay
#' @param calType
#' @param minCalNo
#' @param maxCalNo
#' @param saveAssay
#' @importFrom dplyr bind_cols
#' @export
addLinearRange <- function(exp,
                           assay = "ratio",
                           calType = metadata(exp)$concentration,
                           minCalNo = 2,
                           maxCalNo = 6,
                           saveAssay = "CalRange") {
    if (is.null(metadata(exp)$concentration)) {
        return(exp)
    }

    if (!validateExperiment(exp)) {
        stop("Invalid experiment")
    }

    if (!calType %in% exp$type) {
        return(exp)
    }

    # Per batch, check which compounds ratios are between the minimum and maximum
    batches <- unique(exp$batch)

    # Per batch, return a boolean matrix and bind them together
    # Returns a compounds x (nTypes x nBatches) data.frame of ratios
    ratioDataFrame <- dplyr::bind_cols(lapply(batches, function(batch) {
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

    # Return the updated experiment
    return(exp)
}
