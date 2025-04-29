#' @title Prepare plot data
#' @description placeholder
#' @details placeholder
#' @returns placeholder
#' @param exp
#' @param assay
#' @param batches
#' @param types
#' @param doLog
#' @importFrom dplyr filter arrange %>%
#' @importFrom imputeLCMD impute.QRILC
preparePlotData <- function(exp, assay = assayNames(exp)[1], batches = exp$batch,
                            types = exp$type, doLog = FALSE, doImpute = FALSE) {
    if (grepl("concentration", assay)) {
        exp <- exp[rowData(exp)$hasKnownConcentrations, ]
    }

    df <- reshape2::melt(unlist(t(assay(exp, assay))), value.name = assay)
    rownames(df) <- seq_len(nrow(df))
    df <- cbind(df, colData(exp)) %>%
        filter(.data$batch %in% batches) %>%
        filter(.data$type %in% types) %>%
        arrange(.data$datetime)

    df$color[df$type == "SQC"] <- "black"
    df$color[df$type == "LQC"] <- "gray"

    if (doLog) {
        df[, assay] <- log2(df[, assay])
    }

    if (doImpute) {
        df[, assay] <- imputeLCMD::impute.QRILC(as.matrix(df[, assay]))
    }

    df
}

#' @title Add labels to a plot
#' @description
#' @details
#' @returns
#' @param p
#' @param textcolor
#' @param textSize
#' @export
addLabels <- function(p, textColor = "white", textSize = 3) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_label_repel(
            color = textColor,
            segment.color = "black",
            size = textSize
        )
    }
    return(p)
}

#' @title Facet a ggplot by a column
#' @description placeholder
#' @details placeholder
#' @returns placeholder
#' @param plot A ggplot2 object
#' @param by Character vector that matches one of the properties of the
#' `colData(exp)`. Defaults to `"batch"`
#' @param ncol Number of columns, defaults to 1
#' @param shareY Boolean value, should the y-axis be shared among facets?
#' Defaults to TRUE
#' @param shareX Boolean value, should the x-axis be shared among facets?
#' Defaults to FALSE
#' @importFrom ggplot2 facet_wrap
#' @export
facetPlot <- function(plot, by = "batch", ncol = 1, shareY = TRUE, shareX = FALSE) {
    if (shareY & shareX) {
        scales <- "fixed"
    } else if (shareY) {
        scales <- "free_x"
    } else if (shareX) {
        scales <- "free_y"
    } else {
        scales <- "free"
    }


    plot + facet_wrap(facets = by, ncol = ncol, scales = scales)
}

#' @title Add text bubble to a plot
#' @description
#' @details
#' @export
addTextBubble <- function(exp, df, type = "aliquot") {
    if (type == "aliquot") {
        metadata <- as.data.frame(colData(exp[, df$aliquot]))

        df$text <- vapply(seq_len(nrow(metadata)), function(i) {
            names <- c(type, colnames(metadata))
            values <- c(rownames(metadata)[i], as.vector(metadata[i, ]))

            paste0(sprintf("%s: %s", names, values), collapse = "\n")
        }, character(1))
    } else {
        cols <- vapply(rowData(exp), function(x) !"matrix" %in% class(x), logical(1))
        x <- rowData(exp)[, as.vector(which(cols)), drop = FALSE]

        metadata <- as.data.frame(x)
        rownames(metadata) <- rownames(x)

        res <- vapply(seq_len(nrow(metadata)), function(i) {
            names <- c(type, colnames(metadata))
            values <- c(rownames(metadata)[i], as.vector(metadata[i, ]))
            paste0(sprintf("%s: %s", names, values), collapse = "\n")
        }, character(1))

        names(res) <- rownames(metadata)

        df$text <- res[as.vector(unlist(df[, type]))]
    }

    df
}

#' @title Compound Plot
#' @description
#' Creates a line plot of intensity or ratio values for a
#'   specific compound across injections, colored by sample
#'   type or batch, with optional trendlines and internal
#'   standard visualization.
#' @details This function generates a time-series-like plot for
#'   a specific compound in a `SummarizedExperiment` object. It
#'   includes support for multiple batches, sample types,
#'   internal standards, trend smoothing, and optional use of
#'   short sample names or log transformation.
#' @returns A `ggplot` object showing compound signal
#'   intensities over time, grouped and colored by type or batch.
#' @param exp A `SummarizedExperiment` object containing
#'   metabolomics data.
#' @param compound Integer or character. Row index or name of
#'   the compound to be plotted.
#' @param batches Character vector of batch identifiers to include
#'   (default is all batches in `exp`).
#' @param types Character vector of sample types to include
#'   (default is all types in `exp`).
#' @param assay Character string specifying the assay to use
#'   (default is "ratio_corrected").
#' @param doLog Logical. Whether to log-transform the y-axis
#'   values (default is TRUE).
#' @param withinTrend Logical. If TRUE, trendlines are drawn per
#'   batch within types (default is FALSE).
#' @param shortName Logical. If TRUE, short sample names are used
#'   on the x-axis (default is FALSE).
#' @param trendTypes Character vector of sample types to be used
#'   for trendline fitting (default is `metadata(exp)$QC`).
#' @param addInternalStandards Logical. Whether to add
#'   corresponding internal standards (if available)
#'   (default is FALSE).
#' @param smooth Logical. Whether to apply smoothing with trend
#'   lines (default is TRUE).
#' @importFrom dplyr filter arrange mutate %>%
#' @importFrom ggplot2 ggplot aes scale_color_manual ylab theme_minimal
#'   xlab theme element_text geom_smooth
#' @importFrom methods is
#' @importFrom stats reorder
#' @export
compoundPlot <- function(exp, compound = 1, batches = exp$batch, types = exp$type,
                         assay = "ratio_corrected", doLog = TRUE,
                         withinTrend = FALSE, shortName = FALSE,
                         trendTypes = metadata(exp)$QC,
                         addInternalStandards = FALSE, smooth = TRUE) {
    if (is(compound, "numeric")) {
        compound <- rownames(exp)[compound]
    }

    df <- preparePlotData(exp[compound, ], assay, batches = batches, types = types, doLog = FALSE)


    if (addInternalStandards | assay == "area" & metadata(exp)$hasIS) {
        df_is <- as.data.frame(preparePlotData(exp[compound, ],
            assay = "area_is",
            batches = batches, types = types,
            doLog = FALSE
        ))

        df_is$Var2 <- rowData(exp[compound, ])$compound_is
        colnames(df_is)[colnames(df_is) == "area_is"] <- assay

        df_is$type <- "ISTD"
        df_is$color <- "pink"


        df <- rbind(df, df_is)
    }



    df <- df %>%
        mutate(group = paste(
            .data$type,
            .data[[ifelse(withinTrend, "batch", "type")]],
            ifelse(grepl("CAL", .data$type), .data$injection, 1),
            sep = "_"
        ))



    colors <- setNames(unique(df$color), unique(df$type))

    if (shortName) {
        df$Var1 <- order(df$injection_time)
    }

    df <- df %>%
        arrange(.data$datetime) %>%
        mutate(aliquot = reorder(.data$Var1, .data$injection_time))


    df <- addTextBubble(exp, df, type = "aliquot")


    p <- scatterPlot(df, assay, trendTypes, colors)
    return(p)
}

#' @title Compound plot function
#' @description
#' A short description...
#' @details placeholder
#' @returns placeholder
#' @importFrom dplyr filter arrange mutate
#' @importFrom ggplot2 ggplot aes scale_color_manual ylab theme_minimal xlab
#' theme element_text geom_smooth
#' @importFrom stats lm
#' @export
scatterPlot <- function(df, assay, trendTypes, colors, smooth = TRUE, doLog = TRUE) {
    p <- ggplot(df, aes(x = .data$aliquot, y = .data[[assay]], text = .data$text)) +
        geom_point(aes(fill = .data$type), color = "black", stroke = 0.2, size = 3, shape = 21, na.rm = TRUE) +
        scale_fill_manual(values = colors) +
        ylab(assay) +
        xlab("") +
        theme_minimal() +
        ggplot2::ggtitle(df$Compound[1]) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7))

    trendData <- filter(df, .data$type %in% trendTypes)
    if (smooth && nrow(trendData) > 0) {
        p <- p + scale_color_manual(values = colors) +
            geom_smooth(
                formula = "y ~ x", inherit.aes = TRUE, method = lm, show.legend = FALSE, na.rm = TRUE,
                aes(group = .data$group, color = .data$type), level = .95,
                data = trendData
            )
    }


    if (doLog) {
        p <- p + scale_y_log10()
    }

    return(p)
}

#' @title Aliquot Plot
#' @description
#' A short description...
#' @details placeholder
#' @returns placeholder
#' @param exp
#' @param assay
#' @param batches
#' @param types
#' @param doLog
#' @importFrom dplyr filter arrange mutate
#' @importFrom ggplot2 ggplot aes scale_color_manual ylab theme_minimal xlab
#' theme element_text scale_y_log10
#' @export
aliquotPlot <- function(exp, assay, batches = exp$batch, types = exp$type,
                        doLog = FALSE) {
    df <- preparePlotData(exp, assay, batches, types = types, doLog = doLog)
    vals <- setNames(unique(df$color), unique(df$type))
    p <- df %>%
        arrange(.data$datetime) %>%
        mutate(sample = reorder(.data$Var1, .data$injection_time)) %>%
        ggplot(aes(x = .data$sample, y = .data[[assay]], text = .data$sample)) +
        geom_boxplot(aes(fill = .data$type), na.rm = TRUE) +
        scale_fill_manual(values = vals) +
        ylab(assay) +
        theme_minimal() +
        xlab("") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7))

    if (!doLog) {
        p <- p + scale_y_log10()
    }
    p
}

#' @title Violin Plot
#' @description
#' A short description...
#' @details placeholder
#' @returns placeholder
#' @param exp
#' @param assay
#' @param batches
#' @param types
#' @param doLog
#' @param addMedian
#' @param addConfidenceInterval
#' @param withinTrend
#' @importFrom dplyr filter arrange mutate
#' @importFrom ggplot2 ggplot aes scale_color_manual ylab theme_minimal xlab
#' theme element_text scale_y_log10 stat_summary geom_violin
#' @importFrom stats sd
#' @importFrom dplyr group_by ungroup summarise pull n
#' @export
violinPlot <- function(exp, assay = "ratio_corrected",
                       batches = exp$batch, types = metadata(exp)$QC,
                       doLog = FALSE, addMedian = TRUE,
                       addConfidenceInterval = TRUE,
                       withinTrend = FALSE) {
    df <- preparePlotData(exp, assay, batches, types = types, doLog = doLog) %>%
        group_by(.data$Var2, .data$batch) %>%
        mutate(order = seq_len(n())) %>%
        ungroup()

    vals <- viridis::viridis(
        length(unique(df$batch)),
        begin = 0.2,
        end = 0.8,
        option = "H"
    )
    vals <- setNames(vals, unique(df$batch))

    df$batch <- as.factor(df$batch)



    df$compound <- as.character(df$Var2)
    df <- addTextBubble(exp, df, type = "compound")


    p <- ggplot(df, aes(x = .data$Var1, y = .data[[assay]], fill = .data$batch)) +
        geom_violin(na.rm = TRUE) +
        geom_point(
            aes(text = .data$text),
            size = .2,
            position = position_jitter(0.2),
            show.legend = FALSE
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        xlab("Aliquot")

    if (!doLog) {
        p <- p + scale_y_log10()
    }

    if (addMedian) {
        p <- p + stat_summary(
            fun = "median", geom = "point", shape = 24, size = 3,
            fill = "black", show.legend = FALSE, color = "white", stroke = 0.5
        )
    }

    if (addConfidenceInterval) {
        medians <- df %>%
            group_by(.data$Var1) %>%
            summarise(median = median(.data[[assay]], na.rm = TRUE)) %>%
            pull(.data$median)

        m <- mean(medians, na.rm = TRUE)
        sd <- sd(medians, na.rm = TRUE)

        p <- p + geom_hline(yintercept = m + sd * 3, linetype = "dashed", show.legend = FALSE) +
            geom_hline(yintercept = m - sd * 3, linetype = "dashed", show.legend = FALSE)
    }
    p
}

#' @title Transform data into PCA data
#' @description
#' A short description...
#' @details placeholder
#' @returns placeholder
#' @importFrom stats prcomp var xtabs
#' @importFrom dplyr group_by mutate ungroup if_else
getPcaData <- function(exp, assay = assayNames(exp)[1],
                       batches = exp$batch, types = exp$type,
                       doLog = TRUE) {
    df <- preparePlotData(
        exp = exp, assay = assay, batches = batches,
        types = types, doLog = doLog, doImpute = FALSE
    )

    df <- df[is.finite(df[, assay]), ]

    pc <- df %>%
        group_by(.data$Var1) %>%
        # Impute infinite values
        mutate(value = if_else(condition = is.infinite(.data[[assay]]) | is.na(.data[[assay]]),
            true = min(.data[[assay]][is.finite(.data[[assay]])]) * 0.01,
            false = .data[[assay]]
        )) %>%
        ungroup() %>%
        mutate(Var1 = factor(.data$Var1))

    # Convert long to wide
    pc <- xtabs(value ~ Var1 + Var2, data = pc, subset = NULL) %>%
        as.data.frame.matrix()


    pc <- pc[, which(apply(pc, 2, var) != 0)]
    pc <- prcomp(pc, scale = TRUE, center = TRUE)

    return(pc)
}


#' @title PCA Plot of Metabolomics Data
#' @description Performs principal component analysis (PCA) on
#'   metabolomics data and visualizes selected principal components
#'   with coloring and grouping by sample type or batch.
#' @details This function calculates PCA on selected assay data from a
#'   `SummarizedExperiment` object. Users can choose which principal
#'   components to plot, whether to log-transform the data, and whether
#'   to display confidence intervals. Sample types or batches can be
#'   visualized with custom colors. Optionally, sample names are shown
#'   using label bubbles for clarity.
#' @returns A `ggplot` object representing the PCA plot with options to
#'   customize axes, coloring, and annotation.
#' @param exp A `SummarizedExperiment` object containing metabolomics data.
#' @param assay Character string specifying the assay to use for PCA
#'   (default is the first assay in the object).
#' @param batches Character vector indicating which batches to include.
#'   Defaults to all batches in the experiment.
#' @param types Character vector indicating which sample types to include.
#'   Defaults to all types in the experiment.
#' @param pc1 Integer indicating the principal component to plot on the
#'   x-axis (default is 1).
#' @param pc2 Integer indicating the principal component to plot on the
#'   y-axis (default is 2).
#' @param doLog Logical flag indicating whether to log-transform the assay
#'   data before PCA (default is TRUE).
#' @param sampleAsBatch Logical flag. If TRUE, sample points are labeled
#'   using batch numbers instead of type (default is TRUE).
#' @param addConfidenceInterval Logical flag. If TRUE, a 95% confidence
#'   ellipse is drawn around each group (default is TRUE).
#' @param typeColors Named vector of colors used for each sample type or
#'   batch (default includes black for SQC and grey for LQC).
#' @importFrom ggplot2 guides stat_ellipse
#' @importFrom dplyr select bind_cols arrange all_of mutate filter
#' @importFrom viridis viridis
#' @export
pcaPlot <- function(exp, assay = assayNames(exp)[1],
                    batches = exp$batch, types = exp$type,
                    pc1 = 1, pc2 = 2,
                    doLog = TRUE, sampleAsBatch = TRUE,
                    addConfidenceInterval = TRUE,
                    typeColors = c(SQC = "black", LQC = "grey")) {
    # PCA
    pc <- getPcaData(exp, assay, batches, types, doLog)


    # Create axis labels
    expvar <- (pc$sdev)^2 / sum(pc$sdev^2)

    exp <- exp[, exp$type %in% types & exp$batch %in% batches]
    exp <- exp[, rownames(pc$x)]

    data <- pc$x %>%
        as.data.frame() %>%
        # Select the PCs chosen
        select(all_of(c(pc1, pc2))) %>%
        # Bind column metadata
        bind_cols(as.data.frame(colData(exp))) %>%
        mutate(batch = formatC(
            as.integer(.data$batch),
            max(nchar(.data$batch)),
            format = "d", flag = "0"
        )) %>%
        arrange(.data$batch) %>%
        mutate(type = ifelse(.data$type == "SAMPLE" & sampleAsBatch,
            .data$batch, .data$type
        ))

    data$aliquot <- rownames(data)
    types <- sort(unique(data$type))
    colorVec <- viridis(length(types), begin = 0.2, end = 0.8, option = "H")

    samps <- !types %in% exp$type
    types[samps] <- paste("Batch", types[samps])
    colorVec <- setNames(colorVec, types)
    colorVec <- c(typeColors, colorVec)
    colorVec <- colorVec[!duplicated(names(colorVec))]

    samps <- !data$type %in% exp$type
    data$type[samps] <- paste("Batch", data$type[samps])
    data <- addTextBubble(exp, data, type = "aliquot")


    loadings <- data.frame(pc$rotation, variable = rownames(pc$rotation))
    scale_factor <- max(abs(data[, c(1, 2)])) / max(abs(loadings[, c(1, 2)]))
    loadings[, c(1, 2)] <- loadings[, c(1, 2)] * scale_factor

    # Make plot
    p <- ggplot(
        data = data,
        mapping = aes(
            x = .data$PC1, y = .data$PC2, fill = .data$type,
            group = .data$type, text = .data$text
        )
    ) +
        geom_point(size = 3, stroke = 0.2, shape = 21, color = "black") +
        theme_minimal() +
        xlab(sprintf("PC%s (%s%%)", pc1, format(expvar[pc1] * 100, digits = 4))) +
        ylab(sprintf("PC%s (%s%%)", pc2, format(expvar[pc2] * 100, digits = 4))) +
        scale_color_manual(values = colorVec) +
        scale_fill_manual(values = colorVec) +
        guides(color = "legend")

    if (addConfidenceInterval) {
        p <- p + stat_ellipse(
            data = data,
            mapping = aes(
                group = .data$type, color = .data$type,
                x = .data$PC1, y = .data$PC2
            ),
            inherit.aes = TRUE,
            na.rm = TRUE,
            level = 0.9
        )
    }

    return(p)
}

#' @title Concentration Plot
#' @description placeholder
#' @details placeholder
#' @returns placeholder
#' @param exp
#' @param assay
#' @param compound
#' @param batch
#' @param types
#' @param calType
#' @param removeOutliers
#' @param withinTrend
#' @export
concentrationPlot <- function(exp, assay = "ratio_corrected", compound = 1,
                              batch = 1, types = c("SAMPLE", "SQC"),
                              calType = metadata(exp)$concentration,
                              removeOutliers = TRUE, withinTrend = TRUE,
                              plotOnCalibrationLine = TRUE) {
    cals <- exp$type == calType

    df <- preparePlotData(exp[compound, ], "concentration", batch, types = c(calType, types), doLog = FALSE)
    vals <- setNames(unique(df$color), unique(df$type))

    df$ratio <- preparePlotData(exp[compound, ], assay, batch, types = c(calType, types), doLog = FALSE)[, assay]

    df$batch <- as.factor(df$batch)


    modelData <- df %>%
        filter(.data$type == calType)

    if (!plotOnCalibrationLine) {
        df$concentration[df$type %in% types] <- 0
    }


    modelData$outlier <- c(rowData(exp[compound, ])$concentrationOutliers[, unique(as.vector(modelData$Var1))])

    df$outlier <- FALSE
    df <- dplyr::bind_rows(
        df %>% filter(.data$type != calType),
        modelData
    )

    if (removeOutliers) {
        if (modelData$outlier[1]) {
            modelData <- modelData[-1, ]
        }
        if (modelData$outlier[nrow(modelData)]) {
            modelData <- modelData[-nrow(modelData), ]
        }
    }

    df$aliquot <- as.character(df$Var1)
    df <- addTextBubble(exp, df, type = "aliquot")

    ggplot(df, aes(
        x = .data$concentration, y = .data$ratio, color = .data$type,
        group = .data$batch
    )) +
        geom_point(size = 2, aes(text = .data$text)) +
        scale_color_manual(values = vals) +
        geom_hline(mapping = aes(yintercept = max(.data$ratio)), linetype = "dashed", data = modelData) +
        geom_hline(mapping = aes(yintercept = min(.data$ratio)), linetype = "dashed", data = modelData) +
        geom_smooth(
            formula = "y ~ x", method = lm, inherit.aes = withinTrend, level = .95,
            aes(x = .data$concentration, y = .data$ratio),
            na.rm = TRUE, show.legend = FALSE,
            data = modelData
        ) +
        theme_minimal()
}

#' @title Plot heatmap of RSDQCs with internal standards
#' @description This function plots a heatmap of batch-corrected RSDQCs
#' of all compound - internal standard combinations. This allows for insight
#' if the picked internal standard is the correct one or if possible another
#' internal standard should be picked.
#' @param exp SummarizedExperiment object
#' @returns Plotly heatmap of compound - internal standard RSDQCs
#' @export
#' @examples
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Plot RSDQCs
#' rsdqcPlot(exp)
rsdqcPlot <- function(exp) {
    if (!validateExperiment(exp)) {
        stop("Invalid experiment")
    }

    if (!metadata(exp)$hasIS) {
        return(NULL)
    }

    exp <- exp[!is.na(rowData(exp)$rsdqcCorrected), ]
    rsdqcs <- matrixRSDQCs(exp, returnRSDs = TRUE)

    cols <- viridis::viridis(
        n = 256, alpha = 1, begin = 0.2, end = 0.8,
        option = "H"
    )
    p <- heatmaply::heatmaply(as.matrix(rsdqcs),
        plot_method = "ggplot",
        dendrogram = "none",
        margins = c(0, 0, 0, 0), colors = cols
    )
    return(p)
}

#' @title Plot heatmap of the given assay
#' @description This functon plots a heatmap of the given assay. Constructed
#' using heatmaply, it is an interactive plot that shows the value of the
#' provide assay.
#' @param exp SummarizedExperiment object
#' @param assay Which assay to use for plotting the values
#' @returns Plotly heatmap of values in the provided assay
#' @export
#' @examples
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Plot values in heatmap
#' heatmapPlot(exp)
heatmapPlot <- function(exp, assay = "ratio") {
    if (!validateExperiment(exp)) {
        stop("Invalid experiment")
    }
    df <- assay(exp, assay)
    df[is.infinite(df)] <- NA
    if (all(is.na(df))) {
        return(NULL)
    }

    cols <- viridis::viridis(
        n = 256, alpha = 1, begin = 0.2, end = 0.8,
        option = "H"
    )

    p <- heatmaply::heatmaply(
        df,
        dendrogram = "none", plot_method = "ggplot",
        margins = c(0, 0, 0, 0), colors = cols
    )

    return(p)
}




#' @title Create calibration plot
#' @description The calibration plot shows the CAL Lines, Samples, and QCs of
#' a single compound in a different view to the compoundPlot. The CAL lines
#' are here shown as lines, which aids in determining if the samples are
#' between the lowest and highest CAL line.
#' @param exp SummarizeExperiment object
#' @param assay which assay to plot
#' @param compound Which compound to plot
#' @param batch Batch number
#' @param guides Which guides to add. Options available are c("95% CI",
#' "Limit of Blank", "Limit of Detection", "Limit of Quantification")
#' @returns ggplot object of a scatterplot with "SAMPLE", "QC" and "CAL" types.
#' @importFrom ggplot2 xlim ggplot geom_hline aes
#' @importFrom stats complete.cases median sd
#' @export
#' @examples
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Create scatterplot with calibration lines
#' calibrationPlot(exp)
calibrationPlot <- function(exp, assay = "ratio_corrected",
                            compound = 1, batch = 1, guides = NULL) {
    if (!validateExperiment(exp)) {
        stop("Invalid experiment")
    }

    if (!any(exp$type %in% c("CAL", "ACAL")) || !"calno" %in%
        colnames(colData(exp))) {
        return(NULL)
    }

    qc <- metadata(exp)$QC
    types <- c(qc, "SAMPLE", "CAL", "ACAL")

    sub <- exp[compound, exp$type %in% types & exp$batch == batch]

    dl <- as.data.frame(cbind(
        calno = as.factor(as.integer(sub$calno)),
        type = sub$type, reshape2::melt(assay(sub, assay))
    ))


    dl$order <- seq_len(nrow(dl))
    points <- dl[dl$type %in% c(qc, "SAMPLE"), ]
    points <- points[, -which(colnames(points) == "calno")]
    lines <- dl[dl$type %in% c("CAL", "ACAL"), ]

    lines <- lines[complete.cases(lines), ]
    if (nrow(lines) == 0) {
        return(ggplot() +
            theme_bw())
    }

    dl <- dl[, -which(colnames(dl) == "calno")]
    dl <- dl[complete.cases(dl), ]

    # points$text <- getAliquotPlotLabels(aliquotFilter(sub, type = c(qc, "SAMPLE")))


    p <- ggplot(points, aes(
        x = .data$order, y = .data$value,
        fill = .data$type, text = .data$text
    )) +
        geom_hline(aes(
            yintercept = .data$value,
            color = .data$calno
        ), data = lines) +
        geom_point(shape = 21, size = 3) +
        xlim(min(dl$order,
            na.rm = TRUE
        ), max(dl$order, na.rm = TRUE)) +
        theme_bw() +
        xlab("Injection order") +
        ylab(assay)

    if ("95% CI" %in% guides) {
        df <- dl[dl$type == qc, ]

        p <- p + geom_hline(yintercept = median(df$value) +
            sd(df$value) * 2, linetype = "dashed") +
            geom_hline(yintercept = mean(df$value) -
                sd(df$value) * 2, linetype = "dashed")
    }

    p
}

#' @title Plot the Relative Standard Deviation
#' @description This plot shows the RSD of compounds across batches. This can
#' aid in identifying which compounds have a high variation in the given assay
#' in certain batches.
#' @param assay Which assay to use for calculating the RSDQC
#' @param exp SummarizedExperiment object
#' @param qc type of QC to plot
#' @param number Optional Show the n number of worst RSD
#' @returns ggplot object of a scatter- / lineplot of RSDs per compound
#' and batch
#' @export
#' @importfrom ggplot2 geom_line theme geom_point scale_color_manual scale_fill_manual
#' @examples
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Create scatterplot with Relative Standard Deviation values
#' rsdPlot(exp, number = 10)
rsdPlot <- function(exp, assay = "ratio_corrected",
                    qc = "SQC", number = nrow(exp)) {
    if (!validateExperiment(exp)) {
        stop("Invalid experiment")
    }
    exp <- exp[!is.na(rowData(exp)$rsdqcCorrected), ]

    exp <- exp[utils::tail(
        order(rowData(exp)$rsdqcCorrected),
        as.integer(number)
    ), ]
    df <- do.call(cbind, lapply(unique(exp$batch), function(x) {
        rsdqc(exp[, exp$type == qc & exp$batch == x], assay)
    }))

    colnames(df) <- unique(exp$batch)
    df <- reshape2::melt(df)
    colnames(df) <- c("compound", "batch", "RSD")
    df$batch <- factor(df$batch, levels = unique(df$batch))


    typeColors <- setNames(viridis::viridis(length(unique(df$compound)),
        option = "H"
    ), unique(df$compound))

    df$text <- getCompoundPlotLabels(exp)

    ggplot(df, aes(
        x = .data$batch, y = .data$RSD, fill = .data$compound,
        group = .data$compound, color = .data$compound, text = .data$text
    )) +
        geom_point(pch = 21, stroke = 0.1, size = 3) +
        ggplot2::geom_line() +
        theme_bw() +
        theme(axis.text.x = element_text(
            angle = 45,
            hjust = 1
        )) +
        scale_fill_manual(values = typeColors) +
        scale_color_manual(values = typeColors) +
        ylab("RSDQC") +
        xlab("")
}


#' @title Retrieve compound labels for text hovering
#' @param exp SummarizedExperiment object
#' @noRd
getCompoundPlotLabels <- function(exp) {
    rowData <- as.data.frame(rowData(exp))
    rowData <- rowData[, grep("BLANK|PROC|Suggested", colnames(rowData),
        invert = TRUE, value = TRUE
    )]
    vapply(seq_len(nrow(exp)), function(i) {
        sprintf(
            "\rCompound: %s\n%s", rownames(exp)[i],
            paste0(sprintf(
                "\r%s: %s\n", colnames(rowData),
                rowData[i, ]
            ), collapse = "")
        )
    }, character(1))
}

#' @title Plot an assay per batch with point distribution
#' @param exp SummarizedExperiment object
#' @param assay Name or index of the assay to plot
#' @param compound Compound to plot
#' @returns ggplot object of a Batch Assay Plot
#' @importFrom ggplot2 ggplot aes geom_point geom_boxplot position_jitter
#' @export
#' @examples
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Create Batch Assay Plot
#' batchAssayPlot(exp)
batchAssayPlot <- function(exp, assay = 1, compound = 1) {
    if (!validateExperiment(exp)) {
        stop("Invalid experiment")
    }
    exp <- exp[compound, ]

    df <- reshape2::melt(assay(exp, assay))
    df <- df[order(df$Var2, df$Var1), ]

    labels <- paste(as.data.frame(colData(exp)), collapse = "\n")
    df$text <- sprintf("%s\nValue: %s", labels, round(df$value, 4))

    df$batch <- as.factor(exp$batch)

    ggplot(df, aes(
        x = .data$batch, y = .data$value,
        fill = .data$batch, text = .data$text
    )) +
        geom_boxplot() +
        geom_point() +
        theme_bw() +
        theme(axis.text.x = element_text(
            angle = 90,
            hjust = 1
        ), legend.position = "bottom")
}

#' @title Plot the batch correction factors
#' @description The batch correction plot shows a boxplot of correction factors
#' per batch. These are the calculated multipliers, which ideally should be
#' as close to 1. A value near 1 indicates the batch is similar to the median
#' of all batches. This plot can be used to see if compounds are
#' over-corrected, or if a batch has bad measurements.
#' @param exp SummarizedExperiment Object
#' @returns ggplot object of a boxplot of batch-effect correction factors
#' @importFrom ggplot2 ggplot aes geom_boxplot theme_bw xlab ylab theme
#' element_text .data scale_fill_manual
#' @export
#' @examples
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Create boxplot with batch-correction values
#' batchCorrectionPlot(exp)
batchCorrectionPlot <- function(exp) {
    if (!validateExperiment(exp)) {
        stop("Invalid experiment")
    }

    sub <- exp[, which(!duplicated(exp$batch))]
    df <- assay(sub, "ratio") / assay(sub, "ratio_corrected")

    colnames(df) <- unique(exp$batch)
    df <- reshape2::melt(df)
    df$text <- getCompoundPlotLabels(sub)

    df <- df[complete.cases(df), ]
    df <- df[order(df$Var2), ]
    df$batch <- paste("Batch", df$Var2)

    x <- unique(df$batch)
    df$batch <- factor(df$batch, levels = x, ordered = TRUE)

    typeColors <- setNames(viridis::viridis(length(x),
        begin = 0.2, end = 0.8, option = "H"
    ), x)


    ggplot(df, aes(
        x = .data$batch, y = .data$value,
        fill = .data$batch, text = .data$text
    )) +
        geom_boxplot(outlier.shape = NA, outlier.size = 0) +
        geom_point(position = position_jitter(0.1)) +
        theme_bw() +
        xlab("") +
        ylab("Correction Factor") +
        scale_fill_manual(values = typeColors) +
        theme(axis.text.x = element_text(
            angle = 90,
            hjust = 1
        ), legend.position = "bottom")
}
