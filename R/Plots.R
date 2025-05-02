#' @title Prepare plot data
#' @description Prepares a data frame for plotting by reshaping assay data
#'     from a SummarizedExperiment object into a long format and adding
#'     metadata.
#' @details This function reshapes assay data into a long format, merges it
#'     with metadata, and filters the data based on specified batches and
#'     types. It also supports optional log transformation and imputation of
#'     missing values. A color column is added to distinguish between sample
#'     types (e.g., SQC and LQC).
#' @param exp A SummarizedExperiment object containing the experimental data.
#' @param assay A string specifying the assay to use for plotting. Defaults to
#'     the first assay in `assayNames(exp)`.
#' @param batches A vector specifying the batches to include in the plot.
#'     Defaults to all batches in `exp$batch`.
#' @param types A vector specifying the sample types to include in the plot.
#'     Defaults to all types in `exp$type`.
#' @param logTransform Logical, whether to apply a log2 transformation to the
#' assay values. Defaults to `FALSE`.
#' @return A data frame in long format with columns for aliquots, compounds,
#'     assay values, metadata, and a color column for sample types.
#' @importFrom dplyr filter arrange %>%
#' @noRd
.preparePlotData <- function(
        exp, assay, batches = exp$batch,
        types = exp$type, logTransform = FALSE
) {

    # Convert assay data to long format using base R
    assayData <- assay(exp, assay)


    df <- data.frame(
        aliquot = rep(colnames(assayData), times = nrow(assayData)),
        compound = rep(rownames(assayData), each = ncol(assayData)),
        value = c(unlist(t(assayData)))
    )
    colnames(df)[3] <- assay

    # Add metadata and filter
    metadata <- as.data.frame(colData(exp))
    df <- merge(df, metadata, by.x = "aliquot", by.y = "row.names")
    df <- df[df$batch %in% batches & df$type %in% types, ]
    df <- df[order(df$datetime), ]

    # Apply log transformation if requested
    if (logTransform) {
        df[, assay] <- log2(df[, assay])
    }

    df$batch <- as.factor(df$batch)

    return(df)
}

#' @title Add labels to a plot
#' @description Adds labels to a ggplot object using ggrepel.
#' @details This function uses ggrepel to add labels to a ggplot object. It
#'     requires the ggrepel package to be installed.
#' @returns A ggplot object with added labels.
#' @param p A ggplot object to which labels will be added.
#' @param textColor A character string specifying the color of the text.
#'     Defaults to "white".
#' @param textSize A numeric value specifying the size of the text.
#'     Defaults to 3.
#' @noRd
.addPcaLabels <- function(p, textColor = "white", textSize = 3) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_label_repel(
            mapping = aes(label = .data$aliquot),
            color = textColor,
            segment.color = "black",
            size = textSize
        )
    }
    return(p)
}

#' @title Facet a ggplot by a column
#' @description Facets a ggplot object by a specified column.
#' @details This function uses ggplot2's facet_wrap to create facets for a
#'     ggplot object based on a specified column. It allows customization of
#'     the number of columns and whether axes are shared.
#' @returns A faceted ggplot object.
#' @param plot A ggplot2 object to be faceted.
#' @param by A character vector specifying the column to facet by. Defaults
#'     to "batch".
#' @param ncol An integer specifying the number of columns in the facet.
#'     Defaults to 1.
#' @param shareY A logical value indicating whether the y-axis should be
#'     shared among facets. Defaults to TRUE.
#' @param shareX A logical value indicating whether the x-axis should be
#'     shared among facets. Defaults to FALSE.
#' @importFrom ggplot2 facet_wrap
#' @export
#' @examples
#' # Read data
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Do Analysis
#' exp <- doAnalysis(exp, doAll = TRUE)
#'
#' # Create a compound plot for the first compound
#' p <- compoundPlot(
#'    exp,
#'    compound = 1,
#'    types = c("SQC", "LQC"),
#'    assay = "ratio_corrected",
#'    logTransform = TRUE
#' )
#'
#' facetPlot(p, by = "batch", ncol = 2, shareY = TRUE)
facetPlot <- function(
        plot, by = "batch", ncol = 1, shareY = TRUE, shareX = FALSE
) {
    if (shareY & shareX) {
        scales <- "fixed"
    } else if (shareY) {
        scales <- "free_x"
    } else if (shareX) {
        scales <- "free_y"
    } else {
        scales <- "free"
    }

    plot <- plot +
        facet_wrap(facets = by, ncol = ncol, scales = scales)

    return(plot)
}

#' @title Add text bubbles to a data frame for visualization
#' @description Adds detailed metadata as text bubbles to a data frame for
#'     use in visualizations, based on either aliquot-level or compound-level
#'     metadata from a SummarizedExperiment object.
#' @details This function generates text annotations for each row in the
#'     provided data frame (`df`) by extracting metadata from the
#'     SummarizedExperiment object (`exp`). The annotations are added as a
#'     new column (`text`) in the data frame. The function supports both
#'     aliquot-level and compound-level metadata, depending on the value of
#'     the `type` parameter.
#' @param exp A SummarizedExperiment object containing the experimental data.
#' @param df A data frame to which the text bubbles will be added. It must
#'     include a column matching the `type` parameter (e.g., "aliquot" or
#'     "compound").
#' @param type A string specifying the type of metadata to use for the text
#'     bubbles. Can be `"aliquot"` for aliquot-level metadata or any other
#'     string for compound-level metadata. Defaults to `"aliquot"`.
#' @return A data frame with an additional column, `text`, containing the
#'     generated text bubbles.
#' @noRd
.addTextHover <- function(exp, df, type = "aliquot") {
    if (type == "aliquot") {
        metadata <- as.data.frame(colData(exp[, df$aliquot]))

        df$text <- vapply(seq_len(nrow(metadata)), function(i) {
            names <- c(type, colnames(metadata))
            values <- c(rownames(metadata)[i], as.vector(metadata[i, ]))

            paste0(sprintf("%s: %s", names, values), collapse = "\n")
        }, character(1))
    } else {

        x <- cleanDataframe(rowData(exp))

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
#' @param logTransform Logical. Whether to log-transform the y-axis
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
#' @param addText Logical. Whether to add text bubbles to the
#'  plot, only visible when using plotly (default is FALSE).
#' @importFrom dplyr filter arrange mutate %>%
#' @importFrom ggplot2 ggplot aes scale_color_manual ylab theme_minimal
#'   xlab theme element_text geom_smooth
#' @importFrom methods is
#' @importFrom stats reorder
#' @export
#' @examples
#' # Read data
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Do Analysis
#' exp <- doAnalysis(exp, doAll = TRUE)
#'
#' # Create a compound plot for the first compound
#' compoundPlot(
#'    exp,
#'    compound = 1,
#'    types = c("SQC", "LQC"),
#'    assay = "ratio_corrected",
#'    logTransform = TRUE
#' )
compoundPlot <- function(
        exp, compound = 1, batches = exp$batch, types = exp$type,
        assay = "ratio_corrected", logTransform = TRUE,
        withinTrend = FALSE, shortName = FALSE,
        trendTypes = metadata(exp)$QC,
        addInternalStandards = FALSE, addText = FALSE
) {

    df <- .preparePlotData(
        exp = exp[compound, ],
        assay = assay,
        batches = batches,
        types = types,
        logTransform = FALSE
    )


    if (addInternalStandards & metadata(exp)$hasIS) {
        df_is <- as.data.frame(.preparePlotData(exp[compound, ],
            assay = "area_is",
            batches = batches, types = types,
            logTransform = FALSE
        ))

        df_is$compound <- rowData(exp[compound, ])$compound_is
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
        ))  %>%
        arrange(.data$datetime) %>%
        mutate(aliquot = reorder(.data$aliquot, .data$injection_time))

    if (addText) {
        df <- .addTextHover(exp, df, type = "aliquot")
    }

    p <- scatterPlot(df, assay, trendTypes, logTransform, addText)
    return(p)
}

#' @title Compound plot function
#' @description Creates a scatter plot with optional trend lines and log
#'     scaling.
#' @details This function generates a scatter plot using ggplot2. It allows
#'     for optional trend lines based on specified types and can apply log
#'     scaling to the y-axis.
#' @returns A ggplot object representing the scatter plot.
#' @param df A data frame containing the data to be plotted. Must include
#'     columns for `aliquot`, the assay, `text`, `color`, `type`, and
#'     `Compound`.
#' @param assay A character string specifying the column name of the assay
#'     to be plotted on the y-axis.
#' @param trendTypes A character vector specifying the types of data to be
#'     used for trend lines.
#' @param logTransform A logical value indicating whether to apply log scaling
#'     to the y-axis. Defaults to TRUE.
#' @param addText Logical. Whether to add text bubbles to the
#'  plot, only visible when using plotly. Defaults to FALSE.
#' @importFrom dplyr filter arrange mutate
#' @importFrom ggplot2 ggplot aes scale_color_manual ylab theme_minimal xlab
#'     theme element_text geom_smooth .data ggtitle
#' @importFrom stats lm
#' @noRd
scatterPlot <- function(
        df, assay, trendTypes, logTransform = TRUE, addText = FALSE
) {

    mapping <- aes(x = .data$aliquot, y = .data[[assay]], fill = .data$type)
    if (addText) {
        mapping <- aes(
            x = .data$aliquot, y = .data[[assay]],
            text = .data$text, fill = .data$type
        )
    }

    idx <- !duplicated(df$type)
    colors <- setNames(df$color[idx], df$type[idx])

    p <- ggplot(df, mapping) +
        geom_point(
            aes(fill = .data$type), color = "black", stroke = 0.2, size = 3,
            shape = 21, na.rm = TRUE
        ) +
        scale_fill_manual(values = colors) +
        ylab(assay) +
        xlab("") +
        theme_minimal() +
        ggtitle(df$Compound[1]) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7))

    trendData <- filter(df, .data$type %in% trendTypes)
    if (nrow(trendData) > 0) {
        p <- p +
            geom_smooth(
                formula = "y ~ x", inherit.aes = TRUE, method = lm,
                show.legend = FALSE, na.rm = TRUE,
                aes(group = .data$group, color = .data$type), level = .95,
                data = trendData
            ) +
            scale_color_manual(values = colors)
    }

    if (logTransform) {
        p <- p + scale_y_log10()
    }

    return(p)
}

#' @title Aliquot Plot
#' @description Creates a boxplot for aliquots with optional log scaling.
#' @details This function generates a boxplot for aliquots using ggplot2. It
#'     allows for optional log scaling of the y-axis and customizes the plot
#'     based on the provided experimental data.
#' @returns A ggplot object representing the aliquot plot.
#' @param exp A data frame containing experimental data. Must include columns
#'     for `batch`, `type`, `Aliquot`, and `injection_time`.
#' @param assay A character string specifying the column name of the assay
#'     to be plotted on the y-axis.
#' @param batches A vector specifying the batches to include in the plot.
#'     Defaults to `exp$batch`.
#' @param types A vector specifying the types to include in the plot.
#'     Defaults to `exp$type`.
#' @param logTransform A logical value indicating whether to apply log scaling
#'     to the y-axis. Defaults to FALSE.
#' @importFrom dplyr filter arrange mutate
#' @importFrom ggplot2 ggplot aes scale_color_manual ylab theme_minimal xlab
#'     theme element_text scale_y_log10 .data
#' @export
#' @examples
#' # Read data
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Do Analysis
#' exp <- doAnalysis(exp, doAll = TRUE)
#'
#' # Create an aliquot plot
#' aliquotPlot(
#'    exp,
#'    assay = "ratio_corrected",
#'    types = c("SQC", "LQC"),
#'    logTransform = TRUE
#' )
aliquotPlot <- function(
        exp, assay, batches = exp$batch,
        types = exp$type, logTransform = FALSE
) {
    df <- .preparePlotData(
        exp = exp,
        assay = assay,
        batches = batches,
        types = types,
        logTransform = logTransform
    )

    mapping <- aes(
        x = reorder(.data$aliquot, .data$injection_time),
        y = .data[[assay]],
        text = .data$sample,
        fill = .data$type
    )

    idx <- !duplicated(df$type)
    colors <- setNames(df$color[idx], df$type[idx])

    p <- ggplot(df, mapping) +
        geom_boxplot(color = "black", na.rm = TRUE) +
        scale_fill_manual(values = colors) +
        ylab(assay) +
        theme_minimal() +
        xlab("") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7))

    if (!logTransform) {
        p <- p + scale_y_log10()
    }
    p
}

#' @title Adding confidence lines
#' @description Adds confidence interval lines to a ggplot object.
#' @details This function calculates the mean and standard deviation of
#'     medians for each aliquot and adds dashed horizontal lines to the plot
#'     at ±3 standard deviations from the mean.
#' @returns A ggplot object with added confidence interval lines.
#' @param plot A ggplot object to which confidence interval lines will be
#'     added.
#' @param df A data frame containing the data used to calculate confidence
#'     intervals. Must include columns for `Aliquot` and the assay.
#' @param assay A character string specifying the column name of the assay
#'     to be used for calculating confidence intervals.
#' @importFrom dplyr group_by summarise pull
#' @importFrom ggplot2 geom_hline .data
#' @importFrom stats sd
.addPlotConfidenceInterval <- function(plot, df, assay){
    medians <- df %>%
        group_by(.data$aliquot) %>%
        summarise(median = median(.data[[assay]], na.rm = TRUE)) %>%
        pull(.data$median)

    m <- mean(medians, na.rm = TRUE)
    sd <- sd(medians, na.rm = TRUE)

    plot <- plot +
        geom_hline(
            yintercept = m + sd * 3,
            linetype = "dashed",
            show.legend = FALSE
        ) +
        geom_hline(
            yintercept = m - sd * 3,
            linetype = "dashed",
            show.legend = FALSE
        )

    return(plot)
}

#' @title Violin Plot
#' @description Creates a violin plot with optional median points and
#'     confidence intervals.
#' @details This function generates a violin plot using ggplot2. It allows
#'     for optional log scaling, adding median points, and confidence
#'     intervals.
#' @returns A ggplot object representing the violin plot.
#' @param exp An experimental data object containing the data to be plotted.
#' @param assay A character string specifying the assay column to be plotted
#'     on the y-axis. Defaults to "ratio_corrected".
#' @param batches A vector specifying the batches to include in the plot.
#'     Defaults to `exp$batch`.
#' @param types A vector specifying the types to include in the plot.
#'     Defaults to `metadata(exp)$QC`.
#' @param logTransform A logical value indicating whether to apply log scaling
#'     to the y-axis. Defaults to FALSE.
#' @param addMedian A logical value indicating whether to add median points
#'     to the plot. Defaults to TRUE.
#' @param addConfidenceInterval A logical value indicating whether to add
#'     confidence intervals to the plot. Defaults to TRUE.
#' @param withinTrend A logical value indicating whether to inherit trend
#'     line aesthetics. Defaults to FALSE.
#' @importFrom dplyr filter arrange mutate
#' @importFrom ggplot2 ggplot aes scale_color_manual ylab theme_minimal xlab
#'     theme element_text scale_y_log10 stat_summary geom_violin .data
#' @importFrom stats sd
#' @importFrom dplyr group_by ungroup summarise pull n
#' @importFrom viridis viridis
#' @export
#' @examples
#' # Read data
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Do Analysis
#' exp <- doAnalysis(exp, doAll = TRUE)
#'
#' # Create a compound plot for the first compound
#' violinPlot(
#'    exp,
#'    assay = "ratio_corrected",
#'    logTransform = TRUE,
#'    addMedian = TRUE,
#'    addConfidenceInterval = TRUE
#' )
violinPlot <- function(
        exp, assay = "ratio_corrected",
        batches = exp$batch, types = metadata(exp)$QC,
        logTransform = FALSE, addMedian = TRUE,
        addConfidenceInterval = TRUE, withinTrend = FALSE
) {
    df <- .preparePlotData(
        exp, assay, batches, types = types, logTransform = logTransform
    )

    df <- .addTextHover(exp, df, type = "compound")

    mapping <- aes(x = .data$aliquot, y = .data[[assay]], fill = .data$batch)

    p <- ggplot(df, mapping = mapping) +
        geom_violin(na.rm = TRUE) +
        geom_point(
            aes(text = .data$text), size = .2,
            position = position_jitter(0.2), show.legend = FALSE
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        xlab("Aliquot")

    if (!logTransform) {
        p <- p + scale_y_log10()
    }

    if (addMedian) {
        p <- p + stat_summary(
            fun = "median", geom = "point", shape = 24, size = 3,
            fill = "black", show.legend = FALSE, color = "white", stroke = 0.5
        )
    }

    if (addConfidenceInterval) {
        p <- .addPlotConfidenceInterval(p, df, assay)
    }
    return(p)
}

#' @title Transform Data into PCA Data
#' @description Prepares data for PCA analysis by transforming and imputing
#'     missing or infinite values.
#' @details This function processes experimental data for PCA analysis. It
#'     imputes missing or infinite values, converts data to a wide format,
#'     and performs PCA using prcomp.
#' @returns A prcomp object containing the PCA results.
#' @param exp An experimental data object containing the data to be
#'     transformed.
#' @param assay A character string specifying the assay column to be used
#'     for PCA. Defaults to the first assay in `assayNames(exp)`.
#' @param batches A vector specifying the batches to include in the PCA.
#'     Defaults to `exp$batch`.
#' @param types A vector specifying the types to include in the PCA.
#'     Defaults to `exp$type`.
#' @param logTransform A logical value indicating whether to apply log scaling
#'     to the data. Defaults to TRUE.
#' @importFrom stats prcomp var xtabs
#' @importFrom dplyr group_by mutate ungroup if_else .data
#' @noRd
.getPcaData <- function(
        exp, assay = assayNames(exp)[1],
        batches = exp$batch, types = exp$type, logTransform = TRUE
) {
    df <- .preparePlotData(
        exp = exp, assay = assay, batches = batches,
        types = types, logTransform = logTransform
    )

    df <- df[is.finite(df[[assay]]), ]

    pc <- df %>%
        group_by(.data$aliquot) %>%
        # Impute infinite values
        mutate(value = if_else(
            condition = is.infinite(.data[[assay]]) | is.na(.data[[assay]]),
            true = min(.data[[assay]][is.finite(.data[[assay]])]) * 0.01,
            false = .data[[assay]]
        )) %>%
        ungroup() %>%
        mutate(aliquot = factor(.data$aliquot))

    # Convert long to wide
    pc <- xtabs(value ~ aliquot + compound, data = pc, subset = NULL) %>%
        as.data.frame.matrix()

    pc <- pc[, which(apply(pc, 2, var) != 0)]
    pc <- prcomp(pc, scale = TRUE, center = TRUE)

    return(pc)
}

#' @title Create a PCA plot with optional confidence intervals
#' @description Generates a PCA plot using ggplot2, with options to add
#'     confidence intervals for each group.
#' @details This function creates a PCA plot from the provided data, with
#'     points colored by group and optional confidence intervals (ellipses)
#'     for each group. The x-axis and y-axis labels can be customized, and
#'     colors for the groups are specified using a color vector.
#' @param data A data frame containing the PCA results. It must include
#'     columns `PC1`, `PC2`, `type`, and `text`.
#' @param xAxis A string specifying the label for the x-axis.
#' @param yAxis A string specifying the label for the y-axis.
#' @param addConfidenceInterval Logical, whether to add confidence intervals
#'     (ellipses) for each group. Defaults to `TRUE`.
#' @param addLabels Logical, whether to add labels to the points in the plot.
#'     requires `ggrepel` package. Defaults to `FALSE`.
#' @return A ggplot object representing the PCA plot.
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal xlab ylab
#' scale_color_manual scale_fill_manual guides stat_ellipse .data
#' @return A ggplot object representing the PCA plot.
#' @noRd
.createPcaPlot <- function(
        data, xAxis, yAxis,
        addConfidenceInterval = TRUE, addLabels = FALSE
) {
    pointMapping <- aes(
        x = .data$PC1, y = .data$PC2, fill = .data$type,
        group = .data$type, text = .data$text
    )

    intervalMapping <- aes(
        x = .data$PC1, y = .data$PC2,
        group = .data$type, color = .data$type
    )

    colors <- setNames(data$color, data$type)

    p <- ggplot(data = data, mapping = pointMapping) +
        geom_point(size = 3, stroke = 0.2, shape = 21, color = "black") +
        scale_fill_manual(values = colors) +
        theme_minimal() +
        xlab(xAxis) +
        ylab(yAxis) +
        guides(color = "legend")

    if (addConfidenceInterval) {
        p <- p + stat_ellipse(
            data = data, mapping = intervalMapping, inherit.aes = TRUE,
            na.rm = TRUE, level = 0.9
        ) +
            scale_color_manual(values = colors)

    }

    if (addLabels) {
        p <- .addPcaLabels(p)
    }

    return(p)
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
#' @param logTransform Logical flag indicating whether to log-transform the
#'   assay data before PCA (default is TRUE).
#' @param sampleAsBatch Logical flag. If TRUE, sample points are labeled
#'   using batch numbers instead of type (default is TRUE).
#' @param addConfidenceInterval Logical flag. If TRUE, a 95% confidence
#'   ellipse is drawn around each group (default is TRUE).
#' @param addLabels Logical, whether to add labels to the points in the plot.
#'     requires `ggrepel` package. Defaults to `FALSE`.
#' @importFrom ggplot2 guides stat_ellipse .data
#' @importFrom dplyr select bind_cols arrange all_of mutate filter .data
#' @importFrom viridis viridis
#' @export
#' @examples
#' # Read data
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Do Analysis
#' exp <- doAnalysis(exp, doAll = TRUE)
#'
#' # Create a compound plot for the first compound
#' pcaPlot(
#'    exp,
#'    assay = "ratio_corrected",
#'    logTransform = TRUE,
#'    types = c("SQC", "LQC", "SAMPLE")
#' )
pcaPlot <- function(
        exp, assay = assayNames(exp)[1],
        batches = exp$batch, types = exp$type,
        pc1 = 1, pc2 = 2, logTransform = TRUE,
        sampleAsBatch = TRUE, addConfidenceInterval = TRUE,
        addLabels = FALSE
) {

    pc <- .getPcaData(exp, assay, batches, types, logTransform)

    expvar <- (pc$sdev)^2 / sum(pc$sdev^2)
    xAxis <- sprintf("PC%s (%s%%)", pc1, format(expvar[pc1] * 100, digits = 4))
    yAxis <- sprintf("PC%s (%s%%)", pc2, format(expvar[pc2] * 100, digits = 4))

    exp <- exp[, exp$type %in% types & exp$batch %in% batches][, rownames(pc$x)]

    data <- pc$x %>%
        as.data.frame() %>%
        select(all_of(c(pc1, pc2))) %>%
        bind_cols(as.data.frame(colData(exp))) %>%
        arrange(.data$batch)

    data$aliquot <- rownames(data)

    if (sampleAsBatch) {
        idx <- data$type == "SAMPLE"
        data$type[idx] <- data$batch[idx]
    }

    types <- unique(data$type)
    colorVec <- viridis(length(types), begin = 0.2, end = 0.8, option = "H")
    names(colorVec) <- types

    colorVec <- coalesce(exp$color[types], colorVec)
    data$color <- colorVec[data$type]

    data <- .addTextHover(exp, data, type = "aliquot")

    .createPcaPlot(
        data = data, xAxis = xAxis, yAxis = yAxis,
        addConfidenceInterval = addConfidenceInterval
    )
}

#' @title Get Concentration Model Data
#' @description Prepares model data for concentration plots by filtering and
#'     optionally removing outliers.
#' @details This function filters the input data frame for a specific
#'     calibration type and optionally removes outliers based on the
#'     concentration outliers in the experimental data.
#' @returns A data frame containing the filtered and processed model data.
#' @param df A data frame containing the data to be processed.
#' @param exp An experimental data object containing metadata and outlier
#'     information.
#' @param calType A character string specifying the calibration type to
#'     filter the data.
#' @param removeOutliers A logical value indicating whether to remove
#'     outliers from the model data. Defaults to TRUE.
#' @importFrom dplyr %>% filter .data
#' @noRd
.getConcentrationModelData <- function(df, exp, calType, removeOutliers) {
    modelData <- df %>%
        filter(.data$type == calType)

    if (removeOutliers) {

        if (modelData$outliers[1] == 1) {
            modelData <- modelData[-1, ]
        }

        if (modelData$outliers[nrow(modelData)] == 1) {
            modelData <- modelData[-nrow(modelData), ]
        }
    }

    return(modelData)
}

#' @title Prepare Data for Concentration Plot
#' @description Prepares the data required for generating a concentration plot,
#'     including calibration and sample data, with optional adjustments for
#'     plotting on the calibration line.
#' @details This function processes the experimental data to extract and
#'     structure the necessary information for creating a concentration plot.
#'     It supports filtering by batch and types, and optionally adjusts
#'     concentration values for specific types when not plotting on the
#'     calibration line. Additional annotations, such as text bubbles, are
#'     added to the data for enhanced visualization.
#' @returns A dataframe containing the processed data for the
#'     concentration plot.
#' @param exp A SummarizedExperiment object containing the experimental data.
#' @param assay A character string specifying the assay to be used for
#'     plotting. Defaults to "concentration".
#' @param batch An integer specifying the batch index to filter the data.
#' @param types A character vector specifying the types of data to include
#'     in the plot (e.g., "SAMPLE", "SQC").
#' @param calType A character string specifying the calibration type to include
#'     in the plot. Defaults to the concentration metadata in `exp`.
#' @param plotOnCalibrationLine A logical value indicating whether to plot
#'     data on the calibration line. Defaults to TRUE.
#' @importFrom dplyr select all_of
#' @noRd
.getConcentrationPlotData <- function(
        exp, assay, batch, calType, types, plotOnCalibrationLine = TRUE
){

    df <- .preparePlotData(
        exp = exp,
        assay = sprintf("%s_concentration", calType),
        batches = batch,
        types = c(calType, types),
        logTransform = FALSE
    )

    df$concentration <- df[, sprintf("%s_concentration", calType), drop = TRUE]

    cols <- exp$type %in% c(calType, types) & exp$batch %in% batch

    df$assay <- c(unlist(t(assay(exp[, cols], assay))))

    outlierAssay <- sprintf("%s_Outliers", calType)
    df$outliers <- c(unlist(t(assay(exp[, cols], outlierAssay))))

    df$outliers[is.na(df$outliers)] <- 0
    df$batch <- as.factor(df$batch)

    if (!plotOnCalibrationLine) {
        df$concentration[df$type %in% types] <- 0
    }

    df <- .addTextHover(exp, df, type = "aliquot")
    return(df)
}

#' @title Concentration Plot
#' @description Creates a concentration plot with optional calibration lines
#'     and trend lines.
#' @details This function generates a concentration plot using ggplot2. It
#'     allows for customization of calibration types, outlier removal, and
#'     whether to plot on the calibration line.
#' @returns A ggplot object representing the concentration plot.
#' @param exp An experimental data object containing the data to be plotted.
#' @param assay A character string specifying the assay column to be plotted
#'     on the y-axis. Defaults to "ratio_corrected".
#' @param compound An integer specifying the compound index to be plotted.
#'     Defaults to 1.
#' @param batch An integer specifying the batch index to be plotted.
#'     Defaults to 1.
#' @param types A character vector specifying the types of data to include
#'     in the plot. Defaults to `c("SAMPLE", "SQC")`.
#' @param calType A character string specifying the calibration type.
#'     Defaults to the concentration metadata in `exp`.
#' @param removeOutliers A logical value indicating whether to remove
#'     outliers from the data. Defaults to TRUE.
#' @param withinTrend A logical value indicating whether to inherit trend
#'     line aesthetics. Defaults to TRUE.
#' @param plotOnCalibrationLine A logical value indicating whether to plot
#'     data on the calibration line. Defaults to TRUE.
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_smooth
#'     theme_minimal
#' @export
#' @examples
#' # Read data
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Do Analysis
#' exp <- doAnalysis(exp, doAll = TRUE)
#'
#' # Create a compound plot for the first compound
#' concentrationPlot(
#'    exp,
#'    assay = "ratio_corrected",
#'    types = c("SQC", "LQC", "SAMPLE")
#' )
concentrationPlot <- function(
        exp, assay = "ratio_corrected", compound = 1,
        batch = exp$batch[1], types = c("SAMPLE", "SQC"),
        calType = metadata(exp)$concentration,
        removeOutliers = TRUE, withinTrend = TRUE,
        plotOnCalibrationLine = TRUE
) {
    calType <- toupper(calType)
    types <- toupper(types)
    exp <- exp[rowData(exp)$hasKnownConcentrations, ]

    df <- .getConcentrationPlotData(
        exp = exp[compound, ],
        assay = assay,
        batch = batch,
        calType = calType,
        types = types,
        plotOnCalibrationLine = plotOnCalibrationLine
    )

    modelData <- .getConcentrationModelData(
        df = df,
        exp = exp[compound, ],
        calType = calType,
        removeOutliers = removeOutliers
    )

    ggplot(df, aes(
        x = .data$concentration, y = .data$assay, fill = .data$type,
        color = .data$type, group = .data$batch, text = .data$text
        )) +
        geom_point(pch = 21, stroke = 0.1, size = 3, color = "black") +
        geom_hline(
            mapping = aes(yintercept = max(.data$assay)),
            linetype = "dashed", data = modelData
        ) +
        geom_hline(
            mapping = aes(yintercept = min(.data$assay)),
            linetype = "dashed", data = modelData
        ) +
        geom_smooth(
            formula = "y ~ x", method = lm, inherit.aes = withinTrend,
            aes(x = .data$concentration, y = .data$assay),
            na.rm = TRUE, show.legend = FALSE,
            data = modelData, level = .95
        ) +
        theme_minimal() +
        ylab(assay)
}

#' @title Plot heatmap of RSDQCs with internal standards
#' @description This function plots a heatmap of batch-corrected RSDQCs
#' of all compound - internal standard combinations. This allows for insight
#' if the picked internal standard is the correct one or if possible another
#' internal standard should be picked.
#' @param exp SummarizedExperiment object
#' @returns Plotly heatmap of compound - internal standard RSDQCs
#' @importFrom heatmaply heatmaply
#' @export
#' @examples
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Plot RSDQCs
#' rsdqcPlot(exp)
rsdqcPlot <- function(exp) {
    if (!isValidExperiment(exp)) {
        stop("Invalid experiment")
    }

    if (!metadata(exp)$hasIS) {
        return(NULL)
    }

    exp <- exp[!is.na(rowData(exp)$rsdqcCorrected), ]
    rsdqcs <- matrixRSDQCs(exp)

    cols <- viridis::viridis(
        n = 256, alpha = 1, begin = 0.2, end = 0.8,
        option = "H"
    )
    p <- heatmaply(as.matrix(rsdqcs),
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
#' @importFrom heatmaply heatmaply
#' @importFrom viridis viridis
#' @export
#' @examples
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Plot values in heatmap
#' heatmapPlot(exp)
heatmapPlot <- function(exp, assay = "ratio") {
    if (!isValidExperiment(exp)) {
        stop("Invalid experiment")
    }
    df <- assay(exp, assay)
    df[is.infinite(df)] <- NA
    if (all(is.na(df))) {
        return(NULL)
    }

    cols <- viridis(
        n = 256, alpha = 1, begin = 0.2, end = 0.8,
        option = "H"
    )

    p <- heatmaply(
        df,
        dendrogram = "none", plot_method = "ggplot",
        margins = c(0, 0, 0, 0), colors = cols
    )

    return(p)
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
#' @importFrom ggplot2 geom_line theme geom_point scale_color_manual
#' scale_fill_manual theme_bw
#' @examples
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Create scatterplot with Relative Standard Deviation values
#' rsdPlot(exp, number = 10)
rsdPlot <- function(
        exp, assay = "ratio_corrected",
        qc = "SQC", number = nrow(exp)
) {
    if (!isValidExperiment(exp)) {
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
    df <- data.frame(
        compound = rep(rownames(df), times = ncol(df)),
        batch = rep(colnames(df), each = nrow(df)),
        RSD = as.vector(as.matrix(df))
    )

    df$batch <- factor(df$batch, levels = unique(df$batch))


    typeColors <- setNames(viridis::viridis(length(unique(df$compound)),
        option = "H"
    ), unique(df$compound))

    df$text <- .getCompoundPlotLabels(exp)

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
.getCompoundPlotLabels <- function(exp) {
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
#' theme_bw theme
#' @export
#' @examples
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Create Batch Assay Plot
#' batchAssayPlot(exp)
batchAssayPlot <- function(exp, assay = 1, compound = 1) {
    if (!isValidExperiment(exp)) {
        stop("Invalid experiment")
    }
    exp <- exp[compound, ]

    m <- assay(exp, assay)
    df <- data.frame(
        aliquot = rep(rownames(m), times = ncol(m)),
        compound = rep(colnames(m), each = nrow(m)),
        value = as.vector(m)
    )

    df <- df[order(df$compound, df$aliquot), ]

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
