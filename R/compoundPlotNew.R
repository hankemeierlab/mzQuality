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

#' @title New Compound Plot
#' @description
#' A short description...
#' @details placeholder
#' @returns placeholder
#' @param exp
#' @param compound
#' @param batches
#' @param types
#' @param assay
#' @param doLog
#' @param withinTrend
#' @param shortName
#' @param trendTypes
#' @importFrom dplyr filter arrange mutate %>%
#' @importFrom ggplot2 ggplot aes scale_color_manual ylab theme_minimal xlab
#' theme element_text geom_smooth
#' @importFrom methods is
#' @importFrom stats reorder
#' @export
compoundPlotNew <- function(exp, compound = 1, batches = exp$batch, types = exp$type,
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


    p <- compoundPlot(df, assay, trendTypes, colors)
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
compoundPlot <- function(df, assay, trendTypes, colors, smooth = TRUE, doLog = TRUE) {
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

#' @title New Aliquot Plot
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
aliquotPlotNew <- function(exp, assay, batches = exp$batch, types = exp$type,
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

#' @title New Violin Plot
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
#' @export
violinPlotNew <- function(exp, assay = "ratio_corrected",
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
        mutate(Var1 = factor(.data$Var1)) %>%
        # Convert long to wide
        xtabs(value ~ Var1 + Var2, data = ., subset = NULL) %>%
        as.data.frame.matrix()



    pc <- pc[, which(apply(pc, 2, var) != 0)]
    pc <- prcomp(pc, scale = TRUE, center = TRUE)

    return(pc)
}

pcaLoadingPlot <- function(pc) {
    loadings <- data.frame(pc$rotation, Compound = rownames(pc$rotation))

    ggplot(loadings, aes(x = PC1, y = PC2, label = Compound)) +
        geom_point() +
        geom_text(vjust = -0.5, hjust = 1.1) +
        # xlim(c(-1, 1)) +
        # ylim(c(-1, 1)) +
        theme_minimal()
}


pcaBiPlot <- function(pc) {
    scores <- as.data.frame(pc$x)
    loadings <- data.frame(pc$rotation, variable = rownames(pc$rotation))

    scale_factor <- max(abs(scores[, c(1, 2)])) / max(abs(loadings[, c(1, 2)]))
    loadings[, c(1, 2)] <- loadings[, c(1, 2)] * scale_factor

    ggplot() +
        geom_point(data = scores, aes(x = PC1, y = PC2)) +
        # geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2),
        #             arrow = arrow(length = unit(0.3, "cm")), color = "red") +
        geom_text(data = loadings, aes(x = PC1, y = PC2, label = variable))
}



#' @title New PCA Plot
#' @description placeholder
#' @details placeholder
#' @returns placeholder
#' @param exp
#' @param assay
#' @param batches
#' @param types
#' @param pc1
#' @param pc2
#' @param doLog
#' @param sampleAsBatch
#' @param addConfidenceInterval
#' @param typeColors
#' @importFrom ggplot2 guides stat_ellipse
#' @export
pcaPlotNew <- function(exp, assay = assayNames(exp)[1],
                       batches = exp$batch, types = exp$type,
                       pc1 = 1, pc2 = 2,
                       doLog = TRUE, sampleAsBatch = TRUE,
                       addConfidenceInterval = TRUE,
                       typeColors = c(SQC = "black", LQC = "grey")) {
    # PCA
    pc <- getPcaData(exp, assay, batches, types, doLog)


    # Create axis labels
    expvar <- (pc$sdev)^2 / sum(pc$sdev^2)

    pc1name <- sprintf("PC%s (%s%%)", pc1, format(expvar[pc1] * 100, digits = 4))
    pc2name <- sprintf("PC%s (%s%%)", pc2, format(expvar[pc2] * 100, digits = 4))


    exp <- exp[, exp$type %in% types & exp$batch %in% batches]
    exp <- exp[, rownames(pc$x)]

    data <- pc$x %>%
        as.data.frame() %>%
        # Select the PCs chosen
        select(all_of(c(pc1, pc2))) %>%
        # Bind column metadata
        bind_cols(as.data.frame(colData(exp))) %>%
        mutate(batch = formatC(as.integer(.data$batch), max(nchar(.data$batch)), format = "d", flag = "0")) %>%
        arrange(.data$batch) %>%
        mutate(type = ifelse(.data$type == "SAMPLE" & sampleAsBatch, .data$batch, .data$type))

    data$aliquot <- rownames(data)
    types <- sort(unique(data$type))
    colorVec <- viridis::viridis(length(types), begin = 0.2, end = 0.8, option = "H")

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
    p <- ggplot(data = data, aes(
        x = .data$PC1, y = .data$PC2, fill = .data$type, group = .data$type, text = .data$text)) +
        geom_point(size = 3, stroke = 0.2, shape = 21, color = "black") +
        theme_minimal() +
        xlab(pc1name) +
        ylab(pc2name) +
        scale_color_manual(values = colorVec) +
        scale_fill_manual(values = colorVec) +
        guides(color = "legend")

    if (addConfidenceInterval) {
        p <- p + stat_ellipse(
            data = data, level = 0.95, mapping = aes(group = .data$type, color = .data$type, x = .data$PC1, y = .data$PC2),
            inherit.aes = TRUE, na.rm = TRUE
        )
    }


    return(p)
}

#' @title New Concentration Plot
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
concentrationPlotNew <- function(exp, assay = "ratio_corrected", compound = 1,
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
        geom_hline(mapping = aes(yintercept = max(ratio)), linetype = "dashed", data = modelData) +
        geom_hline(mapping = aes(yintercept = min(ratio)), linetype = "dashed", data = modelData) +
        geom_smooth(
            formula = "y ~ x", method = lm, inherit.aes = withinTrend, level = .95,
            aes(x = .data$concentration, y = .data$ratio),
            na.rm = TRUE, show.legend = FALSE,
            data = modelData
        ) +
        theme_minimal()
}
