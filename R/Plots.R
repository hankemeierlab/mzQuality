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
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality2"))
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
    rsdqcs <- calculateCorrectedRSDQCs2(exp, returnRSDs = TRUE)

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
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality2"))
#'
#' # Plot values in heatmap
#' heatmapPlot(exp)
heatmapPlot <- function(exp, assay = "ratio") {
    if (!validateExperiment(exp)) {
        stop("Invalid experiment")
    }
    df <- assay(exp, assay)
    df[is.infinite(df)] <- NA
    if (all(is.na(df))) return(NULL)

    cols <- viridis::viridis(
        n = 256, alpha = 1, begin = 0.2, end = 0.8,
        option = "H"
    )

    p <- heatmaply::heatmaply(
        df,
        dendrogram = "none",  plot_method = "ggplot",
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
#' @export
#' @examples
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality2"))
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
      return(ggplot() + theme_bw())
    }

    dl <- dl[, -which(colnames(dl) == "calno")]
    dl <- dl[complete.cases(dl), ]

    #points$text <- getAliquotPlotLabels(aliquotFilter(sub, type = c(qc, "SAMPLE")))


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
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality2"))
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
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality2"))
#'
#' # Create Batch Assay Plot
#' batchAssayPlot(exp)
batchAssayPlot <- function(exp, assay = 1, compound = 1){
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
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality2"))
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
