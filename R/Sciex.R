#' @title Parse Sciex Vendor Files
#' @description Processes Sciex vendor files to extract and standardize
#'     aliquot and compound data.
#' @param df A dataframe containing Sciex vendor data.
#' @param regex A regular expression for parsing aliquot identifiers.
#' @returns A standardized dataframe with parsed aliquot and compound data.
#' @noRd
.processSciex <- function(df, regex = NULL) {
    samples <- as.vector(unlist(df[, 1]))
    df <- .whitespaceFix(df)

    aliquot_df <- .parseAliquots(samples, regex)

    IScolumns <- c(
        "IS Name", "IS Retention Time", "IS Area",
        "IS.Name", "IS.Area", "IS.Retention.Time"
    )
    if (!any(IScolumns %in% colnames(df))) {
        df$`IS Area` <- 1
        df$`IS Retention Time` <- 1
        df$`IS Name` <- "NA"
    }

    mandatoryColumns <- c(
        "Acquisition Date & Time", "Acquisition.Date...Time",
        "Component Name", "Component.Name",
        "Retention Time", "Retention.Time",
        "Area",
        "IS Name", "IS.Name",
        "IS Retention Time", "IS.Retention.Time",
        "IS Area", "IS.Area"
    )

    mandatoryData <- do.call(cbind, lapply(mandatoryColumns, function(x) {
        df[, which(colnames(df) == x)]
    }))

    colnames(mandatoryData) <- c(
        "datetime", "compound", "rt", "area",
        "compound_is", "rt_is", "area_is"
    )

    res <- cbind(aliquot_df, mandatoryData)

    otherColumns <- setdiff(colnames(df), mandatoryColumns)
    if (length(otherColumns) > 0) {
        res <- cbind(res, df[, otherColumns, drop = FALSE])
    }

    as.data.frame(res)
}

#' @title Fix Whitespace Issues in Sciex Files
#' @description Resolves whitespace and duplicate line issues in Sciex
#'     vendor files.
#' @param f1 A dataframe with potential whitespace or duplicate issues.
#' @returns A cleaned dataframe with resolved whitespace issues.
#' @importFrom dplyr coalesce
#' @noRd
.whitespaceFix <- function(f1) {

    indexes <- which(trimws(f1[, 1], which = "right") == "")
    if (length(indexes) == 0) {
        return(f1)
    }

    all <- sort(c(indexes, indexes - 1))
    f1[f1 == ""] <- NA
    df <- do.call(rbind, lapply(indexes - 1, function(i) {
        dplyr::coalesce(f1[i, ], f1[i + 1, ])
    }))
    to_remove <- unique(f1[, 1][indexes - 1])
    f1 <- f1[-all, ]
    f1 <- f1[which(f1[, 1] != to_remove), ]
    f1 <- rbind(f1, df)
    f1
}

#' @title Parse Aliquot Names
#' @description Parses aliquot names using a regular expression to extract
#'     metadata such as sample type, calibration number, and replicate.
#' @param aliquots A vector of aliquot names.
#' @param regex A regular expression to parse the aliquot names.
#' @returns A dataframe with parsed aliquot metadata.
#' @noRd
.parseAliquots <- function(aliquots, regex = NULL) {
    if (is.null(regex)) {
        regex <- paste0(
            "([0-9]{4}[A-Z]{3}_",
            "[0-9]{4})([A-Z]+|)([0-9]+|)_([A-Z])([0-9])"
        )
    }

    aliquots <- aliquots[nchar(aliquots) > 0]
    aliquots <- toupper(aliquots)
    matches <- lapply(regmatches(aliquots, gregexec(regex, aliquots)), t)
    df <- as.data.frame(do.call(rbind, matches))
    colnames(df) <- c(
        "aliquot", "sample", "type", "calno", "replicate", "injection"
    )

    df$type[df$type == ""] <- "SAMPLE"
    df$type <- toupper(df$type)
    df$calno[df$calno == ""] <- NA
    df$calno[!is.na(df$calno)] <- as.integer(df$calno[!is.na(df$calno)])

    return(df)
}

#' @title Correct Batch Order by Datetime
#' @description Orders batches by datetime and assigns batch numbers.
#' @param list_of_files A list of dataframes representing batches.
#' @param batch_order An optional vector specifying the batch order.
#' @returns A dataframe with corrected batch order and assigned batch numbers.
#' @noRd
.setBatches <- function(list_of_files, batch_order = NULL) {
    list_of_files <- .fixDates(list_of_files)
    if (is.null(batch_order)) {
        batch_order <- vapply(list_of_files, function(file_df) {
            min(file_df$datetime)
        }, numeric(1))
    }

    list_of_files <- list_of_files[order(batch_order)]
    names(list_of_files) <- seq_len(length(list_of_files))
    df <- do.call(rbind, lapply(names(list_of_files), function(name) {
        list_of_files[[name]]$batch <- name
        list_of_files[[name]]
    }))
    df
}



#' @title Guess the Date Format
#' @description Determines the date format of a vector of datetime strings.
#' @param datetimes A character vector of datetime strings.
#' @returns A string representing the guessed date format.
#' @importFrom lubridate parse_date_time
#' @noRd
.interpretDates <- function(datetimes) {
    dates <- do.call(rbind, strsplit(datetimes, " "))[, 1]
    possible <- c("dmy", "mdy", "ymd")

    # Guess possible formats
    guesses <- lubridate::guess_formats(dates, possible)
    formats <- guesses[names(guesses) %in% possible]
    counts <- table(names(formats))

    # Only keep formats where all guesses could fit. Can result in multiple
    # formats
    formats <- formats[names(counts[counts == length(dates)])]

    # Check if only 1 option remained
    if (length(formats) == 1 | length(dates) == 1) {
        return(paste0(names(formats), "HMS"))
    }
    # If multiple options are possible, check if the difference in days
    # is less than 12 to distinguish between d/m/y and m/d/y
    formats <- formats[unlist(lapply(as.vector(formats), function(format) {
        time_1 <- lubridate::parse_date_time(unique(dates), format)
        all(abs(as.integer(diff(time_1))) <= 12)
    }))]

    # Only 1 format should remain
    if (length(unique(formats)) == 1) {
        return(paste0(names(formats), "HMS"))
    }

    # All formats have been measured on the same day, cannot determine
    # if it would be d/m/y or m/d/y, returning d/m/y as it is the most common
    return("dmyHMS")
}

#' @title Fix Dates in Batch Data
#' @description Converts datetime strings in batch data to proper date
#'     objects.
#' @param list_of_batches A list of dataframes containing batch data.
#' @returns A list of dataframes with corrected datetime formats.
#' @importFrom lubridate parse_date_time as_datetime
#' @noRd
.fixDates <- function(list_of_batches) {
    date_vec <- unlist(lapply(list_of_batches, function(df) {
        as.character(unique(df$datetime))
    }))
    date_format <- .interpretDates(date_vec)
    lapply(list_of_batches, function(x) .fixHMS(x, date_format))
}

#' @title Fix Datetime Format
#' @description Converts datetime strings in a dataframe to proper date
#'     objects, handling both H:M and H:M:S formats.
#' @param df A dataframe containing datetime strings.
#' @param date_format The date format to use for parsing.
#' @returns A dataframe with corrected datetime and injection time columns.
#' @importFrom lubridate parse_date_time as_datetime
#' @noRd
.fixHMS <- function(df, date_format) {
    # Check if the time is hours:minutes or hours:minutes:seconds
    time <- gregexpr(":", df$datetime)
    counts <- vapply(time, function(x) sum(x > 0), double(1))
    no_seconds <- which(counts == 1)
    if (length(no_seconds) > 0) {
        # No seconds, parse as H:M
        df$datetime[no_seconds] <- lubridate::parse_date_time(
            df$datetime[no_seconds],
            substr(date_format, 1, nchar(date_format) - 1)
        )
    }
    with_seconds <- which(counts == 2)
    if (length(with_seconds) > 0) {
        # With seconds, parse as H:M:S
        df$datetime[with_seconds] <- lubridate::parse_date_time(
            df$datetime[with_seconds],
            date_format
        )
    }
    df$datetime <- lubridate::as_datetime(as.double(df$datetime))
    df$injection_time <- as.integer(df$datetime)
    df
}
