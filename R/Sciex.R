#' @title Parse Sciex vendor files
#' @param file Vector of file names
#' @param regex regular expression for parsing aliquot identifiers
#' @noRd
sciex <- function(df, regex = NULL) {
    samples <- as.vector(unlist(df[, 1]))
    if (sum(trimws(samples, which = "right") == "") > 0) {
        df <- whitespaceFix(df)
    }

    aliquot_df <- parseAliquots(samples, regex)

    if (!any(c("IS Retention Time", "IS Area", "IS.Name", "IS.Area", "IS.Retention.Time") %in% colnames(df))) {
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

    res <- cbind(data.frame(
        aliquot = aliquot_df$Aliquot,
        sample = aliquot_df$Sample,
        type = aliquot_df$Type,
        calNo = aliquot_df$CalNo,
        replicate = aliquot_df$Replicate,
        injection = aliquot_df$Injection
    ), mandatoryData)

    otherColumns <- colnames(df)[!colnames(df) %in% c(colnames(df)[1], mandatoryColumns)]
    if (length(otherColumns) > 0) {
        res <- cbind(res, df[, otherColumns, drop = FALSE])
    }

    as.data.frame(res)
}

#' @title Fix whitespaces issues in Sciex files.
#' @param f1 Dataframe with whitespaces and/or duplicate lines.
#' @importFrom dplyr coalesce
#' @noRd
whitespaceFix <- function(f1) {
    message("Found inconsistencies in your file. Trying to fix..")
    indexes <- which(trimws(f1[, 1], which = "right") == "")
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

#' @title Parse aliquot names with regex to obtain aliquot data
#' @param aliquots Vector of aliquot names
#' @param regex Regular expression to parse the names
#' @noRd
parseAliquots <- function(aliquots, regex = NULL) {
    if (is.null(regex)) {
        regex <- paste0(
            "([0-9]{4}[A-Z]{3}_",
            "[0-9]{4})([A-Z]+|)([0-9]+|)_([A-Z])([0-9])"
        )
    }

    aliquots <- aliquots[nchar(aliquots) > 0]
    aliquots <- toupper(aliquots)
    matches <- regmatches(aliquots, gregexec(regex, aliquots))
    conv <- lapply(matches, t)
    df <- as.data.frame(do.call(rbind, conv))
    colnames(df) <- c(
        "Aliquot", "Sample", "Type", "CalNo", "Replicate",
        "Injection"
    )

    df$Type[df$Type == ""] <- "SAMPLE"
    df$Type <- toupper(df$Type)
    df$CalNo[df$CalNo == ""] <- NA
    df$CalNo[!is.na(df$CalNo)] <- as.integer(df$CalNo[!is.na(df$CalNo)])


    df
}

#' @title Correct batch order by datetime
#' @param list_of_files List of dataframes of the files
#' @param batch_order Optional given order of batches
#' @noRd
setBatches <- function(list_of_files, batch_order = NULL) {
    list_of_files <- fixDates(list_of_files)
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



#' @title Guess the date format
#' @param datetimes Character vector of datetimes
#' @importFrom lubridate parse_date_time
#' @noRd
guessDates <- function(datetimes) {
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

#' @title Fix the dates to actual date type
#' @param list_of_batches List object with all batches
#' @importFrom lubridate parse_date_time as_datetime
#' @importFrom stringr str_count
#' @noRd
fixDates <- function(list_of_batches) {
    date_vec <- unlist(lapply(list_of_batches, function(df) {
        as.character(unique(df$datetime))
    }))
    date_format <- guessDates(date_vec)
    lapply(list_of_batches, function(x) fixHMS(x, date_format))
}

fixHMS <- function(df, date_format) {
    # Check if the time is hours:minutes or hours:minutes:seconds
    counts <- stringr::str_count(df$datetime, ":")
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
