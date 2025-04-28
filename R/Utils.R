#' @title Apply jobs multicore
#' @description
#' @details
#' @returns
#' @param func
#' @param jobs
#' @param cores
#' @param progressTitle
#' @param shinyProgress
#' @param verbose
#' @param ...
#' @importFrom foreach foreach %do% %dopar%
#' @export
multiApply <- function(func, jobs, cores = 1, progressTitle = "",
                       shinyProgress = NULL,
                       verbose = FALSE, ...) {
    if (jobs == 0) {
        return(NULL)
    }

    pb <- progressBar(
        title = progressTitle, total = jobs,
        start = 0, shinyProgress
    )



    # Conditionally switch between parallel and serial in a `foreach` loop.
    `%execute%` <- ifelse(cores > 1, `%dopar%`, `%do%`)

    # Create cluster for multi-processing
    cl <- createCluster(cores, seq_len(jobs))


    # Perform parallel MS1/MS2 correlating
    results <- foreach(
        i = seq_len(jobs),
        .verbose = verbose,
        .errorhandling = "pass",
        .options.snow = list(progress = function(n) {
            updateProgress(pb, jobs, n, shinyProgress)
        })
    ) %execute% {
        i <- get("i")
        res <- tryCatch(
            {
                func(i, ...)
            },
            error = function(x) warning(x)
        )

        if (cores == 1) updateProgress(pb, jobs, i, shinyProgress)

        return(res)
    }

    # Stop the clusters
    if (!is.null(cl)) parallel::stopCluster(cl)

    return(results)
}

#' @title Create cluster for multi-processing
#' @description This function aids in registering a cluster for
#' multi-processing.
#' @details Multi-processing is done by creating a SNOW cluster where each core
#' can handle a single task. This speeds up processing of multiple samples.
#' @returns Cluster object from the `parallel` package.
#' @param cores Integer with the amount of cores to use
#' @param files Character vector with each representing one task. This is used
#' to prevent that more cores are assigned than tasks available
#' @importFrom parallel makeCluster
#' @importFrom doSNOW registerDoSNOW
#' @noRd
createCluster <- function(cores, files) {
    cl <- NULL


    if (cores > 1) {
        cores <- ifelse(cores > length(files), length(files), cores)
    } else {
        cores <- 1
    }

    cl <- makeCluster(cores)
    registerDoSNOW(cl)

    return(cl)
}

#' @title Create progress bar for long during tasks
#' @description Progressbars are shown in the console to indicate the progress
#' of a particular task
#' @details Progressbars indicate the progress of a particular task. This
#' shows a title of the current progress, the bar itself, the current iteration,
#' and an estimated time for completion.
#' @returns A progress_bar object from the package `progress`
#' @param title Text shown on the left of the progress bar
#' @param total Total number of tasks to complete
#' @param start Starting number for tasks already completed
#' @param shinyProgress Optional identifier for a progressbar in a shiny application
#' Defaults to `NULL`
#' @noRd
progressBar <- function(title, total, start = 0, shinyProgress = NULL) {
    pb <- progress::progress_bar$new(
        format = sprintf("%s [:bar] | :current/:total (:percent) | Elapsed: :elapsedfull | Remaining: :eta", title),
        total = total, force = TRUE, show_after = 0, clear = FALSE
    )

    if (requireNamespace("shiny", quietly = TRUE) &
        requireNamespace("shinyWidgets", quietly = TRUE)) {
        if (shiny::isRunning()) {
            shinyWidgets::updateProgressBar(
                id = shinyProgress,
                value = 0,
                total = total,
                status = "primary",
                range_value = c(0, total)
            )
        } else {
            pb$tick(0)
        }
    } else {
        pb$tick(0)
    }

    pb
}

#' @title Update the progressbar
#' @description This function aids in updating the progressbar for long-running
#' tasks.
#' @details When a shiny application is run and `shinyProgress` is given, this
#' function will also update the progress bar in the shiny app.
#' @param pb Progressbar object
#' @param total Maximum value that the progressbar can have
#' @param n Current value of the progressbar
#' @param shinyProgress Optional identifier of a progressbar in shiny.
updateProgress <- function(pb, total, n, shinyProgress = NULL) {
    if (requireNamespace("shiny", quietly = TRUE)) {
        if (!is.null(shinyProgress) & shiny::isRunning()) {
            shinyWidgets::updateProgressBar(
                id = shinyProgress,
                value = n,
                total = total,
                status = "primary",
                range_value = c(0, total)
            )
        }
    }
    pb$tick(1)
}


sampleColors <- function(..., overwrite = FALSE) {
    l <- c(...)
    defaults <- c(
        SAMPLE = "red",
        SQC = "black",
        LQC = "gray",
        ISTD = "pink"
    )

    if (overwrite) {
        params <- c(l, defaults)
    } else {
        params <- c(defaults, l)
    }

    unlist(params[!duplicated(names(params))])
}
