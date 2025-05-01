vectorDf <- function(df) {
    df <- as.data.frame(df)
    columns <- vapply(df, function(x) !("matrix" %in% class(x)), logical(1))
    return(df[, columns, drop = FALSE])
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
