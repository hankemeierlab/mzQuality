#' @title Clean a Dataframe by Removing Non-Numeric or Matrix Columns
#' @description Cleans a dataframe by removing columns of class "matrix" and,
#'     optionally, retaining only numeric columns.
#' @param df A dataframe to be cleaned.
#' @param onlyNumeric Logical. If TRUE, retains only numeric columns in the
#'     dataframe. Defaults to FALSE.
#' @returns A cleaned dataframe with the specified columns removed.
#' @export
#' @examples
#' # Example dataframe
#' df <- data.frame(
#'     A = matrix(1:4, nrow = 2),
#'     B = c(1.5, 2.5),
#'     C = c("x", "y")
#' )
#'
#' # Remove matrix columns
#' cleanDataframe(df)
#'
#' # Remove matrix columns and retain only numeric columns
#' cleanDataframe(df, onlyNumeric = TRUE)
cleanDataframe <- function(df, onlyNumeric = FALSE) {
    columns <- vapply(df, function(x) !("matrix" %in% is(x)), logical(1))

    df <- as.data.frame(df[, columns, drop = FALSE])
    if (onlyNumeric) {
        columns <- vapply(df, function(x) ("numeric" %in% is(x)), logical(1))
        as.data.frame(df[, columns, drop = FALSE])
    }

    return(df)
}

#' @title Define Sample Colors
#' @description Creates a named vector of colors for sample types, allowing
#'     customization and optional overwriting of default colors.
#' @param ... Named arguments specifying custom colors for sample types.
#' @param overwrite Logical. If TRUE, custom colors overwrite defaults.
#' @returns A named character vector of sample colors.
#' @examples
#' # Default colors
#' .sampleColors()
#'
#' # Custom colors without overwriting defaults
#' .sampleColors(SAMPLE = "blue", SQC = "green")
#'
#' # Custom colors with overwriting defaults
#' .sampleColors(SAMPLE = "blue", overwrite = TRUE)
#' @noRd
.sampleColors <- function(..., overwrite = FALSE) {
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
