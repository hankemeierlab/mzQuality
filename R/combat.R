#' @title Combat batch correction
#' @description
#' A Bayes Theorem based batch correction
#' @details
#' Workflow:
#' 1. Standardize the data, so that the values (compounds) become relative
#' across samples. By standardizing, you're subtracting the mean and divide
#' by the standard deviation
#' 2. Using the batch labels, multiply the standardized score with the batch
#' occurence so that a matrix of batch x compound contains the weight of each
#' compound for a batch (gamma_hat).
#' 3. Calculate the mean of each batch, so we get the average values of each
#' batch with the standardized compounds (gamma_bar).
#' 4. Calculate the variance per batch (gamma_var)
#' 5. Multiply the batch variance (gamma_var) by the weights (gamma_hat) and
#' add the mean to get the batch-corrected weights. Then divide by the
#' variance + 1
#' 6. For each compound, multiply by the batch occurence. These are the standardized
#' batch effects. Subtract this from the original standardized values to obtain
#' batch-corrected standardized values.
#' 7. Reverse the standardization by multiplying and adding the original
#' standard deviation and means
#' 8. Return the batch-corrected values
#' @param data matrix of compounds x samples
#' @param batch Vector of batch labels, same length as columns in data
#' @export
#' @examples
#' # Simulate some data
#' set.seed(1)
#' n <- 100 # samples
#' p <- 200 #
#' x <- seq(10, 250, length.out = n)
#' batch <- sample(1:5, n, replace = TRUE)
#' data <- sapply(1:p, function(i) 2 * x + rnorm(n, 0, 10) + (batch - 1) * 20)
#' # Transpose data to have samples in columns and features in rows
#' data <- t(data)
#' corrected_data <- combat(data, batch)
#' # Perform PCA
#' pca_original <- prcomp(t(data), scale. = TRUE)
#' pca_corrected <- prcomp(t(corrected_data), scale. = TRUE)
#' # Plot PCA results for original data
#' par(mfrow = c(1, 2))
#' plot(pca_original$x[, 1], pca_original$x[, 2], col = batch, pch = 19,
#'      main = "PCA of Original Data", xlab = "PC1", ylab = "PC2")
#' # Plot PCA results for batch-corrected data
#' plot(pca_corrected$x[, 1], pca_corrected$x[, 2], col = batch, pch = 19,
#'      main = "PCA of Batch-Corrected Data", xlab = "PC1", ylab = "PC2")
combat <- function(data, batch) {

    batch <- as.factor(batch)
    design <- model.matrix(~ 0 + batch)

    # Standardize data
    means <- rowMeans(data, na.rm = TRUE)
    sds <- apply(data, 1, sd, na.rm = TRUE)
    data <- sweep(sweep(data, 1, means, "-"), 1, sds, "/")

    # Estimate batch effects
    t_design <- t(design)
    effects <- solve(t_design %*% design) %*% t_design %*% t(data)


    # Empirical Bayes
    effect_mean <- rowMeans(effects, na.rm = TRUE)
    effect_variance <- apply(effects, 1, var, na.rm = TRUE)
    effect_variance[is.na(effect_variance)] <- 0

    # Calculate weights
    weights <- (effect_variance * effects + effect_mean) / as.vector(effect_variance + 1)


    # Correct data
    data <- data - t(design %*% weights)

    # Reverse the standardization
    data <- sweep(sweep(data, 1, sds, "*"), 1, means, "+")

    # Return the batch-corrected data
    return(data)
}







#
# library(untargetedms)
#
# project <- setProjectFolder(folder = r"(F:\testData\VOILA LUMC)", projectName = "POS")
#
#
#
# exp <- readExperiment(project)
#
# m <- assay(exp, "rt")
#
# peaks <- readPeakData(project, hiveQuery = list(mass = 130)) %>%
#   arrange(aliquot)
#
#
#
#
# idx <- which(peaks$aliquot == unique(peaks$aliquot)[30])[1]:nrow(peaks)
# peaks$batch[idx] <- "batch_04"
#
# idx <- which(peaks$aliquot == unique(peaks$aliquot)[60])[1]:nrow(peaks)
# peaks$batch[idx] <- "batch_05"
# idx <- which(peaks$aliquot == unique(peaks$aliquot)[90])[1]:nrow(peaks)
# peaks$batch[idx] <- "batch_06"
# peaks$rtOrig <- peaks$rt
# peaks$rt[idx] <- peaks$rt[idx] + 2.3 * runif(length(idx), 0.9, 1.1)
#
# peaks <- peaks %>%
#   arrange(targetId)
#
#
# length(unique(peaks$targetId))
#
# m <- matrix(peaks$rt, nrow = length(unique(peaks$targetId)))
#
# batches <- peaks$batch[1:ncol(m)]
#
# m2 <- combat(m, batches)
#
#
#
# batches
# peaks$rtCorrected <- as.vector(m2)
# peaks$rtDiff <- peaks$rt / peaks$rtCorrected
# peaks[idx, ]  %>%
#   select(rtOrig, rt, rtCorrected, rtDiff)
#
# a
# plot(as.vector(m), as.vector(m2) - as.vector(m))
#
# idx <- c(4:6)
#
# m <- matrix(c(t(m[idx, ])), nrow = 1)
# m
#
#
# batches <- rep(1:length(idx), ncol(m) / length(idx))
#
# m
# batches
#
#
# m[is.na(m)] <- rowMeans(m, na.rm = T)
# m2 <- combat(m, batches)
#
# m2
#
#
# matplot(m)
#
# matplot(m2)
# matplot(m - m2)
#
#
# data <- data.frame(
#   value = (c(m / m2) - 1) * 100,
#   compound = rownames(exp),
#   aliquot = rep(colnames(exp), each = nrow(exp)),
#   batch = rep(batches, each = nrow(exp))
# )
#
#
# library(ggplot2)
#
# ggplot(data, aes(x = compound, y = value, group = batch, color = batch)) + geom_point(show.legend = F)
#
#
#
#
# matplot(data - m2)
#
# matplot(m2)
