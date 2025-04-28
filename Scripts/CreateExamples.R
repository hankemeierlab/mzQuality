combined <- buildCombined(system.file("example.tsv", package = "mzQuality2"))

combined$compound <- paste0("Compound", as.integer(as.factor(combined$compound)))
combined$compound_is <- paste0("Standard", as.integer(as.factor(combined$compound_is)))
combined$batch <- combined$batch - 7

exp <- buildExperiment(combined)


# Add concentrations
concentrations <- read.delim(system.file(
    package = "mzQuality",
    "concentrations.txt"
))
concentrations

exp <- addConcentrations(exp, concentrations)

exp
exp <- exp %>%
    doAnalysis()


saveRDS(exp, file.path('inst', "data.RDS"))

library(mzQuality2)




