sapply(c("tidyr", "dplyr", "ggplot2", "ggpubr", "ggalluvial", "purrr", "readxl", "readr", "stringr"), FUN = require,
       character.only = TRUE)

williams <- read_xlsx(path = "data/17Jan2023_Supplemental_File_1.xlsx",
                      sheet = 2, na = "NA") |> 
  select(FlyBaseID, Plei_allImmune_allDev) 
# using the least-strict definition of pleiotropy,
# because their strict threshold was based on 
# expression, which we can make our own judgments on

load("data/ClusteringImmune.RData")
counts_immune <- counts
clustdf_immune <- clustdf
infodf_immune <- tibble("sample_name" = colnames(counts_immune)[-1],
                        "replicate" = if_else(grepl("A", colnames(counts_immune)[-1]),
                                                    true = "A", false = "B"),
                        "hour" = parse_number(colnames(counts_immune)[-1]))
rm(counts, clustdf, robust_cluster_genes)

load("data/ClusteringDev.RData")
counts_dev <- counts
clustdf_dev <- clustdf
infodf_dev <- infodf |> rename("sample_name"="colname")
rm(counts, clustdf, robust_cluster_genes, infodf)

load("data/ClusteringOog.RData")
counts_oog <- counts
clustdf_oog <- clustdf
infodf_oog <- infodf
rm(counts, clustdf, robust_cluster_genes, infodf)

#### Limiting Oogenesis to WholeOvary ####
samples_wholeOvary <- grep("^WholeOvary", colnames(counts_oog), value = TRUE)
counts_oog <- counts_oog[, c("gene_name", samples_wholeOvary)]
infodf_oog <- filter(infodf_oog, sample_name %in% samples_wholeOvary)

#### Combining datasets ####
# combining clustering results
clustdf_immune$environment <- "immune"
clustdf_dev$environment <- "dev"
clustdf_oog$environment <- "oog"
clustdf <- bind_rows(clustdf_immune, clustdf_dev, clustdf_oog) |> 
  select(gene_name, label, environment) |> 
  left_join(y = select(williams, FlyBaseID, 
                       Plei_allImmune_allDev),
            by = c("gene_name"="FlyBaseID"))
clustdf$williams_category <- if_else(is.na(clustdf$Plei_allImmune_allDev),
                                     true = "Other",
                                     false = clustdf$Plei_allImmune_allDev)

# combining counts and infodfs
counts <- list("immune" = counts_immune,
               "dev" = counts_dev,
               "oog" = counts_oog)
infodf <- list("immune" = infodf_immune,
               "dev" = infodf_dev,
               "oog" = infodf_oog)

#### Parsing sample metadata ####
# getter function for any column in the infodf that is needed
# obtained by sample name
getInfo <- function(.sample_vec, .environment, .info_colname) {
  info_env <- infodf[[.environment]]
  mask <- sapply(.sample_vec, \(s) which(info_env$sample_name == s))
  out_vec <- info_env[mask, .info_colname] |> unlist() |> as.numeric()
  return(out_vec)
}
# test for getInfo
test_col <- colnames(counts_immune)[sample(2:ncol(counts_immune), size = 1)]
test_col
getInfo(test_col, .environment = "immune", .info_colname = "hour")
test_cols <- colnames(counts_oog)[sample(2:ncol(counts_oog), size = 5)]
test_cols
getInfo(test_cols, .environment = "oog", .info_colname = "hour")

#### Functions for Signal to Noise ####
# We'll define "signal" as the proportion of the expression variance (ESS,
# explained sum of squares) captured
# in a linear model with time as the predictor with a third degree polynomial fit
# The third degree polynomial should capture expression tending to increase
# or decrease over time, but not necessarily with one single peak or dip
# like it would with a second degree polynomial
# "noise" is therefore the residual sum of squares (RSS) variance
# Very old school simple linear regression: Total SS = ESS + RSS
# So our signal to noise ratio will be ESS/RSS. If we want it to be from 0
# to 1, we can change to R^2: ESS/TSS. But I think we don't really want this
# because we want the genes with more variability to still have that variability

### Helper for modeling, get the counts and time data for the model from
# the gene name and environment (for looping through clustdf)
getData <- function(.gene, .environment) {
  cts <- counts[[.environment]]
  counts_vec <- cts[cts$gene_name == .gene, -1] |> as.numeric()
  hour_vec <- getInfo(colnames(cts)[-1], 
                      .environment = .environment,
                      .info_colname = "hour")
  return(list(counts = counts_vec, hour = hour_vec))
}

### Fit our polynomial linear model
fitModel <- function(.counts_vec, .hour_vec, .degree = 3) {
  x <- data.frame(count = .counts_vec,
                  hour = .hour_vec)
  mod <- lm(count ~ poly(hour, degree = .degree), data = x)
  return(mod)
}
# tests for fitModel and getData
test_data <- getData(.gene = "FBgn0034329", .environment = "immune") # Bomanin Short 1
test_mod <- fitModel(.counts_vec = test_data$counts, 
                     .hour_vec = test_data$hour)
test_anova <- anova(test_mod)
test_anova[1,"Sum Sq"]/test_anova["Residuals", "Sum Sq"]

### get the anova info
getSNR <- function(.mod) {
  mod_anova <- anova(.mod)
  snr <- mod_anova[1,"Sum Sq"]/mod_anova["Residuals", "Sum Sq"]
  return(snr)
}
# tests for getSNR
test_anova <- anova(test_mod)
test_anova[1,"Sum Sq"]/test_anova["Residuals", "Sum Sq"]
getSNR(test_mod) # should be the same as above

#### Calculating Signal to Noise ####
clustdf$signal_to_noise <- map2(clustdf$gene_name, clustdf$environment,
                                \(g, e) {
                                  cat(which(paste(g, e) == paste(clustdf$gene_name, clustdf$environment)),
                                      "/", nrow(clustdf), "\n")
                                  e_degree <- if_else(e == "oog",
                                                      true = 2, false = 3) # Oog only has 3 time categories, polynomial can only be degree 2
                                  g_data <- getData(.gene = g, .environment = e)
                                  mod <- fitModel(.counts_vec = g_data$counts,
                                                  .hour_vec = g_data$hour,
                                                  .degree = e_degree)
                                  return(getSNR(mod))
                                }) |> unlist()
