sapply(c("tidyr", "dplyr", "ggplot2", "ggpubr", "purrr", "readxl", "readr", "stringr"), FUN = require,
       character.only = TRUE)

williams <- read_xlsx(path = "data/17Jan2023_Supplemental_File_1.xlsx",
                      sheet = 2, na = "NA") 

load("data/ClusteringImmune.RData")
counts_immune <- counts
clustdf_immune <- clustdf
robust_immune_genes <- robust_cluster_genes
rm(counts, clustdf, robust_cluster_genes)

load("data/ClusteringDev.RData")
counts_dev <- counts
clustdf_dev <- clustdf
robust_dev_genes <- robust_cluster_genes
rm(counts, clustdf, robust_cluster_genes)

# combining Immune/dev results
# cluster dataframe
clustdf_immune$immune_or_dev <- "immune"
clustdf_immune$robust <- clustdf_immune$gene_name %in% robust_immune_genes
clustdf_dev$immune_or_dev <- "dev"
clustdf_dev$robust <- clustdf_dev$gene_name %in% robust_dev_genes
clustdf <- bind_rows(clustdf_immune, clustdf_dev) |> 
  left_join(y = select(williams, FlyBaseID, 
                       Plei_immuneResponse_embDev,
                       Plei_allImmune_allDev),
            by = c("gene_name"="FlyBaseID"))
clustdf$wiliams_category <- if_else(is.na(clustdf$Plei_allImmune_allDev),
                                    true = "Other",
                                    false = clustdf$Plei_allImmune_allDev)

# count matrix
colnames(counts_immune) <- c("gene_name", paste0(colnames(counts_immune[-1]), "_schlamp"))
colnames(counts_dev) <- c("gene_name", paste0(colnames(counts_dev[-1]), "_becker"))
sum(counts_immune$gene_name == counts_dev$gene_name)/nrow(counts_immune) # should be 1
counts <- bind_cols(counts_immune, counts_dev[,-1])

#### Parsing sample metadata ####
getHour <- function(.name, .info = infodf) {
  if (all(grepl("^[0-9+]*$", .name)))
    return(as.numeric(.name))
  if (grepl("_becker$", .name)) {
    out_vec <- purrr::map(.name, \(nm) {
      nm <- gsub("_becker$", "", nm)
      .info |> filter(colname == nm) |> 
        select(hour) |> pull()
    }) |> unlist()
  }
  if (grepl("_schlamp$", .name)) {
    out_vec <- purrr::map(.name, \(nm) {
      nm <- gsub("_becker$", "", nm)
      parse_number(nm)
    }) |> unlist()
  }
  return(out_vec)
}
# test for getHour
test_col <- colnames(counts)[sample(2:ncol(counts), size = 1)]
test_col
getHour(test_col)

getReplicate <- function(.name, .info = infodf) {
  if (grepl("_becker$", .name)) {
    out_vec <- purrr::map(.name, \(nm) {
      nm <- gsub("_becker$", "", nm)
      .info |> filter(colname == nm) |> 
        select(replicate) |> pull()
    }) |> unlist()
  }
  if (grepl("_schlamp$", .name)) {
    out_vec <- purrr::map(.name, \(nm) {
      nm <- gsub("_becker$", "", nm)
      if_else(grepl("A", nm), true = "A", false = "B")
    }) |> unlist()
  }
  return(out_vec)
}
# test for getReplicate
test_col <- colnames(counts)[sample(2:ncol(counts), size = 1)]
test_col
getReplicate(test_col)

#### Which genes are robust in both sets? ####


