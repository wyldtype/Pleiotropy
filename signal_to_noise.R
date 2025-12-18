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
robust_immune_genes <- robust_cluster_genes
rm(counts, clustdf, robust_cluster_genes)

load("data/ClusteringDev.RData")
counts_dev <- counts
clustdf_dev <- clustdf
robust_dev_genes <- robust_cluster_genes
rm(counts, clustdf, robust_cluster_genes)

load("data/ClusteringOog.RData")
counts_oog <- counts
clustdf_oog <- clustdf
robust_oog_genes <- robust_cluster_genes
rm(counts, clustdf, robust_cluster_genes)