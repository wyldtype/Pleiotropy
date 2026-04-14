sapply(c("tidyr", "dplyr", "ggplot2", "ggpubr", "ggalluvial", "purrr", "readxl", "readr", "stringr", "plotly", "lmtest", "pheatmap"), FUN = require,
       character.only = TRUE)

# loading counts and functions
source("utils.r")

# loading HOG data
dmel <- read_tsv("../phylogeny/215Fly_Orthogroups_tables/Dmel_HOG_association.tsv")
HOG_to_OG <- read_tsv("../phylogeny/215Fly_Orthogroups_tables/HOG_to_OG_table.tsv")
ortho <- read_tsv("../phylogeny/215Fly_Orthogroups_tables/Orthogroups.tsv")
ortho_GeneCounts <- read_tsv("../phylogeny/215Fly_Orthogroups_tables/Orthogroups.GeneCount.tsv")

# TODO: Use the above data to identify duplications that aren't 
# found across all species (I think there's are rows in ortho
# where there are NA values for certain species)

# Why do I think that? 93% of NAs in Dmel or Dsim are NAs in both:
sum(is.na(ortho$DROSOPHILA_MELANOGASTER) & 
      is.na(ortho$DROSOPHILA_SIMULANS))/sum(is.na(ortho$DROSOPHILA_MELANOGASTER) |
                                              is.na(ortho$DROSOPHILA_SIMULANS))

# 
