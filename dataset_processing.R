sapply(c("tidyr", "dplyr", "ggplot2", "ggpubr", "purrr", "readxl", "readr", 
         "stringr", "edgeR"), FUN = require, character.only = TRUE)

#### Reading in datasets ####
immune_counts <- read_csv("data/normalized_counts_log_filt.csv")
colnames(immune_counts) <- c("gene_name", colnames(immune_counts)[-1])

file_path <- "data/GEO_Becker/GSE121160_RAW/"
dev_counts <- list.files(path = file_path, 
                         pattern = "*.tsv", full.names = FALSE) |> 
  map(.f = \(x) {
    full_x <- paste0("data/GEO_Becker/GSE121160_RAW/", x)
    outdf <- read_delim(full_x, delim = "\t", 
                        col_names = c("gene_name", gsub(".tsv", "", x)))
    return(outdf)
  }) |> purrr::reduce(.f = left_join, by = "gene_name")

#### Normalizing to tpm ####
# libsize
plotdf <- bind_rows(tibble(libsize = colSums(immune_counts[,-1]),
                           dataset = "schlamp"),
                    tibble(libsize = colSums(dev_counts[,-1]),
                           dataset = "becker"))
ggplot(filter(plotdf, dataset == "schlamp"), aes(x = libsize)) + 
  geom_density(aes(fill = dataset), alpha = 0.5)
ggplot(filter(plotdf, dataset == "becker"), aes(x = libsize)) + 
  geom_density(aes(fill = dataset), alpha = 0.5)
# repeat with tpm
immune_counts_unlogged <- 2^immune_counts[,-1] - 1
immune_counts_tpm <- apply(immune_counts_unlogged, 2, \(x) {
  libsize <- sum(x)
  out_vec <- round((x/libsize)*1e6)
  return(out_vec)
})
dev_counts_tpm <- apply(dev_counts[,-1], 2, \(x) {
  libsize <- sum(x)
  out_vec <- round((x/libsize)*1e6)
  return(out_vec)
})
dev_counts_tpm <- cbind(dev_counts[,1], dev_counts_tpm)
immune_counts_tpm <- cbind(immune_counts[,1], immune_counts_tpm)

#### Filtering for lowly expressed genes ####
# before filtering lowly expressed in dev
sum(rowMeans(dev_counts_tpm[,-1]) < 1)
sum(rowMeans(immune_counts_tpm[,-1]) < 1)
# filtering lowly expressed
# dev
keep_genes <- rowMeans(dev_counts_tpm[,-1]) >= 1
dev_counts_tpm <- dev_counts_tpm[keep_genes,]
# immune
keep_genes <- rowMeans(immune_counts_tpm[,-1]) >= 1
immune_counts_tpm <- immune_counts_tpm[keep_genes,]
# after filtering lowly expressed in dev
sum(rowMeans(dev_counts_tpm[,-1]) < 1)
sum(rowMeans(immune_counts_tpm[,-1]) < 1)

plotdf <- bind_rows(tibble(libsize = colSums(immune_counts_tpm[,-1]),
                           dataset = "schlamp"),
                    tibble(libsize = colSums(dev_counts_tpm[,-1]),
                           dataset = "becker"))
ggplot(plotdf, aes(x = libsize)) + 
  geom_density(aes(fill = dataset), alpha = 0.5) 
# slight differences in libsize b/c immune already had lowly expressed filtered out,
# but much more comparable

#### Limiting to common set of genes ####
common_genes <- intersect(immune_counts_tpm[,1], dev_counts_tpm[,1])
length(common_genes)
# also removing one annoying gene that becomes 0 in the immune counts rep A
common_genes <- setdiff(common_genes, "FBgn0028519")

immune_counts_tpm <- immune_counts_tpm[immune_counts_tpm[,1] %in% common_genes,]
dim(immune_counts_tpm)
dev_counts_tpm <- dev_counts_tpm[dev_counts_tpm[,1] %in% common_genes,]
dim(dev_counts_tpm)

#### Saving ####
counts <- immune_counts_tpm
save(counts, file = "data/ImmuneCounts.RData")
counts <- dev_counts_tpm
save(counts, file = "data/DevCounts.RData")

