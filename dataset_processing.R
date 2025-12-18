sapply(c("tidyr", "dplyr", "ggplot2", "purrr", "readxl", "readr", 
         "stringr"), FUN = require, character.only = TRUE)

#### Reading in datasets ####
# Imd challenge, Schlamp et al.
immune_counts <- read_csv("data/normalized_counts_log_filt.csv")
colnames(immune_counts) <- c("gene_name", colnames(immune_counts)[-1])

# Embryonic dev, Becker et al. 
file_path <- "data/GEO_Becker/GSE121160_RAW/"
dev_counts <- list.files(path = file_path, 
                         pattern = "*.tsv", full.names = FALSE) |> 
  map(.f = \(x) {
    full_x <- paste0("data/GEO_Becker/GSE121160_RAW/", x)
    outdf <- read_delim(full_x, delim = "\t", 
                        col_names = c("gene_name", gsub(".tsv", "", x)))
    return(outdf)
  }) |> purrr::reduce(.f = left_join, by = "gene_name")

# Oogenesis, Tarikere et al.
colnames_oog_counts_GermSoma <- read_tsv("data/Geo_Tarikere/GSE172015_Germ_Soma_table_counts_wide_ctrl.tsv",
                                         n_max = 0, skip = 0) |> colnames()
oog_counts_GermSoma <- read_tsv("data/Geo_Tarikere/GSE172015_Germ_Soma_table_counts_wide_ctrl.tsv",
                                col_names = FALSE, skip = 1)
colnames(oog_counts_GermSoma) <- c("gene_name", colnames_oog_counts_GermSoma)
colnames_oog_counts_WholeOvary <- read_tsv("data/Geo_Tarikere/GSE172015_WholeOvary_table_counts_wide.tsv",
                                         n_max = 0, skip = 0) |> colnames()
oog_counts_WholeOvary <- read_tsv("data/Geo_Tarikere/GSE172015_WholeOvary_table_counts_wide.tsv",
                                col_names = FALSE, skip = 1)
colnames(oog_counts_WholeOvary) <- c("gene_name", colnames_oog_counts_WholeOvary)
# combining GermSoma and WholeOvary
dim(oog_counts_GermSoma)
dim(oog_counts_WholeOvary)
sum(oog_counts_GermSoma$gene_name == oog_counts_WholeOvary$gene_name)
oog_counts <- left_join(oog_counts_GermSoma, oog_counts_WholeOvary,
                        by = "gene_name")

#### Normalizing to tpm ####
# libsize
plotdf <- bind_rows(tibble(libsize = colSums(immune_counts[,-1]),
                           dataset = "schlamp"),
                    tibble(libsize = colSums(dev_counts[,-1]),
                           dataset = "becker"),
                    tibble(libsize = colSums(oog_counts[,-1]),
                           dataset = "tarikere"))
ggplot(filter(plotdf, dataset == "schlamp"), aes(x = libsize)) + 
  geom_density(aes(fill = dataset), alpha = 0.5)
ggplot(filter(plotdf, dataset == "becker"), aes(x = libsize)) + 
  geom_density(aes(fill = dataset), alpha = 0.5)
ggplot(filter(plotdf, dataset == "tarikere"), aes(x = libsize)) + 
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
oog_counts_tpm <- apply(oog_counts[,-1], 2, \(x) {
  libsize <- sum(x)
  out_vec <- round((x/libsize)*1e6)
  return(out_vec)
})
dev_counts_tpm <- cbind(dev_counts[,1], dev_counts_tpm)
oog_counts_tpm <- cbind(oog_counts[,1], oog_counts_tpm)
immune_counts_tpm <- cbind(immune_counts[,1], immune_counts_tpm)

#### Filtering for lowly expressed genes ####
# before filtering lowly expressed in dev
sum(rowMeans(dev_counts_tpm[,-1]) < 1)
sum(rowMeans(oog_counts_tpm[,-1]) < 1)
sum(rowMeans(immune_counts_tpm[,-1]) < 1)
# filtering lowly expressed
# dev
keep_genes <- rowMeans(dev_counts_tpm[,-1]) >= 1
sum(keep_genes)
dev_counts_tpm <- dev_counts_tpm[keep_genes,]
# oog
keep_genes <- rowMeans(oog_counts_tpm[,-1]) >= 1
sum(keep_genes)
oog_counts_tpm <- oog_counts_tpm[keep_genes,]
# immune
keep_genes <- rowMeans(immune_counts_tpm[,-1]) >= 1
sum(keep_genes)
immune_counts_tpm <- immune_counts_tpm[keep_genes,]
# after filtering lowly expressed in dev
sum(rowMeans(dev_counts_tpm[,-1]) < 1)
sum(rowMeans(oog_counts_tpm[,-1]) < 1)
sum(rowMeans(immune_counts_tpm[,-1]) < 1)

plotdf <- bind_rows(tibble(libsize = colSums(immune_counts_tpm[,-1]),
                           dataset = "schlamp"),
                    tibble(libsize = colSums(dev_counts_tpm[,-1]),
                           dataset = "becker"),
                    tibble(libsize = colSums(oog_counts_tpm[,-1]),
                           dataset = "tarikere"))
ggplot(plotdf, aes(x = libsize)) + 
  geom_density(aes(fill = dataset), alpha = 0.5) 
# slight differences in libsize b/c immune already had lowly expressed filtered out,
# but much more comparable

#### Limiting to common set of genes ####
common_genes <- intersect(immune_counts_tpm[,1], dev_counts_tpm[,1])
common_genes <- intersect(common_genes, oog_counts_tpm[,1])
length(common_genes)
# also removing one annoying gene that becomes 0 in the immune counts rep A
"FBgn0028519" %in% common_genes # if it's in the set, remove it with the line below:
# common_genes <- setdiff(common_genes, "FBgn0028519")

# restricting to common genes
immune_counts_tpm <- immune_counts_tpm[immune_counts_tpm[,1] %in% common_genes,]
dim(immune_counts_tpm)
dev_counts_tpm <- dev_counts_tpm[dev_counts_tpm[,1] %in% common_genes,]
dim(dev_counts_tpm)
oog_counts_tpm <- oog_counts_tpm[oog_counts_tpm[,1] %in% common_genes,]
dim(oog_counts_tpm)

#### Dataset-specific sample processing & saving ####

#### Imd challenge ####
# no getter functions, no infodf for immune (b/c colnames are so interpretable)
counts <- immune_counts_tpm

#### Sample PCA to look for outliers ####
pca <- prcomp(cor(counts[,-1])) # fyi: eigen() and prcomp() produce different eigenvalues and vectors, and I don't know why (scaling and rotation I suspect)
plot(x = c(1:length(pca$sdev)),
     y = pca$sdev/sum(pca$sdev),
     xlab = "PC", ylab = "% var explained") # 4 PCs
plotdf <- tibble(pc1 = pca$x[,"PC1"], pc2 = pca$x[,"PC2"],
                 pc3 = pca$x[,"PC3"], pc4 = pca$x[,"PC4"],
                 label = colnames(counts)[-1])
plotdf$timepoint <- parse_number(plotdf$label)
plotdf$hour <- c(0, 1, 2, 4, 6, 8, 10, 12, 14, 16, 20, 24, 30, 36, 42, 48, 72, 96, 120)[plotdf$timepoint]
# PC1 and PC2
ggplot(plotdf, aes(x = pc1, y = pc2)) +
  geom_text(aes(label = label, color = timepoint)) +
  geom_line(aes(group = timepoint)) +
  xlab(paste0("PC1 ", round(pca$sdev[1]/sum(pca$sdev), digits = 2)*100,
              "% of variance")) +
  ylab(paste0("PC2 ", round(pca$sdev[2]/sum(pca$sdev), digits = 2)*100,
              "% of variance")) +
  theme_classic() +
  theme(legend.position = "none") # 6A is an outlier sample in PC2

# Checking for missing timepoints
colnamesA <- grep("A", colnames(counts), value = TRUE)
colnamesB <- grep("B", colnames(counts), value = TRUE) 
colnamesA
colnamesB
setdiff(parse_number(colnamesA), parse_number(colnamesB)) 
setdiff(parse_number(colnamesB), parse_number(colnamesA))
# conclusion: B is missing timepoint 4

# removing timepoint 6 (outlier), 4 (missing)
counts <- counts[,setdiff(colnames(counts), c("6A", "6B", "4A"))]

save(counts, file = "data/ImmuneCounts.RData")

#### Embryogenesis ####
counts <- dev_counts_tpm

info_matrix <- read_delim(file = "data/GEO_Becker/GSE121160_series_matrix.txt",
                          delim = "\t", col_names = c("name", "content")) # this is not ideal, b/c some rows have tabs delim each samples values and others don't and they're not consecutive
info_list <- as.list(info_matrix$content)
names(info_list) <- info_matrix$name 
infodf <- tibble(condition = info_list$`!Sample_title` |> trimws() |> 
                   str_split(pattern = "\\s") |> 
                   unlist(),
                 sample_name = info_list$`!Series_sample_id` |> trimws() |> 
                   str_split(pattern = "\\s") |> 
                   unlist())
infodf$hour <- str_split_i(infodf$condition, pattern = "_", i = 1) |> 
  parse_number()
infodf$replicate <- str_split_i(infodf$condition, pattern = "_", i = 2)
infodf$colname <- paste(infodf$sample_name, infodf$condition, sep = "_")

# getter functions
getHour <- function(.name, .info = infodf) {
  if (all(grepl("^[0-9+]*$", .name)))
    return(as.numeric(.name))
  out_vec <- purrr::map(.name, \(nm) {
    .info |> filter(colname == nm) |> 
      select(hour) |> pull()
  }) |> unlist()
  return(out_vec)
}
# test for getHour
colnames(counts)[1:10]
getHour(colnames(counts)[1])
getHour(colnames(counts)[10])
getHour(colnames(counts)[2:9])
getHour(c("GSM3427149_18h_2", "GSM3427150_18h_3", "GSM3427153_20h_2"))
getHour(c("0", "3", "4"))

getReplicate <- function(.name, .info = infodf) {
  out_vec <- purrr::map(.name, \(nm) {
    .info |> filter(colname == nm) |> 
      select(replicate) |> pull()
  }) |> unlist()
  return(out_vec)
}
# test for getReplicate
colnames(counts)[1:10]
getReplicate(colnames(counts)[1])
getReplicate(colnames(counts)[9])
getReplicate(colnames(counts)[2:9])
getReplicate(c("GSM3427149_18h_2", "GSM3427150_18h_3", "GSM3427153_20h_2"))

### Checking for missing samples
table(getReplicate(colnames(counts)[-1]),
      getHour(colnames(counts)[-1])) # nothing missing

#### Dev: Sample PCA to look for outliers ####
pca <- prcomp(cor(counts[,-1]))
plot(x = c(1:length(pca$sdev)),
     y = pca$sdev/sum(pca$sdev),
     xlab = "PC", ylab = "% var explained") # 4 PCs
plotdf <- tibble(pc1 = pca$x[,"PC1"], pc2 = pca$x[,"PC2"],
                 pc3 = pca$x[,"PC3"], pc4 = pca$x[,"PC4"],
                 label = colnames(counts)[-1])
plotdf$hour <- getHour(colnames(counts[-1]))
plotdf$replicate <- getReplicate(colnames(counts[-1]))
plotdf$condition <- paste(plotdf$hour, plotdf$replicate, sep = "_")
# PC1 and PC2
ggplot(plotdf, aes(x = pc1, y = pc2)) +
  geom_text(aes(label = condition, color = hour)) +
  geom_line(aes(group = hour)) +
  xlab(paste0("PC1 ", round(pca$sdev[1]/sum(pca$sdev), digits = 2)*100,
              "% of variance")) +
  ylab(paste0("PC2 ", round(pca$sdev[2]/sum(pca$sdev), digits = 2)*100,
              "% of variance")) +
  theme_classic() +
  theme(legend.position = "none") 
# PC2 captures dev time
# 16hr rep1 has a lot of the PC1 variation
# other outliers: 3hr rep2, 1hr rep4
# PC2 and PC3
ggplot(plotdf, aes(x = pc3, y = pc2)) +
  geom_text(aes(label = condition, color = hour)) +
  geom_line(aes(group = hour)) +
  xlab(paste0("PC3 ", round(pca$sdev[3]/sum(pca$sdev), digits = 2)*100,
              "% of variance")) +
  ylab(paste0("PC2 ", round(pca$sdev[2]/sum(pca$sdev), digits = 2)*100,
              "% of variance")) +
  theme_classic() +
  theme(legend.position = "none")

# Based on PCA, reps1 and 2 are no closer related than 3 and 4,
# so should be fine to split it 1-2 and 3-4

#### Dev: Removing outlier samples ####
to_remove <- c("GSM3427103_00h_4", "GSM3427113_03h_2", "GSM3427144_16h_1")
counts <- counts[,setdiff(colnames(counts), to_remove)]

#### Dev: Saving ####
save(counts, infodf, getHour, getReplicate, file = "data/DevCounts.RData")

#### Oogenesis ####
counts <- oog_counts_tpm

# info
info_matrix_GermSoma <- read_delim(file = "data/Geo_Tarikere/GSE172015_Germ_Soma_table_counts_metadata_long_ctrl.tsv",
                                   delim = "\t", col_select = 1:5) |> 
  unique()
info_matrix_WholeOvary <- read_delim(file = "data/Geo_Tarikere/GSE172015_WholeOvary_metadata_long.tsv",
                                     delim = "\t", col_select = 2:6) |> 
  unique()
colnames(info_matrix_GermSoma) <- c("sample_name", "cell_type", "treatment", "stage", "replicate")
colnames(info_matrix_WholeOvary) <- c("sample_name", "cell_type", "treatment", "stage", "replicate")
infodf <- bind_rows(info_matrix_GermSoma, info_matrix_WholeOvary)
infodf$hour <- if_else(infodf$stage == "Early",
                       true = 72, false = if_else(infodf$stage == "Mid",
                                                  true = 96, false = 120)) # Hours After Egg Laying (the egg that is developing into the larva in which oogenesis is measured)

# getter functions
getHour <- function(.name, .info = infodf) {
  parsed_name_list <- purrr::map(.name, \(nm) str_split(nm, pattern = "_"))
  out_vec <- purrr::map(parsed_name_list, \(st) {
    st <- unlist(st)
    if(("Early" %in% st) | (72 %in% st)) {
      return(72)
    }
    if(("Mid" %in% st) | (96 %in% st)) {
      return(96)
    }
    if(("Late" %in% st) | (120 %in% st)) {
      return(120)
    }
    else {
      return(NA)
    }
  }) |> unlist()
  return(out_vec)
}
# test for getHour
colnames(counts)[1:10]
getHour(colnames(counts)[1]) # should be NA
getHour("early") # should be NA
getHour(colnames(counts)[2])
getHour(colnames(counts)[10])
getHour(colnames(counts)[2:9])
getHour(c("G_Ctrl_Early_3", "G_Ctrl_Late_2", "G_Ctrl_Mid_3"))
getHour(c("Early_3", "Late_2", "Mid_3"))

getReplicate <- function(.name, .info = infodf) {
  purrr::map(.name, \(nm) {
    n_segments <- str_split(nm, pattern = "_") |> unlist() |> 
      length()
    rep_name <- str_split_i(nm, pattern = "_", i = n_segments) |> unlist()
    if (rep_name == "name") {
      return(NA)
    }
    return(rep_name)
  }) |> unlist()
}
# test for getReplicate
colnames(counts)[1:10]
getReplicate(colnames(counts)[1]) # should be NA
getReplicate("G_Early_3") # should be 3 b/c replicate is the last item in the colname (unlike getHour and getCellType, which look for specific items
getReplicate(colnames(counts)[2])
getReplicate(colnames(counts)[9])
getReplicate(colnames(counts)[2:9])
getReplicate(c("G_Ctrl_Early_3", "G_Ctrl_Late_2", "G_Ctrl_Mid_1"))

getCellType <- function(.name, .info = infodf) {
  parsed_name_list <- purrr::map(.name, \(nm) str_split(nm, pattern = "_"))
  out_vec <- purrr::map(parsed_name_list, \(st) {
    st <- unlist(st)
    if(("G" %in% st) | ("Germ" %in% st)) {
      return("Germ")
    }
    if(("S" %in% st) | ("Somatic" %in% st)) {
      return("Somatic")
    }
    if("WholeOvary" %in% st) {
      return("WholeOvary")
    }
    else {
      return(NA)
    }
  }) |> unlist()
  return(out_vec)
}
# test for getCellType
colnames(counts)[1:10]
getCellType(colnames(counts)[1]) # should be NA
getCellType(colnames(counts)[2])
getCellType(colnames(counts)[c(2, 13, 20)])
getCellType(c("G_Ctrl_Early_1", "S_Ctrl_Early_3", "WholeOvary_Ctrl_Early_1"))

### Checking for missing samples
table(getReplicate(colnames(counts)[-1]),
      getHour(colnames(counts)[-1]),
      getCellType(colnames(counts)[-1])) # 4th replicate added only in Germ, to 96 and 120 to bring them up to 3 reps per timepoint

#### Oog: Renaming replicates ####
# Seeing as rep1 from the three timepoints isn't the same fly (there are 20-30 ovaries per sample),
# we will rename rep4 in Germ so that all three cell types have reps 1,2,3

# renaming Germ rep4 Mid to Germ rep2 Mid, and rename Germ rep4 Late to Germ rep1 Late
# renaming counts
colnames(counts) <- gsub(x = colnames(counts),
                         pattern = "G_Ctrl_Mid_4", 
                         "G_Ctrl_Mid_2")
colnames(counts) <- gsub(x = colnames(counts),
                         pattern = "G_Ctrl_Late_4", 
                         "G_Ctrl_Late_1")
# renaming infodf
infodf[infodf$replicate == 4 & infodf$stage == "Mid",]$replicate <- 2
infodf[infodf$replicate == 4 & infodf$stage == "Late",]$replicate <- 1

# re-check table
table(getReplicate(colnames(counts)[-1]),
      getHour(colnames(counts)[-1]),
      getCellType(colnames(counts)[-1])) # 3 reps per timepoint/cell type

### Oog: Sample PCA to look for outliers
pca <- prcomp(cor(counts[,-1]))
plot(x = c(1:length(pca$sdev)),
     y = pca$sdev/sum(pca$sdev),
     xlab = "PC", ylab = "% var explained") # 3 PCs
plotdf <- tibble(pc1 = pca$x[,"PC1"], pc2 = pca$x[,"PC2"],
                 pc3 = pca$x[,"PC3"], pc4 = pca$x[,"PC4"],
                 label = colnames(counts)[-1])
plotdf$hour <- getHour(colnames(counts[-1]))
plotdf$replicate <- getReplicate(colnames(counts[-1]))
plotdf$cell_type <- getCellType(colnames(counts[-1]))
plotdf$condition <- paste(plotdf$cell_type, plotdf$hour, plotdf$replicate, sep = "_")
# PC1 and PC2
ggplot(plotdf, aes(x = pc1, y = pc2)) +
  geom_text(aes(label = condition, color = hour)) +
  geom_line(aes(group = interaction(hour, cell_type))) +
  xlab(paste0("PC1 ", round(pca$sdev[1]/sum(pca$sdev), digits = 2)*100,
              "% of variance")) +
  ylab(paste0("PC2 ", round(pca$sdev[2]/sum(pca$sdev), digits = 2)*100,
              "% of variance")) +
  theme_classic() +
  theme(legend.position = "none") # Germ 96h is all over the place
# facet by cell type
ggplot(plotdf, aes(x = pc1, y = pc2)) +
  geom_text(aes(label = condition, color = hour)) +
  geom_line(aes(group = hour)) +
  xlab(paste0("PC1 ", round(pca$sdev[1]/sum(pca$sdev), digits = 2)*100,
              "% of variance")) +
  ylab(paste0("PC2 ", round(pca$sdev[2]/sum(pca$sdev), digits = 2)*100,
              "% of variance")) +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(~cell_type)
# PC1 is basically just variation between Germ replicates. Seeing as this was also
# The cell type that required extra replicates, it was probably the most challenging
# tissue to collect consistently. Germ_120_1 looks fine (one of the rep 4s we renamed).
# Germ_96_2 is a little more out there, but so is Germ_72_1, so I think this is just
# the reality of collecting the Germ cells. We won't remove any samples, but we'll
# take the extensive variation between replicates in the Germ condition into account

#### Oog: Saving ####
save(counts, infodf, getHour, getReplicate, getCellType, file = "data/OogCounts.RData")

