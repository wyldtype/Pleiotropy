sapply(c("tidyr", "dplyr", "ggplot2", "ggpubr", "purrr", "readxl", "readr", 
         "stringr", "WGCNA", "cluster", "DESeq2", "circlize", "stats"), FUN = require,
       character.only = TRUE)

counts <- read_csv("data/normalized_counts_log_filt.csv")
colnames(counts) <- c("gene_name", colnames(counts)[-1])

#### Sample PCA to look for outliers ####
pca <- prcomp(cor(counts[,-1])) # TODO when you care to: eigen() and plotting vectors didn't work the same and I don't know what the difference between eigen and prcomp is (scaling and rotation I think)
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

#### Excluding outlier/missing samples ####
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

# ..and every timepoint after 24 hours (>12)
# counts <- counts[,setdiff(colnames(counts), c("6A", "6B", "4A",
#                                               "13A", "14A", "15A", "16A", "17A", "18A", "19A", "20A", "21A",
#                                               "13B", "14B", "15B", "16B", "17B", "18B", "19B", "20B", "21B"))]

#### Smoothing by spline ####
deg_free <- 5

# applying to entire count matrix
smoothCountMatrix <- function(.cts, .df = deg_free) {
  ngenes <- nrow(.cts)
  outlist <- map(c(1:ngenes), \(g) {
    cat(g, "/", ngenes, "\n")
    gene_vec <- .cts[g,] |> as.numeric()
    smooth_y <- smooth.spline(x = parse_number(colnames(.cts)),
                              y = gene_vec,
                              df = .df)
    return(smooth_y$y)
  })
  if (length(outlist) == 1) {
    outmat <- purrr::reduce(outlist[[1]], .f = cbind) # I have no idea why reducing doesn't work on single-element lists vectors. Reduce is the same problem
    # giving first timepoint its unsmoothed counts back
    first_idx <- which.min(parse_number(colnames(.cts)))
    outmat[, first_idx] <- as.numeric(.cts[, first_idx])
  }
  if (length(outlist) > 1) {
    outmat <- purrr::reduce(outlist, .f = rbind) 
    # giving first timepoint its unsmoothed counts back
    first_idx <- which.min(parse_number(colnames(.cts)))
    outmat[, first_idx] <- as.numeric(unlist(.cts[, first_idx]))
  }
  return(outmat)
}

# tests for smoothCountMatrix
cactus_idx <- "FBgn0000250" # Cactus
dorsal_idx <- "FBgn0260632" # Dorsal
creba_idx <- "FBgn0004396" # CrebA
gene_idx <- creba_idx
test_countsA <- counts[counts$gene_name == gene_idx, 
                       grepl("A", colnames(counts)), 
                       drop = FALSE]
test_countsB <- counts[counts$gene_name == gene_idx, 
                       grepl("B", colnames(counts)), 
                       drop = FALSE]
test_smoothA <- smoothCountMatrix(.cts = test_countsA)
test_smoothB <- smoothCountMatrix(.cts = test_countsB)
plotdf <- bind_rows(tibble(before = as.numeric(test_countsA),
                           after = as.numeric(test_smoothA),
                           timepoint = parse_number(colnames(test_countsA)),
                           replicate = "A"),
                    tibble(before = as.numeric(test_countsB),
                           after = as.numeric(test_smoothB),
                           timepoint = parse_number(colnames(test_countsB)),
                           replicate = "B")) |> 
  pivot_longer(cols = c("before", "after"),
               names_to = "status", values_to = "expr")
ggplot(plotdf, aes(x = timepoint, y = expr)) +
  geom_line(aes(group = interaction(status, replicate), 
                color = status, linetype = replicate))

# Smoothing
colnamesA <- grep("A", colnames(counts), value = TRUE)
colnamesB <- grep("B", colnames(counts), value = TRUE)
setequal(parse_number(colnamesA), parse_number(colnamesB)) # should be equal

# A-B Partition
countsA <- counts[, c("gene_name", colnamesA)]
countsB <- counts[, c("gene_name", colnamesB)]

# smoothing
smoothMatA <- smoothCountMatrix(.cts = countsA[,-1])
smoothMatB <- smoothCountMatrix(.cts = countsB[,-1])
smoothA <- cbind(countsA[,1, drop = FALSE], smoothMatA)
colnames(smoothA) <- c("gene_name", colnamesA)
smoothB <- cbind(countsB[,1, drop = FALSE], smoothMatB)
colnames(smoothB) <- c("gene_name", colnamesB)
smooth_full <- cbind(smoothA, smoothB[,-1])

#### QC: A-B split ####
# kmeans clustering the first x PCs in repA and repB
# then using the same pipeline as the hierarchical clustering
# to first merge clusters based on a second PCA of
# cluster average expression
# then check for cluster overlap between replicates

?kmeans

