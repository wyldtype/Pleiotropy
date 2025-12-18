sapply(c("tidyr", "dplyr", "ggplot2", "ggpubr", "purrr", "readxl", "readr", 
         "stringr", "WGCNA", "cluster", "DESeq2", "circlize", "stats"), FUN = require,
       character.only = TRUE)

# loading counts, infodf, and getter functions
load("data/ImmuneCounts.RData")

# Removing every timepoint after 24 hours (>12) for clustering
# Rationale: Schlamp did the same. The primary expression responses happen
# in the first 24 hours and we want clustering to be based on these
counts24 <- counts[,setdiff(colnames(counts), c("13A", "14A", "15A", "16A", "17A", "18A", "19A", "20A", "21A",
                                                "13B", "14B", "15B", "16B", "17B", "18B", "19B", "20B", "21B"))]

#### QC: Dorsal/Cactus/CrebA ####
# Their expression is very different between rep A and B
# But it seems to be different in the same way
# maybe it's a clue for how to group genes
cactus_idx <- "FBgn0000250" # Cactus
dorsal_idx <- "FBgn0260632" # Dorsal
creba_idx <- "FBgn0004396" # CrebA
plotdf <- tibble(sample = colnames(counts24)[-1],
                 cactus = as.numeric(counts24[counts24$gene_name == cactus_idx, -1]),
                 dorsal = as.numeric(counts24[counts24$gene_name == dorsal_idx, -1]),
                 creba = as.numeric(counts24[counts24$gene_name == creba_idx, -1])) |> 
  pivot_longer(cols = c("cactus", "dorsal", "creba"), names_to = "gene", values_to = "expr")
plotdf$timepoint <- parse_number(plotdf$sample)
plotdf$hours <- c(0, 1, 2, 4, 6, 8, 10, 12, 14, 16, 20, 24, 30, 36, 42, 48, 72, 96, 120)[plotdf$timepoint]
plotdf$replicate <- if_else(grepl("A", plotdf$sample), true = "A", false = "B")
ggplot(plotdf, aes(x = factor(hours), y = expr)) +
  geom_point(aes(color = factor(gene)),
             size = 0.25) +
  geom_line(aes(color = factor(gene), 
                group = interaction(gene, replicate),
                linetype = replicate)) +
  #scale_color_brewer(palette = "Set1") +
  scale_x_discrete(breaks = c(0, 2, 4, 6, 12, 24, 48, 120)) +
  labs(color = "Cluster") +
  theme_classic() +
  ylab("Expression (tpm)") +
  xlab("Hours after injection") +
  facet_wrap(~factor(gene))
# Rep B is a little higher, especially at later timepoints

#### Exploration: smoothing by moving average ####
getMovingAverage <- function(.cts) {
  cts_movavg <- map(.x = colnames(.cts), \(cond) {
    idx <- which(colnames(.cts) == cond)
    idxs <- c(idx - 2,
              idx - 1,
              idx,
              idx + 1,
              idx + 2)[c(idx - 2,
                         idx - 1,
                         idx,
                         idx + 1,
                         idx + 2) > 0 &
                         c(idx - 2,
                           idx - 1,
                           idx,
                           idx + 1,
                           idx + 2) < ncol(.cts)]
    if (ncol(.cts[, idxs, drop = FALSE]) < 2) {
      warning("only one timepoint for", cond, idx, "\n")
      return(.cts[, idxs])
    }
    return(rowMeans(.cts[, idxs, drop = FALSE]))
  }) |> purrr::reduce(.f = cbind)
  colnames(cts_movavg) <- colnames(.cts)
  rownames(cts_movavg) <- rownames(.cts)
  return(cts_movavg)
}
# tests for getMovingAverage
# gene_idx <- sample(counts24$gene_name, 1)
gene_idx <- creba_idx
test_counts24A <- counts24[counts24$gene_name == gene_idx, 
                       grepl("A", colnames(counts24)), 
                       drop = FALSE]
test_counts24B <- counts24[counts24$gene_name == gene_idx, 
                       grepl("B", colnames(counts24)), 
                       drop = FALSE]
test_movavgA <- getMovingAverage(test_counts24A)
test_movavgB <- getMovingAverage(test_counts24B)
plotdf <- bind_rows(tibble(before = as.numeric(test_counts24A),
                           after = as.numeric(test_movavgA),
                           timepoint = parse_number(colnames(test_counts24A)),
                           replicate = "A"),
                    tibble(before = as.numeric(test_counts24B),
                           after = as.numeric(test_movavgB),
                           timepoint = parse_number(colnames(test_counts24B)),
                           replicate = "B")) |> 
  pivot_longer(cols = c("before", "after"),
               names_to = "status", values_to = "expr")
ggplot(plotdf, aes(x = timepoint, y = expr)) +
  geom_line(aes(group = interaction(status, replicate), 
                color = status, linetype = replicate))
# issue: moving averages summarise too much, 
# always finding midpoint obscures real, responsive expression

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
gene_idx <- creba_idx
test_countsA <- counts24[counts24$gene_name == gene_idx, 
                       grepl("A", colnames(counts24)), 
                       drop = FALSE]
test_countsB <- counts24[counts24$gene_name == gene_idx, 
                       grepl("B", colnames(counts24)), 
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
colnamesA <- grep("A", colnames(counts24), value = TRUE)
colnamesB <- grep("B", colnames(counts24), value = TRUE)
setequal(parse_number(colnamesA), parse_number(colnamesB)) # should be equal

# A-B Partition
countsA <- counts24[, c("gene_name", colnamesA)]
countsB <- counts24[, c("gene_name", colnamesB)]

# smoothing
smoothMatA <- smoothCountMatrix(.cts = countsA[,-1])
smoothMatB <- smoothCountMatrix(.cts = countsB[,-1])
smoothA <- cbind(countsA[,1, drop = FALSE], smoothMatA)
colnames(smoothA) <- c("gene_name", colnamesA)
smoothB <- cbind(countsB[,1, drop = FALSE], smoothMatB)
colnames(smoothB) <- c("gene_name", colnamesB)
smooth_full <- cbind(smoothA, smoothB[,-1])

#### Clustering ####
# @input: counts matrix
# and whether to return the heirarchical tree as well as cluster assignments
# @output: dataframe of same nRow as nRow counts matrix assigning each
# gene to a cluster 1, 2, ...nClust (nClust currently determined by cutreeHybrid)
# if .tree_too = TRUE, output is a list with elements tree and df
corCluster <- function(.cts, .tree_too = FALSE) {
  # checking if any genes have entirely NAs for counts (it's ok if they have some)
  if (sum(apply(.cts, 1, \(x) {all(is.na(x))})) != 0) { # if any rows (genes) have all NA values, they will cause cor to fail below
    cat("NA genes in counts matrix, returning counts matrix only\n")
    return(.cts)
  }
  # checking for 0-count genes
  cts_mat <- .cts[,-1] |> as.matrix()
  zero_idxs <- which(rowSums(cts_mat) == 0)
  if (length(zero_idxs) > 0) {
    cat(length(zero_idxs), "genes with all 0 counts, removing them!\n")
    zero_genes <- .cts[zero_idxs,1]
    cts_mat <- cts_mat[-zero_idxs,]
    .cts <- .cts[-zero_idxs,]
  }
  # clustering
  rownames(cts_mat) <- unlist(.cts[,1]) |> as.character()
  cor_mat <- cts_mat |> t() |> cor(use = "pairwise.complete.obs")
  dist_mat <- as.dist(-cor_mat)
  htree <- hclust(dist_mat, method = "average")
  # deciding where to cut tree
  tree_labels <- cutreeHybrid(dendro = htree, distM = as.matrix(dist_mat), minClusterSize = 100, deepSplit = FALSE)
  cat("cutting tree into", length(unique(tree_labels$labels)), "clusters\n")
  labeldf <- tibble(gene_name = htree$labels,
                    label = as.numeric(tree_labels$labels))
  # rearranging genes into their original order
  outdf <- left_join(tibble(gene_name = rownames(cts_mat)),
                     labeldf, by = "gene_name")
  if (!.tree_too) {
    return(outdf)
  }
  if (.tree_too) {
    return(list(tree = htree, df = outdf))
  }
}

# tests for corCluster:
test_clusters <- corCluster(.cts = smooth_full[sample(c(1:nrow(counts)), size = 1000),], .tree_too = TRUE)
plotDendroAndColors(dendro = test_clusters$tree, colors = test_clusters$df$label, dendroLabels = FALSE)

#### QC: A-B split ####
# clustering A and B replicates separately to decide which genes have
# robust clustering (which are presumably the ones that have dynamic and
# reproducible expression responses during the Imd challenge)
clustA <- corCluster(.cts = smoothA, .tree_too = TRUE)
plotDendroAndColors(dendro = clustA$tree, colors = clustA$df$label[clustA$df$label != 0], dendroLabels = FALSE)
clustB <- corCluster(.cts = smoothB, .tree_too = TRUE)
plotDendroAndColors(dendro = clustB$tree, colors = clustB$df$label, dendroLabels = FALSE)
# Use matchLabels to make them as similar as possible
# before matchLabels
table(clustA$df$label, clustB$df$label)

# matching
matchedB <- matchLabels(source = clustB$df$label,
                        reference = clustA$df$label)
# after matching
table(clustA$df$label, matchedB)

### PCA to merge similar modules
MEs_A <- moduleEigengenes(expr = t(smoothA[,-1]),
                          colors = clustA$df$label, nPC = 1) # excludeGrey = TRUE removes unclustered (0) MEs, but we need them for our method
MEs_B <- moduleEigengenes(expr = t(smoothB[,-1]),
                          colors = clustB$df$label, nPC = 1)
# comparing MEs and avgExpr
plotModuleExpressionShapes <- function(.ME_object,
                                       .MEs_or_AvgExpr = "AvgExpr",
                                       .separate_replicates = TRUE) {
  if (.MEs_or_AvgExpr == "MEs") {
    expr_mat <- .ME_object$eigengenes
    ylabel <- "Module eigengene expression"
  }
  else {
    expr_mat <- .ME_object$averageExpr
    ylabel <- "Average cluster expression (log2)"
  }
  plotdf <- expr_mat |>
    mutate(timepoint = rownames(expr_mat)) |> 
    pivot_longer(cols = colnames(expr_mat),
                 names_to = "label",
                 values_to = "expr")
  if (.separate_replicates) {
    plotdf$replicate <- if_else(grepl(x = plotdf$timepoint, pattern = "A"),
                                true = "A", false = "B")
  }
  plotdf$timepoint <- parse_number(plotdf$timepoint)
  plotdf$hours <- c(0, 1, 2, 4, 6, 8, 10, 12, 14, 16, 20, 24, 30, 36, 42, 48, 72, 96, 120)[plotdf$timepoint]
  p <- ggplot(plotdf, aes(x = factor(hours), y = expr)) +
    geom_point(aes(color = factor(label)),
               size = 0.25) +
    #scale_color_brewer(palette = "Set1") +
    scale_x_discrete(breaks = c(0, 2, 4, 6, 12, 24, 48, 120)) +
    labs(color = "Cluster") +
    theme_classic() +
    ylab(ylabel) +
    xlab("Hours after injection") +
    facet_wrap(~factor(label))
  if (.separate_replicates) {
    p <- p + geom_line(aes(color = factor(label), 
                           group = interaction(label, replicate),
                           linetype = replicate))
  }
  if (!.separate_replicates) {
    p <- p + geom_line(aes(color = factor(label), 
                           group = label))
  }
  return(p)
}
# A
p_meA <- plotModuleExpressionShapes(MEs_A, .MEs_or_AvgExpr = "MEs")
p_avgA <- plotModuleExpressionShapes(MEs_A, .MEs_or_AvgExpr = "AvgExpr")
ggarrange(p_meA, p_avgA, nrow = 1, ncol = 2, legend = "none")
# B
p_meB <- plotModuleExpressionShapes(MEs_B, .MEs_or_AvgExpr = "MEs")
p_avgB <- plotModuleExpressionShapes(MEs_B, .MEs_or_AvgExpr = "AvgExpr")
ggarrange(p_meB, p_avgB, nrow = 1, ncol = 2, legend = "none")
# conclusion: ME and AvgExpr are very simlar. Interesting that the first principal component can be very related to the column average
# A vs B clusters by avgExpr:
ggarrange(p_avgA, p_avgB, nrow = 1, ncol = 2, legend = "none")
clustA$df |> filter(gene_name %in% c(dorsal_idx, cactus_idx, creba_idx))
clustB$df |> filter(gene_name %in% c(dorsal_idx, cactus_idx, creba_idx))
# conclusion: fewer clusters in B, but they are similar shapes and very smooth compared to non-smooth

### PCA-based merging
# combining AB clustering results prior to PCA
AvgExpr_A <- MEs_A$averageExpr
rownames(AvgExpr_A) <- parse_number(rownames(AvgExpr_A))
colnames(AvgExpr_A) <- paste0(colnames(AvgExpr_A), "_A")
AvgExpr_B <- MEs_B$averageExpr
rownames(AvgExpr_B) <- parse_number(rownames(AvgExpr_B))
colnames(AvgExpr_B) <- paste0(colnames(AvgExpr_B), "_B")
AvgExpr_AB <- cbind(AvgExpr_A, AvgExpr_B)
### pca
getPCA <- function(.ME_object, .MEs_or_AvgExpr = "AvgExpr") {
  if (.MEs_or_AvgExpr == "AvgExpr") {
    covmat <- .ME_object$averageExpr |> cor()
    pca <- prcomp(covmat)
  }
  if (.MEs_or_AvgExpr == "MEs") {
    covmat <- .ME_object$eigengenes |> cor()
    pca <- prcomp(covmat)
  }
  return(pca)
}
pca_AB <- getPCA(.ME_object = list(averageExpr = AvgExpr_AB),
                 .MEs_or_AvgExpr = "AvgExpr")

# how many PCs should we keep?
plotPCAScree <- function(.PCA) {
  plot(x = c(1:length(.PCA$sdev)),
       y = .PCA$sdev/sum(.PCA$sdev),
       xlab = "PC", ylab = "% var explained")
}
plotPCAScree(pca_AB) # first 4 PCs
nPCs <- 4

# what do they look like?
plotPCs <- function(.PCA, .xPC = 1, .yPC = 2) {
  pcadf <- tibble(pcx = .PCA$x[,paste0("PC", .xPC)], 
                  pcy = .PCA$x[,paste0("PC", .yPC)],
                  label = rownames(.PCA$x))
  p <- ggplot(pcadf, aes(x = pcx, y = pcy)) +
    geom_text(aes(label = label, color = label)) +
    xlab(paste0("PC", .xPC, " ", round(.PCA$sdev[.xPC]/sum(.PCA$sdev), digits = 2)*100,
                "% of variance")) +
    ylab(paste0("PC", .yPC, " ", round(.PCA$sdev[.yPC]/sum(.PCA$sdev), digits = 2)*100,
                "% of variance")) +
    theme_classic() +
    theme(legend.position = "none")
  return(p)
}
plotPCs(pca_AB, .xPC = 1, .yPC = 2)
plotPCs(pca_AB, .xPC = 2, .yPC = 3) 
plotPCs(pca_AB, .xPC = 3, .yPC = 4) 

# Seem to be some decent grouping. With different clusters grouping in different dimensions

### Heatmap of Euclidean distance in the first X PCs, weighted by variance explained
# Helper function to get distance between two clusters
getDist <- function(.clust1, .clust2, .pca, .nPCs = nPCs,
                    .startPC = 1, .weighted = TRUE) {
  dist_vec <- map(c(.startPC:(.startPC+.nPCs)), \(comp) {
    prop_var <- 1
    if (.weighted == TRUE) {
      prop_var <- .pca$sdev[comp]/sum(.pca$sdev)
    }
    return(prop_var*(.pca$x[.clust1, paste0("PC", comp)] - 
                       .pca$x[.clust2, paste0("PC", comp)])^2)
  }) |> unlist()
  return(sqrt(sum(dist_vec)))
}
# tests for getDist
nPCs <- 2
getDist(.clust1 = "AE1_A", .clust2 = "AE1_A", .pca = pca_AB, .nPCs = nPCs) # should be 0
getDist(.clust1 = "AE6_A", .clust2 = "AE10_B", .pca = pca_AB, .nPCs = nPCs) # should be small
getDist(.clust1 = "AE6_A", .clust2 = "AE7_B", .pca = pca_AB, .nPCs = nPCs) # should be small
getDist(.clust1 = "AE6_A", .clust2 = "AE7_B", .pca = pca_AB, .nPCs = 2, .startPC = 3) # should be small
getDist(.clust1 = "AE18_A", .clust2 = "AE5_B", .pca = pca_AB, .nPCs = nPCs) # should be large

# gets all pairwise distances in a PCA
getDistMat <- function(.PCA, .nPCs) {
  plotdf <- expand_grid(x = rownames(.PCA$x),
                        y = rownames(.PCA$x)) # not a typo to repeat x---the pca x is refering to the matrix not an xy plane
  plotdf$dist <- map2(plotdf$x, plotdf$y, 
                      .f = getDist, 
                      .pca = .PCA, .nPCs = .nPCs, 
                      .weighted = TRUE) |> 
    unlist()
  plotmat <- matrix(nrow = sqrt(nrow(plotdf)), 
                    ncol = sqrt(nrow(plotdf)),
                    plotdf$dist, byrow = TRUE) 
  return(plotmat)
}

# Heatmap
nPCs
plotmat <- getDistMat(pca_AB, .nPCs = nPCs)
rownames(plotmat) <- rownames(pca_AB$x)
colnames(plotmat) <- rownames(pca_AB$x)
heatmap(plotmat, symm = TRUE, col = topo.colors(256),
        cexRow = 0.4, cexCol = 0.4,
        hclustfun = \(x) hclust(x , method = "complete"))
# do these cluster groupings make sense?
ggarrange(p_avgA, p_avgB, nrow = 1, ncol = 2, legend = "none")
# defining groups
par(mfrow = c(1,1))
tree <- hclust(dist(plotmat), method = "complete") # same tree in heatmap
plot(tree)
cutHeight <- 6
abline(h = cutHeight, col = "red")
new_cluster_lookup <- cutree(tree, h = cutHeight)
new_cluster_lookup
# grouping
# A
table(clustA$df$label) # before merging
new_labels_A_PCA <- map(clustA$df$label,
                        \(x) {
                          new_cluster_lookup[paste0("AE", x, "_A")] |> as.numeric() 
                        }) |> unlist()
table(new_labels_A_PCA) # after merging
# B
table(clustB$df$label) # before merging
new_labels_B_PCA <- map(clustB$df$label,
                        \(x) {
                          new_cluster_lookup[paste0("AE", x, "_B")] |> as.numeric() 
                        }) |> unlist()
table(new_labels_B_PCA) # after merging

new_labels_A_PCA[counts$gene_name %in% c(cactus_idx, creba_idx, dorsal_idx)]
new_labels_B_PCA[counts$gene_name %in% c(cactus_idx, creba_idx, dorsal_idx)]

## Merging by WGCNA::mergeCloseModules
new_labels_list_A <- mergeCloseModules(exprData = t(as.matrix(smoothA[,-1])),
                                       colors = clustA$df$label,
                                       cutHeight = 0.5)
new_labels_A_WGCNA <- new_labels_list_A$colors
new_labels_list_B <- mergeCloseModules(exprData = t(as.matrix(smoothB[,-1])),
                                       colors = clustB$df$label,
                                       cutHeight = 0.5)
new_labels_B_WGCNA <- new_labels_list_B$colors
# matching labels, b/c the same number is not necessarily the same shape
new_labels_B_WGCNA <- matchLabels(source = new_labels_B_WGCNA, 
                                  reference = new_labels_A_WGCNA)
new_labels_A_WGCNA[counts$gene_name %in% c(cactus_idx, creba_idx, dorsal_idx)]
new_labels_B_WGCNA[counts$gene_name %in% c(cactus_idx, creba_idx, dorsal_idx)]

## Comparing expression shapes
newMEs_A_PCA <- moduleEigengenes(t(smoothA[,-1]),
                                 colors = new_labels_A_PCA)
newMEs_B_PCA <- moduleEigengenes(t(smoothB[,-1]),
                                 colors = new_labels_B_PCA)
newMEs_A_WGCNA <- moduleEigengenes(t(smoothA[,-1]),
                                   colors = new_labels_A_WGCNA)
newMEs_B_WGCNA <- moduleEigengenes(t(smoothB[,-1]),
                                   colors = new_labels_B_WGCNA)
p_A_PCA <- plotModuleExpressionShapes(.ME_object = newMEs_A_PCA, .MEs_or_AvgExpr = "AvgExpr")
p_B_PCA <- plotModuleExpressionShapes(.ME_object = newMEs_B_PCA, .MEs_or_AvgExpr = "AvgExpr")
p_A_WGCNA <- plotModuleExpressionShapes(.ME_object = newMEs_A_WGCNA, .MEs_or_AvgExpr = "AvgExpr")
p_B_WGCNA <- plotModuleExpressionShapes(.ME_object = newMEs_B_WGCNA, .MEs_or_AvgExpr = "AvgExpr")
# A PCA vs A WGCNA
ggarrange(p_A_PCA, p_A_WGCNA, p_B_PCA, p_B_WGCNA,
          nrow = 2, ncol = 2, legend = "none")
# Are the same genes robust in both definitions?
sum((new_labels_A_PCA == new_labels_B_PCA) &
      (new_labels_A_WGCNA == new_labels_B_WGCNA)) # both robust
sum((new_labels_A_PCA == new_labels_B_PCA) |
      (new_labels_A_WGCNA == new_labels_B_WGCNA)) # total robust
# about half and half overlap

## Assigning new labels
clustdfA <- clustA$df
clustdfB <- clustB$df
# we assume all genes are in the same order in both splits:
sum(clustdfA$gene_name == clustdfB$gene_name)/nrow(clustdfA) # should be 1
# updating labels
clustdfA$label <- new_labels_A_WGCNA
clustdfB$label <- new_labels_B_WGCNA
sum(clustdfA$label == clustdfB$label)/nrow(clustdfA) # our number of robust clusterers
# combining
clustdfAB <- select(clustdfA, gene_name)
clustdfAB$label <- map2(clustdfA$label, clustdfB$label, 
                      .f = \(a, b) {
                        if (a == b) {
                          return(a)
                        }
                        else {
                          return(paste(a, b, sep = "_"))
                        }
                      }) |> 
  unlist()
clustdfAB$robust <- map(clustdfAB$label, 
                      .f = \(x) !grepl("_", x)) |> 
  unlist()
sum(clustdfAB$robust)/nrow(clustdfAB) # % robust
clustdfAB |> filter(gene_name %in% c(dorsal_idx, cactus_idx, creba_idx))
robust_cluster_genes <- clustdfAB |> filter(robust) |> select(gene_name) |> pull()

#### Clustering full data ####
clusters <- corCluster(.cts = smooth_full, .tree_too = TRUE)
plotDendroAndColors(dendro = clusters$tree, colors = clusters$df$label, dendroLabels = FALSE)
# initial shapes:
MEs <- moduleEigengenes(expr = t(smooth_full[,-1]), colors = clusters$df$label)
plotModuleExpressionShapes(MEs, .MEs_or_AvgExpr = "AvgExpr") # premerge

### merging clusters using WGCNA::mergeCloseModules 
new_label_list <- mergeCloseModules(exprData = as.matrix(t(smooth_full[,-1])),
                                    colors = clusters$df$label,
                                    cutHeight = 0.5)
length(table(new_label_list$colors))
new_labels_WGCNA <- new_label_list$colors
newMEs_WGCNA <- moduleEigengenes(expr = t(counts24[,-1]), 
                                 colors = new_labels_WGCNA)

plotModuleExpressionShapes(newMEs_WGCNA, .MEs_or_AvgExpr = "AvgExpr")

### merging clusters by PCA
pca <- getPCA(MEs, .MEs_or_AvgExpr = "AvgExpr")
plotPCAScree(pca) # 3 PCs
nPCs <- 2
plotPCs(pca, .xPC = 1, .yPC = 2)
plotPCs(pca, .xPC = 2, .yPC = 3)
# pick a pair that looks close:
getDist("AE4", "AE13", .pca = pca, .nPCs = nPCs, .weighted = TRUE)
# pick a pair that looks far:
getDist("AE4", "AE5", .pca = pca, .nPCs = nPCs, .weighted = TRUE)

plotmat <- getDistMat(pca, .nPCs = nPCs)
rownames(plotmat) <- rownames(pca$x)
colnames(plotmat) <- rownames(pca$x)
heatmap(plotmat, symm = TRUE, col = topo.colors(256),
        cexRow = 0.4, cexCol = 0.4,
        hclustfun = \(x) hclust(x , method = "complete"))
par(mfrow = c(1,1))
tree <- hclust(dist(plotmat), method = "complete") # same tree in heatmap
plot(tree)
cutHeight <- 3
abline(h = cutHeight, col = "red")
new_cluster_lookup <- cutree(tree, h = cutHeight)
new_cluster_lookup
# grouping
table(clusters$df$label) # before merging
new_labels_PCA <- map(clusters$df$label,
                      \(x) {
                        new_cluster_lookup[paste0("AE", x)] |> as.numeric() 
                      }) |> unlist()
table(new_labels_PCA) # after merging

# recalculate MEs
newMEs_PCA <- moduleEigengenes(t(counts24[,-1]), colors = new_labels_PCA)

## comparing merge results
newMEs <- moduleEigengenes(t(counts24[,-1]), colors = clusters$df$label)
# pre-merge:
plotModuleExpressionShapes(newMEs, .MEs_or_AvgExpr = "AvgExpr") # postmerge
# WGCNA::mergeCloseModules:
plotModuleExpressionShapes(newMEs_WGCNA, .MEs_or_AvgExpr = "AvgExpr") # postmerge
# PCA-distance merge:
plotModuleExpressionShapes(newMEs_PCA, .MEs_or_AvgExpr = "AvgExpr") # postmerge
# conclusion: both seem to work well and come up with similar shapes

## setting new labels
clustdf <- clusters$df
table(dense_rank(new_labels_WGCNA)) # WGCNA groups modules without necessarily keeping the values 1,2,3,etc.
clustdf$label <- new_labels_WGCNA |> dense_rank()

# adding robust column and checking if any clusters are over-represented in robust fraction
clustdf$robust <- clustdf$gene_name %in% robust_cluster_genes
table(clustdf$robust, clustdf$label) # 1 is heavily under-represented, 7 also pretty under-represented

# final MEs on all genes
MEs <- moduleEigengenes(t(counts24[,-1]), 
                               colors = clustdf$label)
plotModuleExpressionShapes(MEs, .MEs_or_AvgExpr = "AvgExpr")
# final MEs on just robust genes
MEs_robust <- moduleEigengenes(t(counts24[counts24$gene_name %in% robust_cluster_genes,-1]), 
                        colors = clustdf$label[clustdf$robust])
plotModuleExpressionShapes(MEs_robust, .MEs_or_AvgExpr = "AvgExpr") # not very different
# checking our favorite genes
clustdf |> filter(gene_name %in% c(dorsal_idx, cactus_idx, creba_idx))

#### Saving ####
save(robust_cluster_genes, clustdf, counts, file = "data/ClusteringImmune.RData")

########################### Archive ################################
# #### Exploration: visualizing cluster similarity in a PCA ####
# 
# # In this section, we explore visualizing module expression and
# # merging modules based on Principal Component Analysis of their
# # expression shape, as measured by average expression of all module genes,
# # or the WGCNA term "module eigengene". 
# # We define functions that are useful for the actual clustering
# 
# # Covariance matrix of module average expression
# MEs <- moduleEigengenes(expr = t(counts[,-1]),
#                         colors = clusters$df$label, nPC = 1)
# 
# ### Are Module Eigengenes similar to average expression?
# p_avg <- plotModuleExpressionShapes(.ME_object = MEs, .MEs_or_AvgExpr = "AvgExpr")
# p_me <- plotModuleExpressionShapes(.ME_object = MEs, .MEs_or_AvgExpr = "MEs")
# ggarrange(p_avg, p_me, nrow = 1, ncol = 2, legend = "none")
# 
# ### PCA based on AvgExpr
# pca <- getPCA(MEs, .MEs_or_AvgExpr = "AvgExpr")
# plotPCAScree(pca) # 4 PCs
# # PC1 and PC2
# plotPCs(pca)
# # PC2 and PC3
# plotPCs(pca, .xPC = 2, .yPC = 3)
# 
# plotmat <- getDistMat(.PCA = pca, .nPCs = 3)
# heatmap(plotmat, symm = TRUE, col = topo.colors(256))
# # not a perfect legend:
# # legend(x="bottomright", legend = round(quantile(plotmat, probs = c(0, 0.5, 1)), digits = 2), fill=topo.colors(3))
# 
# # This gives us a sense of the most similar module groups
# # Big groups (pre limiting to 24h):
# # 13, 4, 19, 8 (long early spike by eye)
# # 17, 9 (flat by eye)
# # 1, 7, 10 (down by eye)
# # 2, 3, 15 (up by eye)
# # 16, 11 (weird by eye)
# # 12, 18
# # 5, 14
# # 6
# # Small, more similar groups:
# # 19, 8
# # 9, 17
# # 1, 10
# # 11, 16
# 
# #### Clustering A and B replicates separately ####
# 
# # 1) Splitting rep A and rep B prior to clustering
# # 2) Principal component analysis of average expression of each A and B cluster
# # 3) Manually merging clusters in A and B based on PCA
# # 4) Identifying the robust-clustering genes based on these manually-merged clusters
# 
# colnamesA <- grep("A", colnames(counts), value = TRUE)
# colnamesB <- grep("B", colnames(counts), value = TRUE)
# setequal(parse_number(colnamesA), parse_number(colnamesB)) # should be equal
# 
# # A-B Partition
# countsA <- counts[, c("gene_name", colnamesA)]
# countsB <- counts[, c("gene_name", colnamesB)]
# 
# # clustering
# clustA <- corCluster(.cts = countsA, .tree_too = TRUE)
# clustB <- corCluster(.cts = countsB, .tree_too = TRUE)
# 
# # Use matchLabels to make them as similar as possible
# # before matchLabels
# table(clustA$df$label, clustB$df$label)
# 
# # matching
# matchedB <- matchLabels(source = clustB$df$label,
#                         reference = clustA$df$label)
# # after matching
# table(clustA$df$label, matchedB)
# 
# ### PCA to merge similar modules
# MEs_A <- moduleEigengenes(expr = t(countsA[,-1]),
#                           colors = clustA$df$label, nPC = 1) # excludeGrey = TRUE removes unclustered (0) MEs, but we need them for our method
# MEs_B <- moduleEigengenes(expr = t(countsB[,-1]),
#                           colors = clustB$df$label, nPC = 1)
# # comparing MEs and avgExpr
# # A
# p_meA <- plotModuleExpressionShapes(MEs_A, .MEs_or_AvgExpr = "MEs")
# p_avgA <- plotModuleExpressionShapes(MEs_A, .MEs_or_AvgExpr = "AvgExpr")
# ggarrange(p_meA, p_avgA, nrow = 1, ncol = 2, legend = "none")
# # B
# p_meB <- plotModuleExpressionShapes(MEs_B, .MEs_or_AvgExpr = "MEs")
# p_avgB <- plotModuleExpressionShapes(MEs_B, .MEs_or_AvgExpr = "AvgExpr")
# ggarrange(p_meB, p_avgB, nrow = 1, ncol = 2, legend = "none")
# # Conclusion: MEs and AvgExpr still very similar.
# # A and B avgExpr:
# ggarrange(p_avgA, p_avgB, nrow = 1, ncol = 2, legend = "none")
# # Conclusion: most clusters that have the same number in A and B
# # are not the same expression shape (except 1, which makes sense)
# #
# # slow decreasers: 1A-1B
# # early downspike: 9A-14B, 7A-9B, 6B
# # early upspike: 5A-11B
# # kind of early upspike: 4A-3B-5B-2A
# # late upspike: 10A-11A-3A-4B-8B
# 
# # combining AB clustering results prior to PCA
# AvgExpr_A <- MEs_A$averageExpr
# rownames(AvgExpr_A) <- parse_number(rownames(AvgExpr_A))
# colnames(AvgExpr_A) <- paste0(colnames(AvgExpr_A), "_A")
# AvgExpr_B <- MEs_B$averageExpr
# rownames(AvgExpr_B) <- parse_number(rownames(AvgExpr_B))
# colnames(AvgExpr_B) <- paste0(colnames(AvgExpr_B), "_B")
# AvgExpr_AB <- cbind(AvgExpr_A, AvgExpr_B)
# pca_AB <- getPCA(.ME_object = list(averageExpr = AvgExpr_AB),
#                  .MEs_or_AvgExpr = "AvgExpr")
# plotPCAScree(pca_AB) # first 4 PCs
# plotPCs(pca_AB, .xPC = 1, .yPC = 2)
# plotPCs(pca_AB, .xPC = 2, .yPC = 3)
# plotPCs(pca_AB, .xPC = 3, .yPC = 4)
# 
# # Heatmap
# plotmat <- getDistMat(pca_AB, .nPCs = 4)
# rownames(plotmat) <- rownames(pca_AB$x)
# colnames(plotmat) <- rownames(pca_AB$x)
# heatmap(plotmat, symm = TRUE, col = topo.colors(256),
#         cexRow = 0.2, cexCol = 0.2)
# 
# # By eye, groups to merge (pre 24h filter):
# # new cluster 1: 1A, 1B, 2B, 8A
# # new cluster 2: 2A, 5B
# # new cluster 3: 6A, 13A, 7B, 
# # new cluster 4: 14A, 10B, 11B
# # new cluster 5: 15A, 18B
# # new cluster 6: 17A, 13B, 0A, 5A, 16B
# # new cluster 7: 17B, 11A
# # new cluster 8: 9B, 15B, 4A, 3B, 7A
# # new cluster 9: 16A, 14B,
# # new cluster 10: 12A, 12B
# # new cluster 11: 10A, 4B, 8B, 9A, 6B, 3A
# new_cluster_lookup <- list("1" = c("1A", "1B", "2B", "8A"),
#                            "2" = c("2A", "5B"),
#                            "3" = c("6A", "13A", "7B"),
#                            "4" = c("14A", "10B", "11B"),
#                            "5" = c("15A", "18B"),
#                            "6" = c("17A", "13B", "0A", "5A", "16B"),
#                            "7" = c("17B", "11A"),
#                            "8" = c("9B", "15B", "4A", "3B", "7A"),
#                            "9" = c("16A", "14B", "12A", "12B"),
#                            "10" = c("10A", "4B", "8B", "9A", "6B", "3A"))
# # make sure every old cluster is represented exactly once (There is no 0B or 18A):
# table(unlist(new_cluster_lookup)) |> sort()
# getNewLabel <- function(.old_label, .lookup, .A_or_B) {
#   .old_label <- paste0(.old_label, .A_or_B)
#   new_label <- lapply(new_cluster_lookup, \(x) .old_label %in% x) |> 
#     unlist() |> which(useNames = TRUE) |> names()
#   return(new_label)
# }
# # tests for getNewLabel
# getNewLabel(.old_label = 6, .lookup = new_cluster_lookup,
#             .A_or_B = "A") # should be 3
# getNewLabel(.old_label = "6", .lookup = new_cluster_lookup,
#             .A_or_B = "A")
# 
# # merging AB clustering
# # A
# clustdfA <- clustA$df
# table(clustdfA$label) # before merging
# clustdfA$label <- map(clustdfA$label,
#                       \(x) getNewLabel(.old_label = as.numeric(x), 
#                                        .lookup = new_cluster_lookup,
#                                        .A_or_B = "A")) |> 
#   unlist()
# table(clustdfA$label) # after merging
# # B
# clustdfB <- clustB$df
# table(clustdfB$label) # before merging
# clustdfB$label <- map(clustdfB$label,
#                       \(x) getNewLabel(.old_label = as.numeric(x), 
#                              .lookup = new_cluster_lookup,
#                              .A_or_B = "B")) |> 
#   unlist()
# table(clustdfB$label) # after merging
# # we assume all genes are in the same order in both splits:
# sum(clustdfA$gene_name == clustdfB$gene_name)/nrow(clustdfA) # should be 1
# sum(clustdfA$label == clustdfB$label)/nrow(clustdfA) # our number of robust clusterers
# # combining
# clustdf <- select(clustdfA, gene_name)
# clustdf$label <- map2(clustdfA$label, clustdfB$label, 
#                       .f = \(a, b) {
#                         if (a == b) {
#                           return(a)
#                         }
#                         else {
#                           return(paste(a, b, sep = "_"))
#                         }
#                       }) |> 
#   unlist()
# clustdf$robust <- map(clustdf$label, 
#                       .f = \(x) !grepl("_", x)) |> 
#   unlist()
# sum(clustdf$robust)/nrow(clustdf) # 54%, should be same percent as above
# clustdf_AB <- clustdf
# 
# #### QC: Repeating with random split to compare robustly-clustering fraction ####
# 
# # 1) Splitting dataset in half randomly prior to clustering (random whether each timepoint comes from rep A or rep B)
# # 2) Principal component analysis of average expression of each A and B cluster
# # 3) Manually merging clusters in A and B based on PCA
# # 4) Identifying the robust-clustering genes based on these manually-merged clusters
# 
# # randomly partitioning replicates at each timepoint:
# sample_idxs <- sample(c(1:length(colnamesA)), 
#                       size = length(colnamesA)/2,
#                       replace = FALSE) |> sort(decreasing = FALSE)
# other_idxs <- setdiff(c(1:length(colnamesA)), sample_idxs)
# colnamesA_random <- c(colnamesA[sample_idxs], colnamesB[other_idxs])
# colnamesB_random <- c(colnamesB[sample_idxs], colnamesA[other_idxs])
# 
# # hardcoding the first random split we got:
# #colnamesA_random <- c("1A", "2A", "3A", "5A", "7B", "8A", "9A", "10B", "11B", "12B", "13A",
#                       # "14B", "15A", "16B", "17A", "18B", "19A", "20B", "21B")
# #colnamesB_random <- setdiff(c(colnamesA, colnamesB), colnamesA_random)
# 
# countsA <- counts[, c("gene_name", colnamesA_random)]
# countsB <- counts[, c("gene_name", colnamesB_random)]
# 
# setequal(x = parse_number(colnames(countsA)[-1]), 
#          y = parse_number(colnames(countsB)[-1])) # should be equal
# 
# # clustering
# clustA <- corCluster(.cts = countsA, .tree_too = TRUE)
# clustB <- corCluster(.cts = countsB, .tree_too = TRUE)
# 
# # Use matchLabels to make them as similar as possible
# # before matchLabels
# table(clustA$df$label, clustB$df$label)
# 
# # matching
# matchedB <- matchLabels(source = clustB$df$label,
#                         reference = clustA$df$label)
# # after matching
# table(clustA$df$label, matchedB)
# 
# ### PCA to merge similar modules
# MEs_A <- moduleEigengenes(expr = t(countsA[,-1]),
#                           colors = clustA$df$label, nPC = 1) # excludeGrey = TRUE removes unclustered (0) MEs, but we need them for our method
# MEs_B <- moduleEigengenes(expr = t(countsB[,-1]),
#                           colors = clustB$df$label, nPC = 1)
# 
# # comparing MEs and avgExpr
# # A
# p_meA <- plotModuleExpressionShapes(MEs_A, .MEs_or_AvgExpr = "MEs", .separate_replicates = FALSE)
# p_avgA <- plotModuleExpressionShapes(MEs_A, .MEs_or_AvgExpr = "AvgExpr", .separate_replicates = FALSE)
# ggarrange(p_meA, p_avgA, nrow = 1, ncol = 2, legend = "none")
# # B
# p_meB <- plotModuleExpressionShapes(MEs_B, .MEs_or_AvgExpr = "MEs", .separate_replicates = FALSE)
# p_avgB <- plotModuleExpressionShapes(MEs_B, .MEs_or_AvgExpr = "AvgExpr", .separate_replicates = FALSE)
# ggarrange(p_meB, p_avgB, nrow = 1, ncol = 2, legend = "none")
# # Conclusion: MEs and AvgExpr still very similar.
# # A and B avgExpr:
# ggarrange(p_avgA, p_avgB, nrow = 1, ncol = 2, legend = "none")
# # By eye, these look to be more similar expression shapes than the non-random A-B split
# 
# # combining AB clustering results prior to PCA
# AvgExpr_A <- MEs_A$averageExpr
# rownames(AvgExpr_A) <- parse_number(rownames(AvgExpr_A))
# colnames(AvgExpr_A) <- paste0(colnames(AvgExpr_A), "_A")
# AvgExpr_B <- MEs_B$averageExpr
# rownames(AvgExpr_B) <- parse_number(rownames(AvgExpr_B))
# colnames(AvgExpr_B) <- paste0(colnames(AvgExpr_B), "_B")
# AvgExpr_AB <- cbind(AvgExpr_A, AvgExpr_B)
# pca_AB <- getPCA(.ME_object = list(averageExpr = AvgExpr_AB),
#                  .MEs_or_AvgExpr = "AvgExpr")
# plotPCAScree(pca_AB) # first 4 PCs
# plotPCs(pca_AB, .xPC = 1, .yPC = 2)
# plotPCs(pca_AB, .xPC = 2, .yPC = 3)
# plotPCs(pca_AB, .xPC = 3, .yPC = 4)
# 
# # Heatmap
# plotmat <- getDistMat(pca_AB, .nPCs = 4)
# rownames(plotmat) <- rownames(pca_AB$x)
# colnames(plotmat) <- rownames(pca_AB$x)
# heatmap(plotmat, distfun = dist, symm = TRUE, col = topo.colors(256),
#         cexRow = 0.4, cexCol = 0.4, hclustfun = \(x) hclust(x , method = "average"))
# tree <- hclust(dist(plotmat), method = "average") # same tree in heatmap
# plot(tree)
# cutHeight <- 5
# abline(h = cutHeight, col = "red")
# new_cluster_lookup <- cutree(tree, h = cutHeight)
# table(new_cluster_lookup)
# # merging AB clustering
# # A
# clustdfA <- clustA$df
# table(clustdfA$label) # before merging
# clustdfA$label <- map(clustdfA$label,
#                       \(x) {
#                         new_cluster_lookup[paste0("AE", x, "_A")] |> as.numeric() 
#                       }) |> unlist()
# table(clustdfA$label) # after merging
# # B
# clustdfB <- clustB$df
# table(clustdfB$label) # before merging
# clustdfB$label <- map(clustdfB$label,
#                       \(x) {
#                         new_cluster_lookup[paste0("AE", x, "_B")] |> as.numeric() 
#                       }) |> unlist()
# table(clustdfB$label) # after merging
# # we assume all genes are in the same order in both splits:
# sum(clustdfA$gene_name == clustdfB$gene_name)/nrow(clustdfA) # should be 1
# sum(clustdfA$label == clustdfB$label)/nrow(clustdfA) # our number of robust clusterers
# # combining
# clustdf <- select(clustdfA, gene_name)
# clustdf$label <- map2(clustdfA$label, clustdfB$label, 
#                       .f = \(a, b) {
#                         if (a == b) {
#                           return(a)
#                         }
#                         else {
#                           return(paste(a, b, sep = "_"))
#                         }
#                       }) |> 
#   unlist()
# clustdf$robust <- map(clustdf$label, 
#                       .f = \(x) !grepl("_", x)) |> 
#   unlist()
# sum(clustdf$robust)/nrow(clustdf) # 35%, fewer than with non-random A-B split
# clustdf_random <- clustdf
# rm(clustdf)
# 
# #### Are the same genes robust in both splits? ####
# table(clustdf_AB$robust, clustdf_random$robust) 
# sum(clustdf_AB$robust == clustdf_random$robust)/nrow(clustdf_AB) # 70% still have their same label
# clustdf_AB |> filter(gene_name == "FBgn0004396") # CrebA
# clustdf_random |> filter(gene_name == "FBgn0004396") # CrebA
# clustdf_AB |> filter(gene_name == "FBgn0000250") # Cactus
# clustdf_random |> filter(gene_name == "FBgn0000250") # Cactus
# clustdf_AB |> filter(gene_name == "FBgn0260632") # Dorsal
# clustdf_random |> filter(gene_name == "FBgn0260632") # Dorsal
# 
# #### Merging similar modules ####
# # Issue: WGCNA doesn't have a way to post-hoc merge small modules into their closest bigger module
# is_robust <- counts$gene_name %in% robust_cluster_genes
# pca <- moduleEigengenes(expr = t(counts[is_robust,-1]),
#                         colors = clusters$df$label[is_robust], nPC = 2)
# # Further issue: WGCNA doesn't output more than PC1, so I'll need to make my own
# # function to obtain PC2. But also it might not be what I'm looking for---I think I'm 
# # more looking for covariance between average module expression
# 
# # TODO: genes to monitor during merging close modules part
# "FBgn0004396" # CrebA (non-robust 16)
# "FBgn0000250" # Cactus (robust 16)
# "FBgn0260632" # Dorsal (robust 13, which looks an awful lot like cluster 16)
# 
# # Merging close modules
# newlabels <- mergeCloseModules(exprData = t(counts[,-1]), 
#                                colors = clusters$df$label, 
#                                cutHeight = 0.4) # cutHeight is 1-cor that is our threshold for merging, so cutHeight of 0.2 causes modules to be merged if their eigengenes have >= 0.8 cor
# plotDendroAndColors(dendro = clusters$tree, 
#                     colors = cbind(clusters$df$label, newlabels$colors),
#                     groupLabels = c("unmerged", "merged"),
#                     dendroLabels = FALSE)
# # checking if Dorsal and Cactus are put together
# cactus_idx <- which(counts$gene_name == "FBgn0000250") # Cactus
# dorsal_idx <- which(counts$gene_name == "FBgn0260632") # Dorsal
# # not together in pre-merge:
# clusters$df$label[cactus_idx] |> as.numeric()
# clusters$df$label[dorsal_idx] |> as.numeric()
# # are they together in post-merge?
# newlabels$colors[cactus_idx] |> as.numeric()
# newlabels$colors[dorsal_idx] |> as.numeric() # Even when a lot of merging is happening, they don't get put together. Perhaps they are not as similar expression shapes as they look
# 
# #### QC: how often do genes from the two replicates cluster together? ####
# # A-B Partition
# countsA <- counts[, c("gene_name", colnamesA)]
# countsB <- counts[, c("gene_name", colnamesB)]
# 
# # clustering
# clustA <- corCluster(.cts = countsA, .tree_too = TRUE)
# clustB <- corCluster(.cts = countsB, .tree_too = TRUE)
# 
# # Use matchLabels to make them as similar as possible
# # before matchLabels
# table(clustA$df$label, clustB$df$label)
# 
# # matching
# matchedB <- matchLabels(source = clustB$df$label,
#                         reference = clustA$df$label)
# # after matching
# table(clustA$df$label, matchedB)
# 
# # taking the genes that cluster the same in both replicates
# #robust_cluster_idxs_random <- which(clustA$df$label == matchedB)
# robust_cluster_idxs_AB <- which(clustA$df$label == matchedB)
# 
# # Random Replicate Partition
# # splitting into two halves, randomly partitioning counts from A or B at each timepoint
# sample_idxs <- sample(c(1:length(colnamesA)), 
#                       size = length(colnamesA)/2,
#                       replace = FALSE) |> sort(decreasing = FALSE)
# other_idxs <- setdiff(c(1:length(colnamesA)), sample_idxs)
# countsA <- counts[, c("gene_name", colnamesA[sample_idxs], colnamesB[other_idxs])]
# countsB <- counts[, c("gene_name", colnamesB[sample_idxs], colnamesA[other_idxs])]
# 
# # clustering
# clustA <- corCluster(.cts = countsA, .tree_too = TRUE)
# clustB <- corCluster(.cts = countsB, .tree_too = TRUE)
# 
# # Use matchLabels to make them as similar as possible
# # before matchLabels
# table(clustA$df$label, clustB$df$label)
# 
# # matching
# matchedB <- matchLabels(source = clustB$df$label,
#                         reference = clustA$df$label)
# # after matching
# table(clustA$df$label, matchedB)
# 
# robust_cluster_idxs_random <- which(clustA$df$label == matchedB)
# 
# length(intersect(robust_cluster_idxs_random, robust_cluster_idxs_AB))
# length(c(robust_cluster_idxs_random, robust_cluster_idxs_AB))
# 
# # our core set of genes that cluster together in both the A-B split and in the random replicate partition
# robust_cluster_genes <- counts$gene_name[intersect(robust_cluster_idxs_random, robust_cluster_idxs_AB)]
# # slight differences in exact number of genes, but typically about 2000 genes
# 
# #### Saving ####
# save(robust_cluster_genes, clusters, file = "data/Clustering.RData")
# 
# ###################### Archive ########################
# #### manually setting nClust (k) ####
# # # gap statistic to determine appropriate number of clusters
# # calculateGapStat <- function(.cts) {
# #   # checking if any genes have entirely NAs for counts (it's ok if they have some)
# #   if (sum(apply(.cts, 1, \(x) {all(is.na(x))})) != 0) { # if any rows (genes) have all NA values, they will cause cor to fail below
# #     cat("NA genes in counts matrix, returning counts matrix only\n")
# #     return(.cts)
# #   }
# #   # converting counts df to a matrix
# #   cts_mat <- .cts[,-1] |> as.matrix()
# #   rownames(cts_mat) <- unlist(.cts[,1]) |> as.character()
# #   # calculating gap statistic using clusGap
# #   cluster_fun <- function(x, k) {
# #     list(cluster = cutree(hclust(d = as.dist(-cor(t(x), use = "pairwise.complete.obs")),
# #                                  method = "average"),
# #                           k = k))}
# #   gap_stat <- clusGap(cts_mat, FUNcluster = cluster_fun, K.max = 30, B = 5)
# #   return(gap_stat)
# # }
# # gap_test <- calculateGapStat(counts)
# # optimal_k <- maxSE(gap_test$Tab[, "gap"], 
# #                    gap_test$Tab[, "SE.sim"], 
# #                    method = "globalSEmax")
# # # Depending on what method you use, optimal k ranges from 1 to 2 to 8 to 
# # # 23, so definitely will need to research more what is going on here
# # optimal_k
