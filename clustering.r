sapply(c("tidyr", "dplyr", "ggplot2", "ggpubr", "purrr", "readxl", "readr", 
         "stringr", "WGCNA", "cluster", "DESeq2", "circlize"), FUN = require,
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

# removing timepoint 6 (outlier) and 4 (missing)
counts <- counts[,setdiff(colnames(counts), c("6A", "6B", "4A"))]

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
  # clustering
  cts_mat <- .cts[,-1] |> as.matrix()
  rownames(cts_mat) <- unlist(.cts[,1]) |> as.character()
  cor_mat <- cts_mat |> t() |> cor(use = "pairwise.complete.obs")
  dist_mat <- as.dist(-cor_mat)
  htree <- hclust(dist_mat, method = "average")
  # deciding where to cut tree
  tree_labels <- cutreeHybrid(dendro = htree, distM = as.matrix(dist_mat), minClusterSize = 100, deepSplit = FALSE)
  cat("cutting tree into", length(unique(tree_labels$labels)), "clusters\n")
  labeldf <- tibble(gene_name = htree$labels,
                    label = tree_labels$labels)
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

# testing on the whole thing b/c why not
clusters <- corCluster(.cts = counts, .tree_too = TRUE)
plot(clusters$tree, labels = FALSE) # looks like clustering! Definitely more than 5 clusters
plotDendroAndColors(dendro = clusters$tree, colors = clusters$df$label, dendroLabels = FALSE)

#### Exploration: visualizing cluster similarity in a PCA ####

# In this section, we explore visualizing module expression and
# merging modules based on Principal Component Analysis of their
# expression shape, as measured by average expression of all module genes,
# or the WGCNA term "module eigengene". 
# We define functions that are useful for the actual clustering

# Covariance matrix of module average expression
MEs <- moduleEigengenes(expr = t(counts[,-1]),
                        colors = clusters$df$label, nPC = 1)

### Are Module Eigengenes similar to average expression?
plotModuleExpressionShapes <- function(.ME_object,
                                       .MEs_or_AvgExpr = "AvgExpr") {
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
  plotdf$replicate <- if_else(grepl(x = plotdf$timepoint, pattern = "A"),
                            true = "A", false = "B")
  plotdf$timepoint <- parse_number(plotdf$timepoint)
  plotdf$hours <- c(0, 1, 2, 4, 6, 8, 10, 12, 14, 16, 20, 24, 30, 36, 42, 48, 72, 96, 120)[plotdf$timepoint]
  p <- ggplot(plotdf, aes(x = factor(hours), y = expr)) +
    geom_line(aes(color = factor(label), 
                  group = interaction(label, replicate),
                  linetype = replicate)) +
    geom_point(aes(color = factor(label), 
                   group = interaction(label, replicate)),
               size = 0.25) +
    #scale_color_brewer(palette = "Set1") +
    scale_x_discrete(breaks = c(0, 2, 4, 6, 12, 24, 48, 120)) +
    labs(color = "Cluster") +
    theme_classic() +
    ylab(ylabel) +
    xlab("Hours after injection") +
    facet_wrap(~factor(label))
  return(p)
}
p_avg <- plotModuleExpressionShapes(.ME_object = MEs, .MEs_or_AvgExpr = "AvgExpr")
p_me <- plotModuleExpressionShapes(.ME_object = MEs, .MEs_or_AvgExpr = "MEs")
ggarrange(p_avg, p_me, nrow = 1, ncol = 2, legend = "none")

### PCA based on AvgExpr
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
pca <- getPCA(MEs, .MEs_or_AvgExpr = "AvgExpr")
plotPCAScree <- function(.PCA) {
  plot(x = c(1:length(.PCA$sdev)),
       y = .PCA$sdev/sum(.PCA$sdev),
       xlab = "PC", ylab = "% var explained")
}
plotPCAScree(pca) # 4 PCs
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

# PC1 and PC2
plotPCs(pca)
# PC2 and PC3
plotPCs(pca, .xPC = 2, .yPC = 3)

### Heatmap of Euclidean distance in the first X PCs, weighted by variance explained
# Helper function to get distance between two clusters
getDist <- function(.clust1, .clust2, .pca, .nPCs = 4, .weighted = TRUE) {
  dist_vec <- map(c(1:.nPCs), \(comp) {
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
getDist(.clust1 = "AE1", .clust2 = "AE10", .pca = pca) # should be small
getDist(.clust1 = "AE1", .clust2 = "AE2", .pca = pca) # should be large

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
plotmat <- getDistMat(.PCA = pca, .nPCs = 3)
heatmap(plotmat, symm = TRUE, col = topo.colors(256))
# This gives us a sense of the most similar module groups
# Big groups:
# 13, 4, 19, 8 (long early spike by eye)
# 17, 9 (flat by eye)
# 1, 7, 10 (down by eye)
# 2, 3, 15 (up by eye)
# 16, 11 (weird by eye)
# 12, 18
# 5, 14
# 6
# Small, more similar groups:
# 19, 8
# 9, 17
# 1, 10
# 11, 16

#### Clustering A and B replicates separately ####

# TODO: Use repsplit code from below to do initial clustering, 
# then exploratory pca code above to cluster all the module avgExpr
# and manually group modules

colnamesA <- grep("A", colnames(counts), value = TRUE)
colnamesB <- grep("B", colnames(counts), value = TRUE)
setequal(parse_number(colnamesA), parse_number(colnamesB)) # should be equal

# A-B Partition
countsA <- counts[, c("gene_name", colnamesA)]
countsB <- counts[, c("gene_name", colnamesB)]

# clustering
clustA <- corCluster(.cts = countsA, .tree_too = TRUE)
clustB <- corCluster(.cts = countsB, .tree_too = TRUE)

# Use matchLabels to make them as similar as possible
# before matchLabels
table(clustA$df$label, clustB$df$label)

# matching
matchedB <- matchLabels(source = clustB$df$label,
                        reference = clustA$df$label)
# after matching
table(clustA$df$label, matchedB)

#### PCA to merge similar modules ####
MEs_A <- moduleEigengenes(expr = t(countsA[,-1]),
                          colors = clustA$df$label, nPC = 1) # excludeGrey = TRUE removes unclustered (0) MEs, but we need them for our method
MEs_B <- moduleEigengenes(expr = t(countsB[,-1]),
                          colors = clustB$df$label, nPC = 1)
# comparing MEs and avgExpr
# A
p_meA <- plotModuleExpressionShapes(MEs_A, .MEs_or_AvgExpr = "MEs")
p_avgA <- plotModuleExpressionShapes(MEs_A, .MEs_or_AvgExpr = "AvgExpr")
ggarrange(p_meA, p_avgA, nrow = 1, ncol = 2, legend = "none")
# B
p_meB <- plotModuleExpressionShapes(MEs_B, .MEs_or_AvgExpr = "MEs")
p_avgB <- plotModuleExpressionShapes(MEs_B, .MEs_or_AvgExpr = "AvgExpr")
ggarrange(p_meB, p_avgB, nrow = 1, ncol = 2, legend = "none")
# Conclusion: MEs and AvgExpr still very similar.
# A and B avgExpr:
ggarrange(p_avgA, p_avgB, nrow = 1, ncol = 2, legend = "none")
# Conclusion: most clusters that have the same number in A and B
# are not the same expression shape (except 1, which makes sense)
#
# slow decreasers: 1A-1B
# early downspike: 9A-14B, 7A-9B, 6B
# early upspike: 5A-11B
# kind of early upspike: 4A-3B-5B-2A
# late upspike: 10A-11A-3A-4B-8B

# combining AB clustering results prior to PCA
AvgExpr_A <- MEs_A$averageExpr
rownames(AvgExpr_A) <- parse_number(rownames(AvgExpr_A))
colnames(AvgExpr_A) <- paste0(colnames(AvgExpr_A), "_A")
AvgExpr_B <- MEs_B$averageExpr
rownames(AvgExpr_B) <- parse_number(rownames(AvgExpr_B))
colnames(AvgExpr_B) <- paste0(colnames(AvgExpr_B), "_B")
AvgExpr_AB <- cbind(AvgExpr_A, AvgExpr_B)
pca_AB <- getPCA(.ME_object = list(averageExpr = AvgExpr_AB),
                 .MEs_or_AvgExpr = "AvgExpr")
plotPCAScree(pca_AB) # first 4 PCs
plotPCs(pca_AB, .xPC = 1, .yPC = 2)
plotPCs(pca_AB, .xPC = 2, .yPC = 3)
plotPCs(pca_AB, .xPC = 3, .yPC = 4)

# Heatmap
plotmat <- getDistMat(pca_AB, .nPCs = 4)
rownames(plotmat) <- rownames(pca_AB$x)
colnames(plotmat) <- rownames(pca_AB$x)
heatmap(plotmat, symm = TRUE, col = topo.colors(256),
        cexRow = 0.4, cexCol = 0.4)
# By eye, groups to merge:
# new cluster 1: 1A, 1B, 2B, 8A
# new cluster 2: 2A, 5B
# new cluster 3: 6A, 13A, 7B, 
# new cluster 4: 14A, 10B, 11B
# new cluster 5: 15A, 18B
# new cluster 6: 17A, 13B, 0A, 5A, 16B
# new cluster 7: 17B, 11A
# new cluster 8: 9B, 15B, 4A, 3B, 7A
# new cluster 9: 16A, 14B,
# new cluster 10: 12A, 12B
# new cluster 11: 10A, 4B, 8B, 9A, 6B, 3A
new_cluster_lookup <- list("1" = c("1A", "1B", "2B", "8A"),
                           "2" = c("2A", "5B"),
                           "3" = c("6A", "13A", "7B"),
                           "4" = c("14A", "10B", "11B"),
                           "5" = c("15A", "18B"),
                           "6" = c("17A", "13B", "0A", "5A", "16B"),
                           "7" = c("17B", "11A"),
                           "8" = c("9B", "15B", "4A", "3B", "7A"),
                           "9" = c("16A", "14B", "12A", "12B"),
                           "10" = c("10A", "4B", "8B", "9A", "6B", "3A"))
# make sure every old cluster is represented exactly once (There is no 0B or 18A):
table(unlist(new_cluster_lookup)) |> sort()
getNewLabel <- function(.old_label, .lookup, .A_or_B) {
  .old_label <- paste0(.old_label, .A_or_B)
  new_label <- lapply(new_cluster_lookup, \(x) .old_label %in% x) |> 
    unlist() |> which(useNames = TRUE) |> names()
  return(new_label)
}
# tests for getNewLabel
getNewLabel(.old_label = 6, .lookup = new_cluster_lookup,
            .A_or_B = "A") # should be 3
getNewLabel(.old_label = "6", .lookup = new_cluster_lookup,
            .A_or_B = "A")

# merging AB clustering
# A
clustdfA <- clustA$df
table(clustdfA$label) # before merging
clustdfA$label <- map(clustdfA$label,
                      \(x) getNewLabel(.old_label = as.numeric(x), 
                                       .lookup = new_cluster_lookup,
                                       .A_or_B = "A")) |> 
  unlist()
table(clustdfA$label) # after merging
# B
clustdfB <- clustB$df
table(clustdfB$label) # before merging
clustdfB$label <- map(clustdfB$label,
                      \(x) getNewLabel(.old_label = as.numeric(x), 
                             .lookup = new_cluster_lookup,
                             .A_or_B = "B")) |> 
  unlist()
table(clustdfB$label) # after merging
# we assume all genes are in the same order in both splits:
sum(clustdfA$gene_name == clustdfB$gene_name)/nrow(clustdfA) # should be 1
sum(clustdfA$label == clustdfB$label)/nrow(clustdfA) # our number of robust clusterers
# combining
clustdf <- select(clustdfA, gene_name)
clustdf$label <- map2(clustdfA$label, clustdfB$label, 
                      .f = \(a, b) {
                        if (a == b) {
                          return(a)
                        }
                        else {
                          return(paste(a, b, sep = "_"))
                        }
                      }) |> 
  unlist()
clustdf$robust <- map(clustdf$label, 
                      .f = \(x) !grepl("_", x)) |> 
  unlist()
sum(clustdf$robust)/nrow(clustdf) # should be same percent as above
clustdf_AB <- clustdf

#### Repeating with random split to compare ####

# TODO: see if you get more robust clusterers with a random split

#### Merging similar modules ####
# Issue: WGCNA doesn't have a way to post-hoc merge small modules into their closest bigger module
is_robust <- counts$gene_name %in% robust_cluster_genes
pca <- moduleEigengenes(expr = t(counts[is_robust,-1]),
                        colors = clusters$df$label[is_robust], nPC = 2)
# Further issue: WGCNA doesn't output more than PC1, so I'll need to make my own
# function to obtain PC2. But also it might not be what I'm looking for---I think I'm 
# more looking for covariance between average module expression

# TODO: genes to monitor during merging close modules part
"FBgn0004396" # CrebA (non-robust 16)
"FBgn0000250" # Cactus (robust 16)
"FBgn0260632" # Dorsal (robust 13, which looks an awful lot like cluster 16)

# Merging close modules
newlabels <- mergeCloseModules(exprData = t(counts[,-1]), 
                               colors = clusters$df$label, 
                               cutHeight = 0.4) # cutHeight is 1-cor that is our threshold for merging, so cutHeight of 0.2 causes modules to be merged if their eigengenes have >= 0.8 cor
plotDendroAndColors(dendro = clusters$tree, 
                    colors = cbind(clusters$df$label, newlabels$colors),
                    groupLabels = c("unmerged", "merged"),
                    dendroLabels = FALSE)
# checking if Dorsal and Cactus are put together
cactus_idx <- which(counts$gene_name == "FBgn0000250") # Cactus
dorsal_idx <- which(counts$gene_name == "FBgn0260632") # Dorsal
# not together in pre-merge:
clusters$df$label[cactus_idx] |> as.numeric()
clusters$df$label[dorsal_idx] |> as.numeric()
# are they together in post-merge?
newlabels$colors[cactus_idx] |> as.numeric()
newlabels$colors[dorsal_idx] |> as.numeric() # Even when a lot of merging is happening, they don't get put together. Perhaps they are not as similar expression shapes as they look

#### QC: how often do genes from the two replicates cluster together? ####
# A-B Partition
countsA <- counts[, c("gene_name", colnamesA)]
countsB <- counts[, c("gene_name", colnamesB)]

# clustering
clustA <- corCluster(.cts = countsA, .tree_too = TRUE)
clustB <- corCluster(.cts = countsB, .tree_too = TRUE)

# Use matchLabels to make them as similar as possible
# before matchLabels
table(clustA$df$label, clustB$df$label)

# matching
matchedB <- matchLabels(source = clustB$df$label,
                        reference = clustA$df$label)
# after matching
table(clustA$df$label, matchedB)

# taking the genes that cluster the same in both replicates
#robust_cluster_idxs_random <- which(clustA$df$label == matchedB)
robust_cluster_idxs_AB <- which(clustA$df$label == matchedB)

# Random Replicate Partition
# splitting into two halves, randomly partitioning counts from A or B at each timepoint
sample_idxs <- sample(c(1:length(colnamesA)), 
                      size = length(colnamesA)/2,
                      replace = FALSE) |> sort(decreasing = FALSE)
other_idxs <- setdiff(c(1:length(colnamesA)), sample_idxs)
countsA <- counts[, c("gene_name", colnamesA[sample_idxs], colnamesB[other_idxs])]
countsB <- counts[, c("gene_name", colnamesB[sample_idxs], colnamesA[other_idxs])]

# clustering
clustA <- corCluster(.cts = countsA, .tree_too = TRUE)
clustB <- corCluster(.cts = countsB, .tree_too = TRUE)

# Use matchLabels to make them as similar as possible
# before matchLabels
table(clustA$df$label, clustB$df$label)

# matching
matchedB <- matchLabels(source = clustB$df$label,
                        reference = clustA$df$label)
# after matching
table(clustA$df$label, matchedB)

robust_cluster_idxs_random <- which(clustA$df$label == matchedB)

length(intersect(robust_cluster_idxs_random, robust_cluster_idxs_AB))
length(c(robust_cluster_idxs_random, robust_cluster_idxs_AB))

# our core set of genes that cluster together in both the A-B split and in the random replicate partition
robust_cluster_genes <- counts$gene_name[intersect(robust_cluster_idxs_random, robust_cluster_idxs_AB)]
# slight differences in exact number of genes, but typically about 2000 genes

#### Saving ####
save(robust_cluster_genes, clusters, file = "data/Clustering.RData")

###################### Archive ########################
#### manually setting nClust (k) ####
# # gap statistic to determine appropriate number of clusters
# calculateGapStat <- function(.cts) {
#   # checking if any genes have entirely NAs for counts (it's ok if they have some)
#   if (sum(apply(.cts, 1, \(x) {all(is.na(x))})) != 0) { # if any rows (genes) have all NA values, they will cause cor to fail below
#     cat("NA genes in counts matrix, returning counts matrix only\n")
#     return(.cts)
#   }
#   # converting counts df to a matrix
#   cts_mat <- .cts[,-1] |> as.matrix()
#   rownames(cts_mat) <- unlist(.cts[,1]) |> as.character()
#   # calculating gap statistic using clusGap
#   cluster_fun <- function(x, k) {
#     list(cluster = cutree(hclust(d = as.dist(-cor(t(x), use = "pairwise.complete.obs")),
#                                  method = "average"),
#                           k = k))}
#   gap_stat <- clusGap(cts_mat, FUNcluster = cluster_fun, K.max = 30, B = 5)
#   return(gap_stat)
# }
# gap_test <- calculateGapStat(counts)
# optimal_k <- maxSE(gap_test$Tab[, "gap"], 
#                    gap_test$Tab[, "SE.sim"], 
#                    method = "globalSEmax")
# # Depending on what method you use, optimal k ranges from 1 to 2 to 8 to 
# # 23, so definitely will need to research more what is going on here
# optimal_k
