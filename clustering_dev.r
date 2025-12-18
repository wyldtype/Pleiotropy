sapply(c("tidyr", "dplyr", "ggplot2", "ggpubr", "purrr", "readxl", "readr", 
         "stringr", "WGCNA", "cluster", "DESeq2", "circlize", "stats"), FUN = require,
       character.only = TRUE)

# loading counts, infodf, and getter functions
load("data/DevCounts.RData")

# TODO: Filtering out timepoints before 3hr AEL (pre-zygotic control)

#### Checking for 0-count genes ####

# also removing 0-count genes in split counts (which will also remove 
# 0-count genes in full count matrix)
colnamesA_idxs <- which(getReplicate(colnames(counts)[-1]) %in% c(1, 2))
colnamesB_idxs <- which(getReplicate(colnames(counts)[-1]) %in% c(3, 4))
colnamesA <- colnames(counts)[-1][colnamesA_idxs]
colnamesB <- colnames(counts)[-1][colnamesB_idxs]
colnamesA
colnamesB

# A-B Partition
countsA <- counts[, c("gene_name", colnamesA)]
countsB <- counts[, c("gene_name", colnamesB)]

# Any genes all 0s now?
sum(rowSums(countsA[,-1]) == 0)
sum(rowSums(countsB[,-1]) == 0) # if not, don't need to remove
# keep_idxs <- rowSums(countsA[,-1]) > 0 & rowSums(countsB[,-1]) > 0

# # removing
# counts <- counts[keep_idxs,]
# countsA <- countsA[keep_idxs,]
# countsB <- countsB[keep_idxs,]

# all hours are still represented in A-B split?
setdiff(x = getHour(colnames(countsA)[-1]),
        y = getHour(colnames(countsB)[-1]))
setdiff(x = getHour(colnames(countsB)[-1]),
        y = getHour(colnames(countsA)[-1])) # should be empty

#### QC: Dorsal/Cactus/CrebA ####
# Their expression is very different between rep A and B
# But it seems to be different in the same way
# maybe it's a clue for how to group genes
cactus_idx <- "FBgn0000250" # Cactus
dorsal_idx <- "FBgn0260632" # Dorsal
creba_idx <- "FBgn0004396" # CrebA
plotdf <- tibble(sample = colnames(counts)[-1],
                 hour = getHour(colnames(counts)[-1]),
                 replicate = getReplicate(colnames(counts)[-1]),
                 cactus = as.numeric(counts[counts$gene_name == cactus_idx, -1]),
                 dorsal = as.numeric(counts[counts$gene_name == dorsal_idx, -1]),
                 creba = as.numeric(counts[counts$gene_name == creba_idx, -1])) |> 
  pivot_longer(cols = c("cactus", "dorsal", "creba"), names_to = "gene", values_to = "expr")
ggplot(plotdf, aes(x = factor(hour), y = expr)) +
  geom_point(aes(color = factor(gene)),
             size = 0.25) +
  geom_line(aes(color = factor(gene), 
                group = interaction(gene, replicate),
                linetype = replicate)) +
  #scale_color_brewer(palette = "Set1") +
  labs(color = "cluster") +
  theme_classic() +
  ylab("expression (transcripts per million)") +
  xlab("developmental stage (hours)") +
  facet_wrap(~factor(gene))

#### Smoothing by spline ####
deg_free <- 7

# applying to entire count matrix
smoothCountMatrix <- function(.cts_x, .cts_y, .df = deg_free) {
  ngenes <- nrow(.cts_x)
  hour_vec_x <- getHour(colnames(.cts_x)) |> as.numeric()
  hour_vec_y <- getHour(colnames(.cts_y)) |> as.numeric()
  common_hours <- sort(union(hour_vec_x, hour_vec_y), decreasing = FALSE)
  outlist <- map(c(1:ngenes), \(g) {
    cat(g, "/", ngenes, "\n")
    # mean count between replciates, if the hour is present in both reps,
    # if not just the single count from one replicate
    gene_vec <- purrr::map(common_hours, \(hr) {
      if (hr %in% hour_vec_x & hr %in% hour_vec_y)
        return(as.numeric(.cts_x[g, hour_vec_x == hr] + 
                            .cts_y[g, hour_vec_y == hr])/2)
      if (hr %in% hour_vec_x)
        return(as.numeric(.cts_x[g, hour_vec_x == hr]))
      if (hr %in% hour_vec_y)
        return(as.numeric(.cts_y[g, hour_vec_y == hr]))
    }) |> unlist()
    smooth_y <- smooth.spline(x = common_hours,
                              y = gene_vec,
                              df = .df,
                              keep.stuff = TRUE)
    return(smooth_y$y)
  })
  if (length(outlist) == 1) {
    outmat <- purrr::reduce(outlist[[1]], .f = cbind) # I have no idea why reducing doesn't work on single-element lists vectors. Reduce is the same problem
    # giving first timepoint its unsmoothed counts back
    first_idx <- which.min(common_hours)
    mean_vec <- unlist((.cts_x[, first_idx] + .cts_y[, first_idx])/2) |> 
      as.numeric()
    outmat[, first_idx] <- mean_vec
    colnames(outmat) <- common_hours
  }
  if (length(outlist) > 1) {
    outmat <- purrr::reduce(outlist, .f = rbind)
    # giving first timepoint its unsmoothed counts back
    first_idx <- which.min(common_hours)
    mean_vec <- unlist((.cts_x[, first_idx] + .cts_y[, first_idx])/2) |> 
      as.numeric()
    outmat[, first_idx] <- mean_vec
    colnames(outmat) <- common_hours
  }
  
  return(outmat)
}

# tests for smoothCountMatrix
colnames1_idxs <- which(getReplicate(colnames(counts)[-1]) == 1) + 1
colnames2_idxs <- which(getReplicate(colnames(counts)[-1]) == 2) + 1
gene_idx <- creba_idx
test_counts_1 <- counts[counts$gene_name %in% gene_idx, 
                        colnames1_idxs, 
                        drop = FALSE]
test_counts_2 <- counts[counts$gene_name %in% gene_idx,
                        colnames2_idxs,
                       drop = FALSE]
test_smooth <- smoothCountMatrix(.cts_x = test_counts_1,
                                 .cts_y = test_counts_2)
plotdf <- bind_rows(tibble(status = "before1",
                           expr = as.numeric(test_counts_1),
                           hour = getHour(colnames(test_counts_1))),
                    tibble(status = "before2",
                           expr = as.numeric(test_counts_2),
                           hour = getHour(colnames(test_counts_2))),
                    tibble(status = "after",
                           expr = as.numeric(test_smooth),
                           hour = as.numeric(colnames(test_smooth))))
ggplot(plotdf, aes(x = hour, y = expr)) +
  geom_line(aes(group = status, 
                color = grepl("before", status),
                linetype = status))

# smoothing
deg_free
colnames1_idxs <- colnames(counts)[-1][which(getReplicate(colnames(counts)[-1]) == 1)]
colnames2_idxs <- colnames(counts)[-1][which(getReplicate(colnames(counts)[-1]) == 2)]
colnames3_idxs <- colnames(counts)[-1][which(getReplicate(colnames(counts)[-1]) == 3)]
colnames4_idxs <- colnames(counts)[-1][which(getReplicate(colnames(counts)[-1]) == 4)]
smoothA <- smoothCountMatrix(.cts_x = counts[,colnames1_idxs],
                             .cts_y = counts[,colnames2_idxs],
                             .df = deg_free)
smoothA <- bind_cols(counts$gene_name, smoothA)
colnames(smoothA)[1] <- "gene_name"
smoothB <- smoothCountMatrix(.cts_x = counts[,colnames3_idxs],
                             .cts_y = counts[,colnames4_idxs],
                             .df = deg_free)
smoothB <- bind_cols(counts$gene_name, smoothB)
colnames(smoothB)[1] <- "gene_name"

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

# testing on the whole thing to optimize cutting height
test_clusters <- corCluster(.cts = counts[sample(c(1:nrow(counts)), size = 1000),], .tree_too = TRUE)
plotDendroAndColors(dendro = test_clusters$tree, colors = test_clusters$df$label, dendroLabels = FALSE)

#### QC: A-B split ####
# clustering A and B replicates separately to decide which genes have
# robust clustering (which are presumably the ones that have dynamic and
# reproducible expression responses during Embrogenesis)
clustA <- corCluster(.cts = smoothA, .tree_too = TRUE)
plotDendroAndColors(dendro = clustA$tree, colors = clustA$df$label, dendroLabels = FALSE)
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
                                       .separate_replicates = FALSE) {
  if (.MEs_or_AvgExpr == "MEs") {
    expr_mat <- .ME_object$eigengenes
    ylabel <- "Module eigengene expression"
  }
  else {
    expr_mat <- .ME_object$averageExpr
    ylabel <- "Average cluster expression (log2)"
  }
  plotdf <- expr_mat |>
    mutate(sample_id = rownames(expr_mat)) |> 
    pivot_longer(cols = colnames(expr_mat),
                 names_to = "label",
                 values_to = "expr")
  if (.separate_replicates) {
    plotdf$replicate <- getReplicate(plotdf$sample_id)
  }
  plotdf$hour <- getHour(plotdf$sample_id)
  p <- ggplot(plotdf, aes(x = factor(hour), y = expr)) +
    geom_point(aes(color = factor(label)),
               size = 0.25) +
    #scale_color_brewer(palette = "Set1") +
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
# conclusion: most clusters seem to be capturing real expression patterns
# consistent in both replicates, but a few have very dramatic spikes at single
# timepoints in single replicates. I could delve into what genes these are
# exactly that is causing this, but for now I'm guessing the problem
# will resolve when separating out the non-robust clusterers

# A vs B clusters by avgExpr:
ggarrange(p_avgA, p_avgB, nrow = 1, ncol = 2, legend = "none")
clustA$df |> filter(gene_name %in% c(dorsal_idx, cactus_idx, creba_idx))
clustB$df |> filter(gene_name %in% c(dorsal_idx, cactus_idx, creba_idx))

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
plotPCAScree(pca_AB) # first 3 PCs
nPCs <- 3

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
nPCs <- 3
getDist(.clust1 = "AE1_A", .clust2 = "AE1_A", .pca = pca_AB, .nPCs = nPCs) # should be 0
getDist(.clust1 = "AE1_A", .clust2 = "AE3_B", .pca = pca_AB, .nPCs = nPCs) # should be small
getDist(.clust1 = "AE2_A", .clust2 = "AE1_B", .pca = pca_AB, .nPCs = nPCs) # should be small
getDist(.clust1 = "AE1_A", .clust2 = "AE2_A", .pca = pca_AB, .nPCs = nPCs) # should be large

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
cutHeight <- 5
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
                                       cutHeight = 0.4)
new_labels_A_WGCNA <- new_labels_list_A$colors
new_labels_list_B <- mergeCloseModules(exprData = t(as.matrix(smoothB[,-1])),
                                       colors = clustB$df$label,
                                       cutHeight = 0.4)
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
# a lot of overlap

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
clusters <- corCluster(.cts = counts, .tree_too = TRUE)
plotDendroAndColors(dendro = clusters$tree, colors = clusters$df$label, dendroLabels = FALSE)

### merging clusters using WGCNA::mergeCloseModules 
new_label_list <- mergeCloseModules(exprData = t(counts[,-1]),
                                    colors = clusters$df$label,
                                    cutHeight = 0.5)
length(table(new_label_list$colors))
new_labels_WGCNA <- new_label_list$colors
newMEs_WGCNA <- moduleEigengenes(expr = t(counts[,-1]), 
                                 colors = new_labels_WGCNA)

plotModuleExpressionShapes(newMEs_WGCNA, .MEs_or_AvgExpr = "AvgExpr",
                           .separate_replicates = TRUE)

### merging clusters by PCA
MEs <- moduleEigengenes(expr = t(counts[,-1]), colors = clusters$df$label)
plotModuleExpressionShapes(MEs, .MEs_or_AvgExpr = "AvgExpr",
                           .separate_replicates = TRUE) # premerge

pca <- getPCA(MEs, .MEs_or_AvgExpr = "AvgExpr")
plotPCAScree(pca) # 3 PCs
nPCs <- 3
plotPCs(pca, .xPC = 1, .yPC = 2)
plotPCs(pca, .xPC = 2, .yPC = 3)
getDist("AE4", "AE10", .pca = pca, .nPCs = nPCs, .weighted = TRUE)
getDist("AE2", "AE1", .pca = pca, .nPCs = nPCs, .weighted = TRUE)

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
newMEs_PCA <- moduleEigengenes(t(counts[,-1]), colors = new_labels_PCA)

## comparing merge results
newMEs <- moduleEigengenes(t(counts[,-1]), colors = clusters$df$label)
# pre-merge:
plotModuleExpressionShapes(newMEs, .MEs_or_AvgExpr = "AvgExpr",
                           .separate_replicates = TRUE) # postmerge
# WGCNA::mergeCloseModules:
plotModuleExpressionShapes(newMEs_WGCNA, .MEs_or_AvgExpr = "AvgExpr",
                           .separate_replicates = TRUE) # postmerge
# PCA-distance merge:
plotModuleExpressionShapes(newMEs_PCA, .MEs_or_AvgExpr = "AvgExpr",
                           .separate_replicates = TRUE) # postmerge
# conclusion: both seem to work well and come up with similar shapes

## setting new labels
clustdf <- clusters$df
table(dense_rank(new_labels_WGCNA)) # WGCNA groups modules without necessarily keeping the values 1,2,3,etc.
clustdf$label <- new_labels_WGCNA |> dense_rank()

# recalculate MEs on un-smoothed counts
MEs <- moduleEigengenes(t(counts[,-1]), colors = clustdf$label)
plotModuleExpressionShapes(MEs, .MEs_or_AvgExpr = "AvgExpr",
                           .separate_replicates = TRUE) # postmerge

# adding robust column and checking if any clusters are over-represented in robust fraction
clustdf$robust <- clustdf$gene_name %in% robust_cluster_genes
table(clustdf$robust, clustdf$label) # 3 is the least represented

# MEs on just robust genes
MEs_robust <- moduleEigengenes(t(counts[counts$gene_name %in% robust_cluster_genes,-1]), 
                        colors = clustdf$label[clustdf$robust])
plotModuleExpressionShapes(MEs_robust, .MEs_or_AvgExpr = "AvgExpr",
                           .separate_replicates = TRUE) # not very different

clustdf |> filter(gene_name %in% c(dorsal_idx, cactus_idx, creba_idx))

#### Saving ####
save(robust_cluster_genes, clustdf, counts, infodf, file = "data/ClusteringDev.RData")

########################### Archive ################################
