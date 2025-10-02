sapply(c("tidyr", "dplyr", "ggplot2", "purrr", "readxl", "readr", "stringr", "WGCNA", "cluster"), FUN = require,
       character.only = TRUE)

counts <- read_csv("data/normalized_counts_log_filt.csv")
colnames(counts) <- c("gene_name", colnames(counts)[-1])

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
  cat("cutting tree into", length(unique(tree_labels)), "clusters\n")
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

# testing on the whole dang thing b/c why not
clusters <- corCluster(.cts = counts, .tree_too = TRUE)
plot(clusters$tree, labels = FALSE) # looks like clustering! Definitely more than 5 clusters

#### QC: how often do genes from the two replicates cluster together? ####
colnamesA <- grep("A", colnames(counts), value = TRUE)
colnamesB <- grep("B", colnames(counts), value = TRUE) 
colnamesA
colnamesB # B is missing timepoint 4
colnamesA <- setdiff(colnamesA, "4A") # removing timepoint 4 from A for comparability

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
