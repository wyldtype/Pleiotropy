sapply(c("tidyr", "dplyr", "ggplot2", "ggpubr", "purrr", "readxl", "readr", "stringr"), FUN = require,
       character.only = TRUE)
load("data/ClusteringDev.RData")
williams <- read_xlsx(path = "data/17Jan2023_Supplemental_File_1.xlsx",
                      sheet = 2, na = "NA") 

#### Pie charts of cluster membership for Immune, Developmental, and Pleiotropic genes ####

clustdf <- left_join(x = clustdf,
                     y = select(williams, FlyBaseID, 
                                Plei_immuneResponse_embDev,
                                Plei_allImmune_allDev),
                     by = c("gene_name"="FlyBaseID")) |> 
  mutate(williams_category = if_else(is.na(Plei_allImmune_allDev),
                                     true = "Other",
                                     false = Plei_allImmune_allDev))
clustdf$robust <- clustdf$gene_name %in% robust_cluster_genes

# plotting cluster membership
# minClusterSize <- 7 # I initially was filtering all clusters with less than 11 genes but turns out Dorsal is in a cluster with 7 genes
plei_clusters <- clustdf |> filter(robust & Plei_allImmune_allDev == "Pleiotropic") |> 
  select(label) |> unique()
plotdf <- clustdf |> filter(robust) |> 
  select(label, williams_category)
plotdf_counts <- plotdf |> group_by(label) |> 
  summarise(clusterSize = n())
plotdf <- left_join(plotdf, plotdf_counts, by = "label") #|> filter(clusterSize >= minClusterSize)
table(plotdf$williams_category) # to get counts
p <- ggplot(plotdf, aes(x = factor(1), fill = factor(label))) + 
  geom_bar(position = "fill") +
  coord_polar(theta = "y") +
  facet_wrap(~factor(williams_category,
                     levels = c("Developmental_Non_Pleiotropic", "Immune_Non_Pleiotropic", "Pleiotropic", "Other"),
                     labels = c("Developmental (n = 1426)", "Immune (n = 206)", "Pleiotropic (n = 202)", "Other (n = 4394)"))) +
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ylab("") + xlab("") +
  scale_fill_brewer(palette = "Set1") +
  labs(fill = "Cluster")
p
pdf(file = "figures/pies_dev.pdf",
    width = 5, height = 5)
p
dev.off()

#### Cluster expression shapes ####

# pairing cluster assignments to expression data
avgdf <- clustdf |> filter(gene_name %in% robust_cluster_genes) |> 
  left_join(counts, by = "gene_name") |> 
  pivot_longer(cols = colnames(counts)[-1],
               names_to = "sample_id",
               values_to = "expr")
# taking average expression of each cluster in each sample
avgdf <- avgdf |> group_by(label, sample_id) |> 
  summarise(avg_expr = mean(expr))
# taking average of each cluster among replicates of the same hour
# Taking sample means before the hour mean to speed up computation 
# (each sample mean came from the same number of genes, 
# so taking the mean of the 3-4 sample means is
# equivalent to taking one mean of all genes in all replicates)
avgdf$replicate <- getReplicate(avgdf$sample_id)
avgdf$hour <- getHour(avgdf$sample_id)
# transcripts per million:
p <- ggplot(avgdf, 
       aes(x = factor(hour), y = avg_expr)) +
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
  ylab("Average cluster expression (tpm)") +
  xlab("Hours after injection") +
  facet_wrap(~factor(label))
p
# conclusion: no difference in mean cluster expr level here, unlike immune
pdf("figures/lines_dev.pdf",
    height = 4, width = 7)
p
dev.off()

# scaled expression:
scaledf <- avgdf |> group_by(label, replicate) |> 
  summarise(mean = mean(avg_expr),
            sd = sd(avg_expr))
avgdf <- left_join(avgdf, scaledf, by = c("label", "replicate"))
avgdf$scaled_expr <- (avgdf$avg_expr - avgdf$mean)/avgdf$sd
p <- ggplot(avgdf, 
            aes(x = factor(hour), y = scaled_expr)) +
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
  ylab("Average cluster expression\n(centered and scaled)") +
  xlab("Hours after injection") +
  facet_wrap(~factor(label))
#theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p
pdf("figures/lines_dev_scaled.pdf",
    height = 4, width = 7)
p
dev.off()

#### Do non-robust genes have less-responsive expression? ####
vardf <- counts |> select(gene_name)
vardf$mean_expr <- apply(as.matrix(counts[,-1]), 1, FUN = mean)
vardf$sd_expr <- apply(as.matrix(counts[,-1]), 1, FUN = sd)
vardf$cv <- vardf$sd_expr/vardf$mean_expr
plotdf <- clustdf |> left_join(vardf, by = "gene_name")
plot_counts <- plotdf |> group_by(robust, williams_category) |> 
  summarise(n = n())
p <- ggplot(plotdf, aes(y = cv, x = robust)) +
  geom_boxplot(aes(fill = robust)) +
  stat_compare_means(method = "t.test") +
  facet_wrap(~williams_category) +
  theme_classic() +
  ylab("coefficient of variation") +
  geom_text(data = plot_counts, aes(label = n, 
                                    x = robust,
                                    y = 3))
p # robust have more responsive expression
pdf("figures/boxes_dev.pdf", width = 7, height = 7)
p
dev.off()

#### Saving ####
# writing to table to extend annotated bib
outdf <- clustdf |> filter(robust) |> 
  select(gene_name, label, Plei_allImmune_allDev) |> 
  drop_na() |> 
  arrange(Plei_allImmune_allDev, label)
write_csv(outdf, file = "data/clusters_dev_BLANK.csv", col_names = TRUE)

# writing table where columns are pleiotropy category,
# for metascape's multiple gene list feature
outdf <- clustdf |> filter(robust) |> 
  select(gene_name, Plei_allImmune_allDev) |> 
  drop_na() |> 
  pivot_wider(names_from = "Plei_allImmune_allDev",
              values_from = "gene_name",
              values_fn = list)
outdf
# oh how my tidyverse heart breaks to put these in this format
# filling in the gene lists for non-Dev genes with blanks so they will be the same length
dev_list <- unlist(outdf$Developmental_Non_Pleiotropic)
immune_list <- unlist(outdf$Immune_Non_Pleiotropic)
immune_list <- c(immune_list, rep("", times = length(dev_list) - length(immune_list)))
plei_list <- unlist(outdf$Pleiotropic)
plei_list <- c(plei_list, rep("", times = length(dev_list) - length(plei_list)))
metascape <- bind_cols("Developmental_Non_Pleiotropic" = dev_list,
                       "Immune_Non_Pleiotropic" = immune_list,
                       "Pleiotropic" = plei_list)
write_csv(metascape, file = "data/metascape/metascape_dev.csv", col_names = TRUE)

# repeating with the non-robust genes, which presumably don't have as responsive expression
outdf <- clustdf |> filter(!robust) |> 
  select(gene_name, Plei_allImmune_allDev) |> 
  drop_na() |> 
  pivot_wider(names_from = "Plei_allImmune_allDev",
              values_from = "gene_name",
              values_fn = list)
outdf
# filling in the gene lists for non-Dev genes with blanks so they will be the same length
dev_list <- unlist(outdf$Developmental_Non_Pleiotropic)
immune_list <- unlist(outdf$Immune_Non_Pleiotropic)
immune_list <- c(immune_list, rep("", times = length(dev_list) - length(immune_list)))
plei_list <- unlist(outdf$Pleiotropic)
plei_list <- c(plei_list, rep("", times = length(dev_list) - length(plei_list)))
metascape <- bind_cols("Developmental_Non_Pleiotropic" = dev_list,
                       "Immune_Non_Pleiotropic" = immune_list,
                       "Pleiotropic" = plei_list)
write_csv(metascape, file = "data/metascape/metascape_dev_nonRobust.csv", col_names = TRUE)
