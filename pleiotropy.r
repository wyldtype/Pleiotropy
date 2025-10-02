sapply(c("tidyr", "dplyr", "ggplot2", "purrr", "readxl", "readr", "stringr"), FUN = require,
       character.only = TRUE)
load("data/Clustering.RData")
williams <- read_xlsx(path = "data/17Jan2023_Supplemental_File_1.xlsx",
                      sheet = 2, na = "NA") 

#### Pie charts of cluster membership for Immune, Developmental, and Pleiotropic genes ####

clustdf <- left_join(x = clusters$df,
                     y = select(williams, FlyBaseID, 
                                Plei_immuneResponse_embDev,
                                Plei_allImmune_allDev),
                     by = c("gene_name"="FlyBaseID"))
clustdf$robust <- clustdf$gene_name %in% robust_cluster_genes

# plotting cluster membership
plotdf <- clustdf |> filter(robust) |> 
  select(label, Plei_allImmune_allDev)
plotdf_counts <- plotdf |> group_by(label) |> 
  summarise(clusterSize = n())
plotdf <- left_join(plotdf, plotdf_counts, by = "label") |> 
  filter(clusterSize > 10)
plotdf$class <- if_else(is.na(plotdf$Plei_allImmune_allDev),
                        true = "Non_Developmental_Non_Immune",
                        false = plotdf$Plei_allImmune_allDev)
table(plotdf$class) # to get counts
p <- ggplot(plotdf, aes(x = factor(1), fill = factor(label))) + 
  geom_bar(position = "fill") +
  coord_polar(theta = "y") +
  facet_wrap(~factor(class,
                     levels = c("Developmental_Non_Pleiotropic", "Immune_Non_Pleiotropic", "Pleiotropic", "Non_Developmental_Non_Immune"),
                     labels = c("Developmental (n = 252)", "Immune (n = 110)", "Pleiotropic (n = 33)", "Other (n = 2509)"))) +
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ylab("") + xlab("") +
  scale_fill_brewer(palette = "Set1") +
  labs(fill = "Cluster")

pdf(file = "figures/pies.pdf",
    width = 5, height = 5)
p
dev.off()

#### Cluster expression shapes ####

# need count data now
counts <- read_csv("data/normalized_counts_log_filt.csv")
colnames(counts) <- c("gene_name", colnames(counts)[-1])
avgdf <- clusters$df |> filter(gene_name %in% robust_cluster_genes) |> 
  left_join(counts, by = "gene_name") |> 
  pivot_longer(cols = colnames(counts)[-1],
               names_to = "timepoint",
               values_to = "expr") |> 
  group_by(label, timepoint) |> 
  summarise(avg_expr = mean(expr))
avgdf <- left_join(avgdf, plotdf_counts, by = "label") |> 
  filter(clusterSize > 10)
avgdf$replicate <- if_else(grepl(x = avgdf$timepoint, pattern = "A"),
                           true = "A", false = "B")
avgdf$timepoint <- parse_number(avgdf$timepoint)
avgdf$hours <- c(0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 20, 24, 30, 36, 42, 48, 72, 96, 120)[avgdf$timepoint]
p <- ggplot(avgdf, aes(x = factor(hours), y = avg_expr)) +
  geom_line(aes(color = factor(label), 
                group = interaction(label, replicate),
                linetype = replicate)) +
  geom_point(aes(color = factor(label), 
                 group = interaction(label, replicate)),
             size = 0.25) +
  scale_color_brewer(palette = "Set1") +
  scale_x_discrete(breaks = c(0, 2, 4, 6, 12, 24, 48, 120)) +
  labs(color = "Cluster") +
  theme_classic() +
  ylab("Average cluster expression (log2)") +
  xlab("Hours after injection") +
  facet_wrap(~factor(label))
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pdf("figures/lines.pdf",
    height = 4, width = 7)
p
dev.off()

