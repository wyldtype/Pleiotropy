sapply(c("tidyr", "dplyr", "ggplot2", "ggpubr", "ggalluvial", "purrr", "readxl", "readr", "stringr"), FUN = require,
       character.only = TRUE)

williams <- read_xlsx(path = "data/17Jan2023_Supplemental_File_1.xlsx",
                      sheet = 2, na = "NA") |> 
  select(FlyBaseID, Plei_allImmune_allDev) 
# using the least-strict definition of pleiotropy,
# because their strict threshold was based on 
# expression, which we can make our own judgments on

load("data/ClusteringImmune.RData")
counts_immune <- counts
clustdf_immune <- clustdf
robust_immune_genes <- robust_cluster_genes
rm(counts, clustdf, robust_cluster_genes)

load("data/ClusteringDev.RData")
counts_dev <- counts
clustdf_dev <- clustdf
robust_dev_genes <- robust_cluster_genes
rm(counts, clustdf, robust_cluster_genes)

load("data/ClusteringOog.RData")
counts_oog <- counts
clustdf_oog <- clustdf
robust_oog_genes <- robust_cluster_genes
rm(counts, clustdf, robust_cluster_genes)

# combining Immune/dev results
# cluster dataframe
clustdf_immune$immune_or_dev <- "immune"
clustdf_immune$robust <- clustdf_immune$gene_name %in% robust_immune_genes
clustdf_dev$immune_or_dev <- "dev"
clustdf_dev$robust <- clustdf_dev$gene_name %in% robust_dev_genes
clustdf <- bind_rows(clustdf_immune, clustdf_dev) |> 
  left_join(y = select(williams, FlyBaseID, 
                       Plei_allImmune_allDev),
            by = c("gene_name"="FlyBaseID"))
clustdf$williams_category <- if_else(is.na(clustdf$Plei_allImmune_allDev),
                                     true = "Other",
                                     false = clustdf$Plei_allImmune_allDev)

# count matrix
colnames(counts_immune) <- c("gene_name", paste0(colnames(counts_immune[-1]), "_schlamp"))
colnames(counts_dev) <- c("gene_name", paste0(colnames(counts_dev[-1]), "_becker"))
sum(counts_immune$gene_name == counts_dev$gene_name)/nrow(counts_immune) # should be 1
counts <- bind_cols(counts_immune, counts_dev[,-1])

#### Parsing sample metadata ####
getHour <- function(.name, .info = infodf) {
  if (all(grepl("^[0-9+]*$", .name)))
    return(as.numeric(.name))
  if (grepl("_becker$", .name)) {
    out_vec <- purrr::map(.name, \(nm) {
      nm <- gsub("_becker$", "", nm)
      .info |> filter(colname == nm) |> 
        select(hour) |> pull()
    }) |> unlist()
  }
  if (grepl("_schlamp$", .name)) {
    out_vec <- purrr::map(.name, \(nm) {
      nm <- gsub("_becker$", "", nm)
      parse_number(nm)
    }) |> unlist()
  }
  return(out_vec)
}
# test for getHour
test_col <- colnames(counts)[sample(2:ncol(counts), size = 1)]
test_col
getHour(test_col)

getReplicate <- function(.name, .info = infodf) {
  if (grepl("_becker$", .name)) {
    out_vec <- purrr::map(.name, \(nm) {
      nm <- gsub("_becker$", "", nm)
      .info |> filter(colname == nm) |> 
        select(replicate) |> pull()
    }) |> unlist()
  }
  if (grepl("_schlamp$", .name)) {
    out_vec <- purrr::map(.name, \(nm) {
      nm <- gsub("_becker$", "", nm)
      if_else(grepl("A", nm), true = "A", false = "B")
    }) |> unlist()
  }
  return(out_vec)
}
# test for getReplicate
test_col <- colnames(counts)[sample(2:ncol(counts), size = 1)]
test_col
getReplicate(test_col)

#### Which genes are robust in both sets? ####
plotDevImmuneExpr <- function(.gene_idxs, .title,
                              .gene_names = NULL) {
  plotdf <- counts[counts$gene_name %in% .gene_idxs,] |> 
    pivot_longer(cols = colnames(counts)[-1],
                 names_to = "name", values_to = "expr")
  if (!is.null(.gene_names)) {
    namedf <- bind_cols(idxs = .gene_idxs,
                        names = .gene_names)
    plotdf$gene_name <- map(plotdf$gene_name, \(g) {
      gname <- namedf |> filter(idxs == g) |> 
        select(names) |> pull()
      return(gname[1])
    }) |> unlist()
  }
  plotdf$dataset <- if_else(grepl("_schlamp$", plotdf$name),
                            true = "Imd challenge",
                            false = "embryonic dev")
  plotdf$hour <- map(plotdf$name, .f = getHour) |> unlist()
  plotdf$replicate <- map(plotdf$name, .f = getReplicate) |> unlist()
  p <- ggplot(plotdf, aes(x = factor(hour), y = expr)) +
    geom_point(aes(color = factor(gene_name)),
               size = 0.25) +
    geom_line(aes(color = factor(gene_name), 
                  group = interaction(gene_name, replicate),
                  linetype = replicate)) +
    #scale_color_brewer(palette = "Set1") +
    labs(color = "Cluster") +
    theme_classic() +
    ylab("Expression (tpm)") +
    xlab("Timepoint") +
    facet_wrap(~factor(dataset)) +
    ggtitle(.title)
  return(p)
}
genedf <- clustdf |> pivot_wider(id_cols = c("gene_name", "williams_category"),
                                 names_from = "immune_or_dev",
                                 values_from = c("label", "robust"))
# How many are pleiotropic?
genedf |> filter(robust_immune & robust_dev) |> 
  select(williams_category) |> table()
plotdf <- filter(genedf, robust_immune & robust_dev)
plotdf$cat <- factor(plotdf$williams_category,
                     labels = c("Developmental", "Other", "Immune", "Pleiotropic"),
                     levels = c("Developmental_Non_Pleiotropic", "Other", "Immune_Non_Pleiotropic", "Pleiotropic"))
ggplot(plotdf, 
       aes(x = cat)) +
  geom_bar(aes(fill = cat)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.title = element_blank()) +
  xlab("") + ylab("Number of genes")

#### Visualizing gene expression ####
# getting list of pleiotropic immune/dev responsive genes:
plei_genes <- genedf |> filter(robust_immune & robust_dev) |> 
  filter(williams_category == "Pleiotropic") |> 
  select(gene_name) |> pull()

# 12 Genes found when using the strictest definition of pleiotropy: Plei_immuneResponse_embDev
# FBgn0000251 = Caudal, promotes constitutive AMP production in barrier epithelia
# FBgn0001291 = Jra, JNK pathway, dimer with kayak to form dAP-1
# FBgn0001311 = krotzkopf verkehrt, chitin synthase for embryonic cuticle and wound healing
# FBgn0003205 = Ras85D, oncogene
# FBgn0004864 = Hopscotch, JAK/STAT
# FBgn0010303 = Hemipterous, kinase for basket, JNK pathway
# FBgn0026404 = Dronc, apoptosis initiator caspase
# FBgn0030018 = Slipper, MAP kinase kinase kinase, JNK pathway
# FBgn0035049 = Mmp1, ECM remodeling, possibly recruited by JNK pathway
# FBgn0040068 = Vav, GEF downstream of Rac1
# FBgn0086357 = Sec61-alpha, one of three subunits of protein channel to import nascent proteins to the ER
# FBgn0265464 = Traf6 (aka Traf2), phosphorylated by Pelle after Toll activation, also associted with JNK pathway

# Toll:
cactus_idx <- "FBgn0000250" # Cactus
dorsal_idx <- "FBgn0260632" # Dorsal
# ER:
creba_idx <- "FBgn0004396" # CrebA
sec61a_idx <- "FBgn0086357" # Sec61-alpha
sec61b_idx <- "FBgn0010638" # Sec61-beta (not flagged as pleiotropic)
sec61g_idx <- "FBgn0031049" # Sec61-gamma (not flagged as pleiotropic)
# JNK pathway:
jra_idx <- "FBgn0001291" # Jra
slipper_idx <- "FBgn0030018" # Slipper
hemipterous_idx <- "FBgn0010303" # Hemipterous
# other cool ones
mmp1_idx <- "FBgn0035049"
# ones that are not variable in one dataset (controls to reference for CV check)
krotzkopf_idx <- "FBgn0001311"
ras85_idx <- "FBgn0003205"
# 3-1 vs 4-1 genes
jra_idx <- "FBgn0001291" # Jra
bom1 <- "FBgn0034329" # bomanin short 1
bom2 <- "FBgn0025583" # bomanin short 2
pepLF <- "FBgn0035977" # Peptidoglycan recognition protein LF
pepSD <- "FBgn0035806" # Peptidoglycan recognition protein SD


plotDevImmuneExpr(.gene_idxs = c(krotzkopf_idx, ras85_idx),
                  .title = "Not Responsive",
                  .gene_names = c("krotzkopf verkehrt", "Ras85D"))
plotDevImmuneExpr(.gene_idxs = c(creba_idx, sec61a_idx, sec61b_idx, sec61g_idx),
                  .title = "ER",
                  .gene_names = c("CrebA", "Sec61-alpha", "Sec61-beta", "Sec61-gamma"))
plotDevImmuneExpr(.gene_idxs = c(jra_idx, slipper_idx, hemipterous_idx),
                  .gene_names = c("Jra", "Slipper", "Hemipterous"),
                  .title = "JNK Pathway")
plotDevImmuneExpr(.gene_idxs = c(cactus_idx, dorsal_idx),
                  .gene_names = c("Cactus", "Dorsal"),
                  .title = "Toll Pathway")

plotDevImmuneExpr(.gene_idxs = c(jra_idx, bom1, bom2),
                  .gene_names = c("Jra", "Bomanin Short 1", "Bomanin Short 2"),
                  .title = "")
plotDevImmuneExpr(.gene_idxs = c(jra_idx),
                  .gene_names = c("Jra"),
                  .title = "Jra alone")

plotDevImmuneExpr(.gene_idxs = c(pepLF, pepSD),
                  .gene_names = c("Pep LF (Pleiotropic)", "Pep SD (Immune)"),
                  .title = "Peptidoglycan recognition proteins")

#### How variable is gene expression in each dataset? ####

# motivation: All the pleiotropic examples above seem to be more
# highly expressed in development, is this common for dual-robust genes?

# Every dual-robust gene, all expression data
dual_robust_idxs <- genedf |> filter(robust_immune & robust_dev) |> 
  select(gene_name) |> pull()
megadf <- counts[counts$gene_name %in% dual_robust_idxs,] |> 
  pivot_longer(cols = colnames(counts)[-1],
               names_to = "name", values_to = "expr")
megadf$dataset <- if_else(grepl("_schlamp$", megadf$name),
                          true = "Imd challenge",
                          false = "embryonic dev")
# TODO: these take a moment, seeing as they're just lookup functions,
# we could re-factor this to left_join the results one-to-many instead
megadf$hour <- map(megadf$name, .f = getHour) |> unlist()
megadf$replicate <- map(megadf$name, .f = getReplicate) |> unlist()

vardf <- megadf |> group_by(gene_name, dataset) |> 
  summarise(mean_expr = mean(expr),
            sd_expr = sd(expr)) |> 
  mutate(cv = sd_expr/mean_expr)

vardf <- left_join(vardf, williams, 
                   by = c("gene_name"="FlyBaseID"))
vardf$williams_category <- if_else(is.na(vardf$Plei_allImmune_allDev),
                                   true = "Other", false = vardf$Plei_allImmune_allDev)
# mean
ggplot(vardf, aes(x = mean_expr)) +
  geom_density(aes(fill = dataset)) +
  facet_wrap(~williams_category)
ggplot(vardf, aes(x = log2(mean_expr))) +
  geom_density(aes(fill = dataset), alpha = 0.75) +
  facet_wrap(~williams_category)
# cv
ggplot(vardf, aes(x = cv)) +
  geom_density(aes(fill = dataset))
ggplot(vardf, aes(x = cv)) +
  geom_density(aes(fill = dataset)) +
  facet_wrap(~williams_category)

# Conclusion: Immune/Dev datasets have similar mean expression,
# But dev has more variable expression

#### Alleuvial: Which cluster shapes are most common among pleiotropic, immune, and developmental genes? ####
plotdf <- clustdf |> select(gene_name, williams_category, immune_or_dev, label, robust) |> 
  pivot_wider(id_cols = c("gene_name", "williams_category"),
              names_from = immune_or_dev, values_from = c("label", "robust"))
plotdf_counts <- plotdf |> 
  filter(robust_immune & robust_dev) |> 
  group_by(label_immune, label_dev, williams_category) |> 
  summarise(ngenes = n()) |> 
  ungroup()
category_counts <- plotdf |> 
  filter(robust_immune & robust_dev) |>
  group_by(williams_category) |> 
  summarise(ncategory = n()) |> 
  ungroup()
plotdf_counts <- left_join(plotdf_counts, category_counts,
                           by = "williams_category",
                           relationship = "many-to-one")
plotdf_counts$pctcategory <- plotdf_counts$ngenes/plotdf_counts$ncategory
is_alluvia_form(plotdf_counts, silent = TRUE)
p_pct <- ggplot(plotdf_counts,
       aes(y = pctcategory, axis1 = label_immune, axis2 = label_dev)) +
  geom_alluvium(aes(fill = williams_category), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Imd challenge", "Embryonic Development"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  #ggtitle("Cluster membership in Imd challenge and embryonic development, by gene category") +
  theme_classic() +
  ylab("percent of genes,\nby gene category")
# separated by category:
p_n <- ggplot(plotdf_counts,
       aes(y = ngenes, axis1 = label_immune, axis2 = label_dev)) +
  geom_alluvium(aes(fill = williams_category), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Imd challenge", "Embryonic Development"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  #ggtitle("Cluster membership in Imd challenge and embryonic development, by gene category") +
  theme_classic() +
  ylab("number of genes")

pdf("figures/alluvial_pct.pdf",
    width = 8, height = 16)
p_pct
dev.off()

pdf("figures/alluvial_n.pdf",
    width = 8, height = 16)
p_n
dev.off()

#### Investigating alleuvial examples ####
# 3-1 (Plei) vs 4-1 (Immune)
clustdf |> filter(robust) |>
  select(gene_name, williams_category, label, immune_or_dev) |> 
  filter(williams_category %in% c("Immune_Non_Pleiotropic", "Pleiotropic")) |> 
  pivot_wider(id_cols = c("gene_name", "williams_category"),
              names_from = immune_or_dev, values_from = label) |> 
  filter(immune %in% c(3, 4) & dev == 1) |> View()
