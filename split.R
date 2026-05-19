sapply(c("tidyr", "dplyr", "ggplot2", "ggpubr", "ggalluvial", "purrr", "readxl", "readr", "stringr", "plotly", "lmtest", "pheatmap"), FUN = require,
       character.only = TRUE)

# loading counts and functions
source("utils.r")
load("data/CountsList.RData")

# williams dataset of gene functions
williams <- read_xlsx(path = "data/17Jan2023_Supplemental_File_1.xlsx",
                      sheet = 2, na = "NA") |> 
  select(FlyBaseID, Plei_allImmune_allDev) 

# immune gene names (would be nice to just make this a comprehensive gene name lookup table)
drosophila_lookup <- read_delim("data/best_gene_summary_fb_2026_01.tsv", delim = "\t",
                                comment = "##")
colnames(drosophila_lookup) <- c("FlyBaseID", "Gene_Symbol", "Summary_Source", "Summary")
getSymbol <- function(.flybase) {
  drosophila_lookup$Gene_Symbol[drosophila_lookup$FlyBaseID == .flybase]
}
getSymbol("FBgn0014018")
getFlyBaseID <- function(.symbol) {
  drosophila_lookup$FlyBaseID[drosophila_lookup$Gene_Symbol == .symbol]
}
getFlyBaseID("PGRP-SB1")

# loading snr values and coefficients
load("data/SignalToNoise.RData")
SNR_thresh # decided in signal_to_noise.R

# loading enhancers
load("data/Enhancers.RData")

# loading HOG gene family info
load("data/HOG.RData")

### Can we assign pleiotropic genes to Development or Immune based on enhancers/expression/evolution?
# Our enhancer exploration has shown that pleiotropic genes have similar rates of enhancers
# being assigned to genes in both development and immune---
# 1) Similar pct of pleiotropic genes have dev enhancers compared to non-pleiotropic dev genes
# 2) Similar pct of pleiotropic genes have hksm enhancers compared to non-pleiotropic immune genes
# 3) But, the pleiotropic genes getting these dual enhancer assignments tend to be different genes.
# 4) Similar expression plasticity for plei and dev genes
# 5) Similar rates of gene family expansion for plei and dev genes

# So, based on these observations, can we split the Pleiotropic genes to be 
# "predominantly developmental" or "predominantly immune"?
# Basically any plei gene with an HKSM enhancer and/or high HKSM plasticity and/or low rates of gene family expansion
# Can go in the "immune" category, and any genes with 

#### expression/enhancer/evolution heatmap ####
# Instead of deciding which of these properties is the most discriminating, maybe it'll be obvious in
# a heatmap form..
# Really wasn't expecting to see that big an effect
plotdf <- genedf |> filter(williams_category %in% c("Developmental_Non_Pleiotropic",
                                                    "Immune_Non_Pleiotropic",
                                                    "Pleiotropic")) |> 
  select(gene_name, signal_to_noise_immune, signal_to_noise_dev, williams_category) |> 
  left_join(y = nEnhancerdf |> filter(condition %in% c("dev", "HKSM_only")) |> 
              pivot_wider(id_cols = "FlyBaseID", 
                          names_from = "condition", 
                          values_from = "nEnhancers",
                          names_prefix = "nEnhancers_"),
            by = c("gene_name"="FlyBaseID")) |> 
  mutate(nEnhancers_dev = coalesce(nEnhancers_dev, 0),
         nEnhancers_HKSM_only = coalesce(nEnhancers_HKSM_only, 0)) |> 
  left_join(y = HOGdf |> select(Dmel_FBgnID, OG_cons, nOG, nOG_dmel),
            by = c("gene_name"="Dmel_FBgnID")) |> 
  mutate(gene_family_pct_dmel = nOG_dmel/nOG) |> 
  mutate(williams_code = dense_rank(williams_category)) # 1 = dev, 2 = immune, 3 = plei
plot_mat <- plotdf |> select(williams_code,
                             signal_to_noise_immune, signal_to_noise_dev,
                             nEnhancers_dev, nEnhancers_HKSM_only, 
                             gene_family_pct_dmel) |> 
  as.matrix()
rownames(plot_mat) <- plotdf$gene_name

# pheatmap, default doesn't scale, I have it scaling by column
pheatmap(plot_mat, cluster_cols = FALSE, cluster_rows = TRUE, 
         show_rownames = FALSE, scale = "column", 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-2, 2, length.out = 101))

# base R heatmap, default scales by row
heatmap(plot_mat, Colv = NA, labRow = NA)

# TODO: Scaling by row makes little sense. But it is showing a pattern
# That I want to make sure I understand. When SNR_immune > SNR_dev, both
# values are likely below 1 (which is why gene_family_pct_dmel looks high
# for these genes). So I think that's not very interesting, provided these
# cases do tend to have gene_family_pct_dmel == 1. If they have <1, it might
# be a pattern worth looking into further. I could force the pct to be 0 or 1million
# or something to guarantee it's the highest value. But that would probably
# wash out every row's color

#### Specific gene set heatmaps ####
# Pleiotropic genes
plotdf_plei <- genedf |> filter(williams_category == "Pleiotropic") |> 
  select(gene_name, signal_to_noise_immune, signal_to_noise_dev) |> 
  left_join(y = nEnhancerdf |> filter(condition %in% c("dev", "HKSM_only")) |> 
              pivot_wider(id_cols = "FlyBaseID", 
                          names_from = "condition", 
                          values_from = "nEnhancers",
                          names_prefix = "nEnhancers_"),
            by = c("gene_name"="FlyBaseID")) |> 
  mutate(nEnhancers_dev = coalesce(nEnhancers_dev, 0),
         nEnhancers_HKSM_only = coalesce(nEnhancers_HKSM_only, 0)) |> 
  left_join(y = HOGdf |> select(Dmel_FBgnID, OG_cons, nOG, nOG_dmel),
            by = c("gene_name"="Dmel_FBgnID")) |> 
  mutate(has_gene_family_gain = nOG > nOG_dmel,
         has_gene_family_loss = if_else(is.na(OG_cons),
                                        true = NA,
                                        false = if_else(OG_cons == "dmel",
                                                        true = 1,
                                                        false = 0)))
plot_mat_plei <- plotdf_plei |> select(signal_to_noise_immune, signal_to_noise_dev,
                                       nEnhancers_dev, nEnhancers_HKSM_only, 
                                       has_gene_family_gain, has_gene_family_loss) |> 
  arrange(desc(signal_to_noise_immune)) |> 
  as.matrix()
rownames(plot_mat_plei) <- plotdf_plei$gene_name
pheatmap(plot_mat_plei[,c(1, 2, 5)], cluster_cols = FALSE, cluster_rows = TRUE, 
         show_rownames = FALSE, scale = "column", 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-2, 2, length.out = 101))
png("figures/heatmap_plei.png", width = 700, height = 900)
print(pheatmap(plot_mat_plei, cluster_cols = FALSE, cluster_rows = TRUE, 
         show_rownames = FALSE, scale = "column", 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-2, 2, length.out = 101)))
dev.off()

# what are those plei genes with gene family expansions and high dev SNR?
plotdf_plei |> filter(has_gene_family_gain) |> 
  arrange(desc(signal_to_noise_dev)) |> 
  left_join(drosophila_lookup, by = c("gene_name"="FlyBaseID")) |> 
  View()

# Dev genes
plotdf_dev <- genedf |> filter(williams_category == "Developmental_Non_Pleiotropic") |> 
  select(gene_name, signal_to_noise_immune, signal_to_noise_dev) |> 
  left_join(y = nEnhancerdf |> filter(condition %in% c("dev", "HKSM_only")) |> 
              pivot_wider(id_cols = "FlyBaseID", 
                          names_from = "condition", 
                          values_from = "nEnhancers",
                          names_prefix = "nEnhancers_"),
            by = c("gene_name"="FlyBaseID")) |> 
  mutate(nEnhancers_dev = coalesce(nEnhancers_dev, 0),
         nEnhancers_HKSM_only = coalesce(nEnhancers_HKSM_only, 0)) |> 
  left_join(y = HOGdf |> select(Dmel_FBgnID, OG_cons, nOG, nOG_dmel),
            by = c("gene_name"="Dmel_FBgnID")) |> 
  mutate(has_gene_family_gain = nOG > nOG_dmel,
         has_gene_family_loss = if_else(is.na(OG_cons),
                                        true = NA,
                                        false = if_else(OG_cons == "dmel",
                                                        true = 1,
                                                        false = 0)))
plot_mat_dev <- plotdf_dev |> select(signal_to_noise_immune, signal_to_noise_dev,
                                     nEnhancers_dev, nEnhancers_HKSM_only, 
                                     has_gene_family_gain, has_gene_family_loss) |> 
  as.matrix()
rownames(plot_mat_dev) <- plotdf_dev$gene_name
pheatmap(plot_mat_dev[,c(1, 2, 5)], cluster_cols = FALSE, cluster_rows = TRUE, 
         show_rownames = FALSE, scale = "column", 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-2, 2, length.out = 101))
png("figures/heatmap_dev.png", width = 700, height = 900)
pheatmap(plot_mat_dev, cluster_cols = FALSE, cluster_rows = TRUE, 
         show_rownames = FALSE, scale = "column", 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-2, 2, length.out = 101))
dev.off()

# Immune genes
plotdf_immune <- genedf |> filter(williams_category == "Immune_Non_Pleiotropic") |> 
  select(gene_name, signal_to_noise_immune, signal_to_noise_dev) |> 
  left_join(y = nEnhancerdf |> filter(condition %in% c("dev", "HKSM_only")) |> 
              pivot_wider(id_cols = "FlyBaseID", 
                          names_from = "condition", 
                          values_from = "nEnhancers",
                          names_prefix = "nEnhancers_"),
            by = c("gene_name"="FlyBaseID")) |> 
  mutate(nEnhancers_dev = coalesce(nEnhancers_dev, 0),
         nEnhancers_HKSM_only = coalesce(nEnhancers_HKSM_only, 0)) |> 
  left_join(y = HOGdf |> select(Dmel_FBgnID, OG_cons, nOG, nOG_dmel),
            by = c("gene_name"="Dmel_FBgnID")) |> 
  mutate(has_gene_family_gain = nOG > nOG_dmel,
         has_gene_family_loss = if_else(is.na(OG_cons),
                                        true = NA,
                                        false = if_else(OG_cons == "dmel",
                                                        true = 1,
                                                        false = 0)))
plot_mat_immune <- plotdf_immune |> select(signal_to_noise_immune, signal_to_noise_dev,
                                           nEnhancers_dev, nEnhancers_HKSM_only, 
                                           has_gene_family_gain, has_gene_family_loss) |> 
  as.matrix()
rownames(plot_mat_immune) <- plotdf_immune$gene_name
pheatmap(plot_mat_immune[,c(1, 2, 5)], cluster_cols = FALSE, cluster_rows = TRUE, 
         show_rownames = FALSE, scale = "column", 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-2, 2, length.out = 101))
png("figures/heatmap_immune.png", width = 700, height = 900)
pheatmap(plot_mat_immune, cluster_cols = FALSE, cluster_rows = TRUE, 
         show_rownames = FALSE, scale = "column", 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-2, 2, length.out = 101))
dev.off()

