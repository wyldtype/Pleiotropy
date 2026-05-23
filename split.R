sapply(c("tidyr", "dplyr", "ggplot2", "ggpubr", "ggalluvial", "purrr", "readxl", "readr", "stringr", "plotly", "lmtest", "pheatmap"), FUN = require,
       character.only = TRUE)

# loading counts and functions
source("utils.r")
source("global.r")
load("data/CountsList.RData")

# loading snr values and coefficients
load("data/SignalToNoise.RData")
rm(SNR_thresh)
load("data/VariabilityThresholds.RData")
# adding variability categories to genedf
genedf$variability <- map2(genedf$signal_to_noise_dev, genedf$signal_to_noise_immune,
                           \(snr_d, snr_i) {
                             if_else(snr_d <= SNR_thresh_dev,
                                     true = if_else(snr_i <= SNR_thresh_immune,
                                                    true = "var_neither",
                                                    false = "var_immune"),
                                     false = if_else(snr_i <= SNR_thresh_immune,
                                                     true = "var_dev",
                                                     false = "var_both"))
                           }) |> unlist() 

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
plotdf <- genedf |> 
  # # Un-comment below lines to look at conservative/embryonic set:
  # left_join(select(williams, FlyBaseID, Plei_allImmune_embDev), 
  #           by = c("gene_name"="FlyBaseID")) |> 
  # select(-williams_category) |> 
  # dplyr::rename("williams_category"="Plei_allImmune_embDev") |> 
  # filter(williams_category %in% c("Developmental_Non_Pleiotropic",
  #                                                   "Immune_Non_Pleiotropic",
  #                                                   "Pleiotropic")) |> 
  left_join(y = nEnhancerdf |> filter(condition %in% c("dev", "s2_constit", "HKSM_only")) |> 
              pivot_wider(id_cols = "FlyBaseID", 
                          names_from = "condition", 
                          values_from = "nEnhancers",
                          names_prefix = "nEnhancers_"),
            by = c("gene_name"="FlyBaseID")) |> 
  mutate(nEnhancers_dev = coalesce(nEnhancers_dev, 0),
         nEnhancers_HKSM_only = coalesce(nEnhancers_HKSM_only, 0),
         nEnhancers_s2_constit = coalesce(nEnhancers_s2_constit, 0)) |> 
  left_join(y = HOGdf |> select(Dmel_FBgnID, OG_cons, nOG, nOG_dmel),
            by = c("gene_name"="Dmel_FBgnID")) |> 
  mutate(gene_family_expanded = nOG_dmel < nOG,
         hasEnhancers = nEnhancers_dev > 0 | nEnhancers_HKSM_only > 0 | nEnhancers_s2_constit,
         hasEnhancers_dev = nEnhancers_dev > 0,
         hasEnhancers_HKSM_only = nEnhancers_HKSM_only > 0,
         hasEnhancers_s2_constit = nEnhancers_s2_constit > 0,
         williams_code = dense_rank(williams_category), # 1 = dev, 2 = immune, 3 = plei
         variability_code = dense_rank(variability)) |> # 1 = var_both, 2 = var_dev, 3 = var_immune, 4 = var_neither
  filter(!is.na(gene_family_expanded)) |> 
  arrange(hasEnhancers, gene_family_expanded, williams_code)

# heatmap matrix
plot_mat <- plotdf |> select(variability_code,
                             gene_family_expanded,
                             hasEnhancers) |> 
  as.matrix()
rownames(plot_mat) <- plotdf$gene_name

# pheatmap, default doesn't scale, I have it scaling by column
pheatmap(plot_mat*1, cluster_cols = FALSE, cluster_rows = TRUE, 
         show_rownames = FALSE, scale = "none", 
         #color = c("grey80", "#D7191C"))
         color = c("grey80", "#D7191C", "#FDAE61", "#ABDDA4", "#2B83BA"))

# pretty wild un-enrichment of enhancers among gene family expansions
# Is that significant unenrichment?
plotdf |> group_by(hasEnhancers, 
                   gene_family_expanded) |> 
  summarise(n = n()) # input that to fisher test, couldn't be less significant

#### do gene family expansions have higher SNR? ####

### SNR, dev
# anova, yes gene family expansion
datadf <- filter(plotdf, gene_family_expanded & signal_to_noise_dev > 0) |> 
  mutate(isDev = williams_category == "Developmental_Non_Pleiotropic",
         isImmune = williams_category == "Immune_Non_Pleiotropic",
         isOther = williams_category == "Other",
         isPlei = williams_category == "Pleiotropic")
mod <- glm(signal_to_noise_dev ~ isPlei, #+ isDev + isImmune + isOther, 
                  data = datadf,
                  family = Gamma(link = "log"))
anova(mod)
# anova, no gene family expansion
datadf <- filter(plotdf, !gene_family_expanded & signal_to_noise_dev > 0) |> 
  mutate(isDev = williams_category == "Developmental_Non_Pleiotropic",
         isImmune = williams_category == "Immune_Non_Pleiotropic",
         isOther = williams_category == "Other",
         isPlei = williams_category == "Pleiotropic")
mod <- glm(signal_to_noise_dev ~ isPlei, #+ isDev + isImmune + isOther, 
           data = datadf,
           family = Gamma(link = "log"))
anova(mod)
# boxplot
p_devSNR <- ggplot(filter(plotdf, signal_to_noise_dev > 0), 
                   aes(x = gene_family_expanded, 
                       y = log2(signal_to_noise_dev))) +
  geom_boxplot(aes(fill = factor(williams_category,
                                 levels = c("Other", 
                                            "Developmental_Non_Pleiotropic",
                                            "Immune_Non_Pleiotropic",
                                            "Pleiotropic")))) +
  scale_fill_manual(values = ColorValues) +
  xlab("Gene family expanded outside of D. mel") +
  ylab("Signal to noise ratio\nembryonic development (log2)") +
  theme_classic() +
  theme(legend.position = "none")
p_devSNR

### SNR, immune
# anova, yes gene family expansion
datadf <- filter(plotdf, gene_family_expanded & signal_to_noise_immune > 0) |> 
  mutate(isDev = williams_category == "Developmental_Non_Pleiotropic",
         isImmune = williams_category == "Immune_Non_Pleiotropic",
         isOther = williams_category == "Other",
         isPlei = williams_category == "Pleiotropic")
mod <- glm(signal_to_noise_immune ~ isPlei, #+ isDev + isImmune + isOther, 
           data = datadf,
           family = Gamma(link = "log"))
anova(mod) # should be significant
# anova, no gene family expansion
datadf <- filter(plotdf, !gene_family_expanded & signal_to_noise_immune > 0) |> 
  mutate(isDev = williams_category == "Developmental_Non_Pleiotropic",
         isImmune = williams_category == "Immune_Non_Pleiotropic",
         isOther = williams_category == "Other",
         isPlei = williams_category == "Pleiotropic")
mod <- glm(signal_to_noise_immune ~ isPlei, #+ isDev + isImmune + isOther, 
           data = datadf,
           family = Gamma(link = "log"))
anova(mod) # shouldn't be significant
# boxplot
p_immuneSNR <- ggplot(filter(plotdf, signal_to_noise_immune > 0), 
                   aes(x = gene_family_expanded, 
                       y = log2(signal_to_noise_immune))) +
  geom_boxplot(aes(fill = factor(williams_category,
                                 levels = c("Other", 
                                            "Developmental_Non_Pleiotropic",
                                            "Immune_Non_Pleiotropic",
                                            "Pleiotropic")))) +
  scale_fill_manual(values = ColorValues) +
  xlab("Gene family expanded outside of D. mel") +
  ylab("Signal to noise ratio\nImd challenge (log2)") +
  theme_classic() +
  theme(legend.position = "none")
p_immuneSNR

ggarrange(p_devSNR, p_immuneSNR)
ggsave(filename = "figures/geneFamilyExpansions.png", 
       width = 7, height = 3)

### Repeat with embryonic/conservative set
plotdf <- genedf |> 
  # Un-comment below lines to look at conservative/embryonic set:
  left_join(select(williams, FlyBaseID, Plei_allImmune_embDev),
            by = c("gene_name"="FlyBaseID")) |>
  select(-williams_category) |>
  dplyr::rename("williams_category"="Plei_allImmune_embDev") |>
  filter(williams_category %in% c("Developmental_Non_Pleiotropic",
                                                    "Immune_Non_Pleiotropic",
                                                    "Pleiotropic")) |>
  left_join(y = nEnhancerdf |> filter(condition %in% c("dev", "s2_constit", "HKSM_only")) |> 
              pivot_wider(id_cols = "FlyBaseID", 
                          names_from = "condition", 
                          values_from = "nEnhancers",
                          names_prefix = "nEnhancers_"),
            by = c("gene_name"="FlyBaseID")) |> 
  mutate(nEnhancers_dev = coalesce(nEnhancers_dev, 0),
         nEnhancers_HKSM_only = coalesce(nEnhancers_HKSM_only, 0),
         nEnhancers_s2_constit = coalesce(nEnhancers_s2_constit, 0)) |> 
  left_join(y = HOGdf |> select(Dmel_FBgnID, OG_cons, nOG, nOG_dmel),
            by = c("gene_name"="Dmel_FBgnID")) |> 
  mutate(gene_family_expanded = nOG_dmel < nOG,
         hasEnhancers = nEnhancers_dev > 0 | nEnhancers_HKSM_only > 0 | nEnhancers_s2_constit,
         hasEnhancers_dev = nEnhancers_dev > 0,
         hasEnhancers_HKSM_only = nEnhancers_HKSM_only > 0,
         hasEnhancers_s2_constit = nEnhancers_s2_constit > 0,
         williams_code = dense_rank(williams_category), # 1 = dev, 2 = immune, 3 = plei
         variability_code = dense_rank(variability)) |> # 1 = var_both, 2 = var_dev, 3 = var_immune, 4 = var_neither
  filter(!is.na(gene_family_expanded))

### SNR, dev
# just boxplots for now
p_devSNR <- ggplot(filter(plotdf, signal_to_noise_dev > 0), 
                   aes(x = gene_family_expanded, 
                       y = log2(signal_to_noise_dev))) +
  geom_boxplot(aes(fill = factor(williams_category,
                                 levels = c("Other", 
                                            "Developmental_Non_Pleiotropic",
                                            "Immune_Non_Pleiotropic",
                                            "Pleiotropic")))) +
  scale_fill_manual(values = ColorValues) +
  xlab("Gene family expanded outside of D. mel") +
  ylab("Signal to noise ratio\nembryonic development (log2)") +
  theme_classic() +
  theme(legend.position = "none")
p_devSNR # same direction as full data

# boxplot
p_immuneSNR <- ggplot(filter(plotdf, signal_to_noise_immune > 0), 
                      aes(x = gene_family_expanded, 
                          y = log2(signal_to_noise_immune))) +
  geom_boxplot(aes(fill = factor(williams_category,
                                 levels = c("Other", 
                                            "Developmental_Non_Pleiotropic",
                                            "Immune_Non_Pleiotropic",
                                            "Pleiotropic")))) +
  scale_fill_manual(values = ColorValues) +
  xlab("Gene family expanded outside of D. mel") +
  ylab("Signal to noise ratio\nImd challenge (log2)") +
  theme_classic() +
  theme(legend.position = "none")
p_immuneSNR # somewhat same direction as full data

### Means just for fun:
### mean, dev
# anova, yes gene family expansion
datadf <- filter(plotdf, gene_family_expanded & mean_dev > 2) |> 
  mutate(isDev = williams_category == "Developmental_Non_Pleiotropic",
         isImmune = williams_category == "Immune_Non_Pleiotropic",
         isOther = williams_category == "Other",
         isPlei = williams_category == "Pleiotropic")
mod <- glm(log2(mean_dev) ~ isPlei, 
           data = datadf,
           family = Gamma(link = "log"))
anova(mod)
# anova, no gene family expansion
datadf <- filter(plotdf, !gene_family_expanded & mean_dev > 2) |> 
  mutate(isDev = williams_category == "Developmental_Non_Pleiotropic",
         isImmune = williams_category == "Immune_Non_Pleiotropic",
         isOther = williams_category == "Other",
         isPlei = williams_category == "Pleiotropic")
mod <- glm(mean_dev ~ isPlei, 
           data = datadf,
           family = Gamma(link = "log"))
anova(mod)
# boxplot
p_devMean <- ggplot(filter(plotdf, mean_dev > 2), 
                    aes(x = gene_family_expanded, 
                        y = log2(mean_dev))) +
  geom_boxplot(aes(fill = factor(williams_category,
                                 levels = c("Other", 
                                            "Developmental_Non_Pleiotropic",
                                            "Immune_Non_Pleiotropic",
                                            "Pleiotropic")))) +
  scale_fill_manual(values = ColorValues) +
  xlab("Gene family expanded outside of D. mel") +
  ylab("Mean expression level\nembryonic development (log2)") +
  theme_classic() +
  theme(legend.position = "none")
p_devMean

### mean, immune
# anova, yes gene family expansion
datadf <- filter(plotdf, gene_family_expanded & mean_immune > 2) |> 
  mutate(isDev = williams_category == "Developmental_Non_Pleiotropic",
         isImmune = williams_category == "Immune_Non_Pleiotropic",
         isOther = williams_category == "Other",
         isPlei = williams_category == "Pleiotropic")
mod <- glm(log2(mean_immune) ~ isPlei + isDev + isImmune + isOther, 
           data = datadf,
           family = Gamma(link = "log"))
anova(mod)
# anova, no gene family expansion
datadf <- filter(plotdf, !gene_family_expanded & mean_immune > 2) |> 
  mutate(isDev = williams_category == "Developmental_Non_Pleiotropic",
         isImmune = williams_category == "Immune_Non_Pleiotropic",
         isOther = williams_category == "Other",
         isPlei = williams_category == "Pleiotropic")
mod <- glm(mean_immune ~ isPlei + isDev + isImmune + isOther, 
           data = datadf,
           family = Gamma(link = "log"))
anova(mod)
# boxplot
p_immuneMean <- ggplot(filter(plotdf, mean_immune > 2), 
                       aes(x = gene_family_expanded, 
                           y = log2(mean_immune))) +
  geom_boxplot(aes(fill = factor(williams_category,
                                 levels = c("Other", 
                                            "Developmental_Non_Pleiotropic",
                                            "Immune_Non_Pleiotropic",
                                            "Pleiotropic")))) +
  scale_fill_manual(values = ColorValues) +
  xlab("Gene family expanded outside of D. mel") +
  ylab("Mean expression level\nImd challenge (log2)") +
  theme_classic() +
  theme(legend.position = "none")
p_immuneMean