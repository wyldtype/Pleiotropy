sapply(c("tidyr", "dplyr", "ggplot2", "ggpubr", "ggalluvial", "purrr", "readxl", "readr", "stringr", "plotly", "lmtest", "pheatmap"), FUN = require,
       character.only = TRUE)

# loading counts and functions
source("utils.r")
source("global.r")
load("data/CountsList.RData")

# getter functions
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

#### Example genes from Williams et al. ####

### Expression of Imd components from Williams et al. figure in 
# Imd challenge, Embyonic Dev, and Oogenesis

# Our colors
immune_color <- ColorValues["Immune_Non_Pleiotropic"] |> as.character()
plei_color <- ColorValues["Pleiotropic"] |> as.character()
dev_color <- ColorValues["Developmental_Non_Pleiotropic"] |> as.character()
none_color <- ColorValues["Other"] |> as.character()

### Example expression and SNR values for BAD poster 2026
# Ideal: 1 dev, 1 immune, 2 plei genes in dev and immune environment
# dev gene is more variable in dev, immune in immune
# 1 plei more vairable in dev, 1 in immune
# filtering for example genes (change filters to select different genes):
genedf |> left_join(select(williams, FlyBaseID, Plei_allImmune_embDev),
                    by = c("gene_name"="FlyBaseID")) |> 
  filter(signal_to_noise_immune > 0.5 &
           Plei_allImmune_embDev == "Immune_Non_Pleiotropic" &
           mean_dev > 0, mean_immune > 30) |> 
  left_join(select(drosophila_lookup, FlyBaseID, Gene_Symbol, Summary),
            by = c("gene_name"="FlyBaseID")) |> 
  View()
# example genes:
imm_gene <- getFlyBaseID("glob1")
dev_gene <- getFlyBaseID("fkh")
plei_gene <- getFlyBaseID("Rap1")
# corresponding SNRs
imm_gene_title <- paste("globin 1, SNR_imd =", 
                        round(genedf$signal_to_noise_immune[genedf$gene_name == imm_gene], digits = 2),
                        "SNR_dev =", 
                        round(genedf$signal_to_noise_dev[genedf$gene_name == imm_gene], digits = 2))
dev_gene_title <- paste("forkhead, SNR_imd =", 
                        round(genedf$signal_to_noise_immune[genedf$gene_name == dev_gene], digits = 2),
                        "SNR_dev =", 
                        round(genedf$signal_to_noise_dev[genedf$gene_name == dev_gene], digits = 2))
plei_gene_title <- paste("Rap1, SNR_imd =", 
                        round(genedf$signal_to_noise_immune[genedf$gene_name == plei_gene], digits = 2),
                        "SNR_dev =", 
                        round(genedf$signal_to_noise_dev[genedf$gene_name == plei_gene], digits = 2))

# plot for poster
ggarrange(plotlist = plotExpressionProfileAllEnvironments(.gene_idxs = c(imm_gene, dev_gene, pleiImm_gene), 
                                     .gene_names = c(imm_gene_title, dev_gene_title, plei_gene_title),
                                     .environments = c("dev", "immune"),
                                     .color_vec = c(immune_color, dev_color, plei_color)),
          nrow = 2, ncol = 1, legend = "none")
ggsave(filename = "figures/example_SNRs.png", width = 10, height = 7)

### Plots to match Williams et al. pathways figure
# # Williams et al. colors
# immune_color <- "#F04D23"
# plei_color <- "#0083C3"
# none_color <- "grey60"
### Extracellular
pgrp_sd <- getFlyBaseID("PGRP-SD")
pgrp_sb1 <- getFlyBaseID("PGRP-SB1")
pgrp_lb <- getFlyBaseID("PGRP-LB")
pgrp_le <- getFlyBaseID("PGRP-LE")
plotlist <- plotExpressionProfileAllEnvironments(.gene_idxs = c(pgrp_sd, pgrp_sb1, pgrp_lb, pgrp_le), 
                                                 .gene_names = c("PGRP-SD", "PGRP-SB1", "PGRP-LB", "PGRP-LE"),
                                                 .environments = c("oog", "dev", "immune"),
                                                 .color_vec = c(immune_color, immune_color, immune_color, plei_color))
p <- ggarrange(plotlist = plotlist,
               nrow = 2, ncol = 2, legend = "none")
p
ggsave("figures/Imd_expression_extracellular.png",
       width = 6, height = 6, plot = p)

# SNR ranks
genedf |> filter(gene_name == pgrp_sd) |> select(immune_rank, dev_rank, oog_rank)
genedf |> filter(gene_name == pgrp_sb1) |> select(immune_rank, dev_rank, oog_rank)
genedf |> filter(gene_name == pgrp_lb) |> select(immune_rank, dev_rank, oog_rank)
genedf |> filter(gene_name == pgrp_le) |> select(immune_rank, dev_rank, oog_rank)

### Membrane bound
pgrp_lc <- getFlyBaseID("PGRP-LC")
pgrp_la <- getFlyBaseID("PGRP-LA")
pgrp_lf <- getFlyBaseID("PGRP-LF")
plotlist <- plotExpressionProfileAllEnvironments(.gene_idxs = c(pgrp_lc, pgrp_la, pgrp_lf), 
                                                 .gene_names = c("PGRP-LC", "PGRP-LA", "PGRP-LF"),
                                                 .environments = c("oog", "dev", "immune"),
                                                 .color_vec = c(immune_color, immune_color, plei_color))
p <- ggarrange(plotlist = plotlist,
               nrow = 3, ncol = 1, legend = "none")
p
ggsave("figures/Imd_expression_membrane.png",
       width = 5, height = 6, plot = p)
# Intracellular pinwheel complex (ubiquitin ligase activity):
uev1a <- getFlyBaseID("Uev1A")
bendless <- getFlyBaseID("ben")
effete <- getFlyBaseID("eff")
diap2 <- getFlyBaseID("Diap2")
diap1 <- getFlyBaseID("Diap1")
p <- ggarrange(plotlist = plotExpressionProfileAllEnvironments(.gene_idxs = c(uev1a, bendless, effete, diap2), 
                                                               .gene_names = c("Uev1A", "bendless", "effete", "Diap2"),
                                                               .environments = c("oog", "dev", "immune"),
                                                               .color_vec = c(immune_color, plei_color, plei_color, plei_color)),
               nrow = 2, ncol = 2, legend = "none")
p
ggsave("figures/Imd_expression_ubiquitin.png",
       width = 6, height = 6, plot = p)

# Diap2 vs Diap1 (cause Diap1 has the double enhancer situation FBgn0260635 = Diap1)
ggarrange(plotlist = plotExpressionProfileAllEnvironments(.gene_idxs = c(diap1, diap2), 
                                                          .gene_names = c("Diap1", "Diap2"),
                                                          .environments = c("oog", "dev", "immune"),
                                                          .color_vec = c(plei_color, plei_color)),
          nrow = 1, ncol = 3, legend = "none")

# Intracellular pentagram complex (all plei so not that helpful):
imd <- getFlyBaseID("imd")
dredd <- getFlyBaseID("Dredd")
fadd <- getFlyBaseID("Fadd")
p <- ggarrange(plotlist = plotExpressionProfileAllEnvironments(.gene_idxs = c(imd, dredd, fadd), 
                                                               .gene_names = c("Imd", "Dredd", "Fadd"),
                                                               .environments = c("oog", "dev", "immune"),
                                                               .color_vec = c(plei_color, plei_color, plei_color, plei_color,
                                                                              plei_color, plei_color, plei_color, plei_color)),
               nrow = 3, ncol = 1, legend = "none")
p
ggsave("figures/Imd_expression_ImdComplex.png",
       width = 5, height = 6, plot = p)

# Tak1 and Tab2
tab2 <- getFlyBaseID("Tab2")
tak1 <- getFlyBaseID("Tak1")
plotlist <- plotExpressionProfileAllEnvironments(.gene_idxs = c(tab2, tak1), 
                                                 .gene_names = c("Tab2", "Tak1"),
                                                 .environments = c("oog", "dev", "immune"),
                                                 .color_vec = c(immune_color, plei_color))
p <- ggarrange(plotlist = plotlist,
               nrow = 3, ncol = 1, legend = "none")
p
ggsave("figures/Imd_expression_Tak1Tab2.png",
       width = 5, height = 6, plot = p)

# IKKB, kenny, Casper, and Relish:
ikkb <- getFlyBaseID("IKKbeta")
kenny <- getFlyBaseID("ken")
caspar <- getFlyBaseID("casp")
relish <- getFlyBaseID("Rel")
genedf |> filter(gene_name %in% c(ikkb, kenny, caspar, relish)) |> 
  select(gene_name, signal_to_noise_immune) |> 
  left_join(drosophila_lookup, by = c("gene_name"="FlyBaseID"))
plotlist <- plotExpressionProfileAllEnvironments(.gene_idxs = c(ikkb, kenny, caspar, relish), 
                                                 .gene_names = c("IKKB", "kenny", "Caspar", "Relish"),
                                                 .environments = c("oog", "dev", "immune"),
                                                 .color_vec = c(none_color, immune_color, immune_color, plei_color))
p <- ggarrange(plotlist = plotlist,
               nrow = 2, ncol = 2, legend = "none")
p
ggsave("figures/Imd_expression_NFKB.png",
       width = 6, height = 6, plot = p)

# Hep, basket, kayak
hemipterous <- getFlyBaseID("hep")
basket <- getFlyBaseID("bsk")
kayak <- getFlyBaseID("kay")
p <- ggarrange(plotlist = plotExpressionProfileAllEnvironments(.gene_idxs = c(hemipterous, basket, kayak), 
                                                               .gene_names = c("hemipterous/JNKK", "basket/JNK", "kayak"),
                                                               .environments = c("oog", "dev", "immune"),
                                                               .color_vec = c(none_color, plei_color, plei_color)),
               nrow = 3, ncol = 1, legend = "none")
p
ggsave("figures/Imd_expression_JNK.png",
       width = 5, height = 6, plot = p)

### Expression of Toll components
# with the caveat that the Imd challenge does not directly
# activate Toll, but there is cross-talk
immune_color <- "#F04D23"
plei_color <- "#0083C3"
none_color <- "grey60"

### Extracellular
pgrp_sa <- getFlyBaseID("PGRP-SA")
gnbp1 <- getFlyBaseID("GNBP1")
gnbp3 <- getFlyBaseID("GNBP3")
modsp <- getFlyBaseID("modSP")
grass <- getFlyBaseID("grass")
sphinx1 <- getFlyBaseID("sphinx1")
sphinx2 <- getFlyBaseID("sphinx2")
spheroide <- getFlyBaseID("sphe")
spirit <- getFlyBaseID("spirit")
plotExpressionProfileAllEnvironments(.gene_idxs = c(pgrp_sa, gnbp1, gnbp3, modsp, grass,
                                                    sphinx1, sphinx2, spheroide, spirit), 
                                     .gene_names = c("pgrp_sa", "gnbp1", "gnbp3", "modsp", "grass",
                                                     "sphinx1", "sphinx2", "spheroide", "spirit"),
                                     .environments = "immune",
                                     .color_vec = c(immune_color, immune_color, immune_color, 
                                                    immune_color, immune_color, immune_color,
                                                    immune_color, immune_color, none_color))

### MMPs
mmp1 <- getFlyBaseID("Mmp1")
mmp2 <- getFlyBaseID("Mmp2")
genedf |> filter(gene_name %in% c(mmp1, mmp2)) |> 
  select(williams_category)
ggarrange(plotlist = plotExpressionProfileAllEnvironments(.gene_idxs = c(mmp1, mmp2), 
                                     .gene_names = c("MMP1", "MMP2"),
                                     .environments = c("dev", "immune"),
                                     .color_vec = c(plei_color, plei_color)),
          nrow = 2, ncol = 1, legend = "none")
plotExpressionProfileAllEnvironments(.gene_idxs = c(mmp1, mmp2), 
                                     .gene_names = c("MMP1", "MMP2"),
                                     .environments = "immune",
                                     .color_vec = c(plei_color, plei_color))

### Intracellular, Toll components
spatzle <- getFlyBaseID("spz")
toll7 <-  getFlyBaseID("Toll-7")
toll6 <-  getFlyBaseID("Toll-6")
toll9 <-  getFlyBaseID("Toll-9")
toll18w <- getFlyBaseID("18w")
plotExpressionProfileAllEnvironments(.gene_idxs = c(spatzle, toll7, toll6, toll9,
                                                    toll18w), 
                                     .gene_names = c("spatzle", "toll7", "toll6", "toll9", "18-wheeler"),
                                     .environments = "immune",
                                     .color_vec = c(plei_color, plei_color, plei_color, plei_color,
                                                    plei_color))

### Gatekeepers of Death
diap1 <- getFlyBaseID("Diap1")
diap2 <- getFlyBaseID("Diap2")
dronc <- getFlyBaseID("Dronc")
drice <- getFlyBaseID("Drice")
bruce <- getFlyBaseID("Bruce")
plotlist <- plotExpressionProfileAllEnvironments(.gene_idxs = c(diap1, diap2, dronc,
                                                                drice, bruce), 
                                                 .gene_names = c("diap1", "diap2", "dronc",
                                                                 "drice", "bruce"),
                                                 .color_vec = c(plei_color, plei_color, plei_color, immune_color,
                                                                immune_color),
                                                 .environments = c("oog", "dev", "immune"))
ggarrange(plotlist = plotlist, nrow = 3, ncol = 1)

# just the Diaps
plotlist <- plotExpressionProfileAllEnvironments(.gene_idxs = c(diap1, diap2), 
                                                 .gene_names = c("diap1", "diap2"),
                                                 .color_vec = c(plei_color, plei_color),
                                                 .environments = c("dev", "immune"))
ggarrange(plotlist = plotlist, nrow = 2, ncol = 1, legend = "none")


#### Determining immune/dev variable quadrants ####
### Non-conservative/all-dev
plotdf <- genedf |> filter(williams_category != "Other")

# picking variability thresholds
threshdf <- plotdf |> select(gene_name, williams_category, 
                 signal_to_noise_immune, signal_to_noise_dev) |> 
  pivot_longer(cols = c("signal_to_noise_immune", "signal_to_noise_dev"),
               names_to = "environment", values_to = "signal_to_noise", 
               names_prefix = "signal_to_noise_") |>
  group_by(environment, williams_category) |> 
  summarise(quant50 = quantile(signal_to_noise, probs = 0.5))
SNR_thresh_immune <- threshdf |> 
  filter(environment == "immune" & 
           williams_category == "Immune_Non_Pleiotropic") |> 
  select(quant50) |> pull()
SNR_thresh_dev <- threshdf |> 
  filter(environment == "dev" & 
           williams_category == "Developmental_Non_Pleiotropic") |> 
  select(quant50) |> pull()
# QC: should be same number of genes less than and equal to threshold
sum(genedf$signal_to_noise_immune[genedf$williams_category == "Immune_Non_Pleiotropic"] > SNR_thresh_immune)
sum(genedf$signal_to_noise_immune[genedf$williams_category == "Immune_Non_Pleiotropic"] < SNR_thresh_immune)
sum(genedf$signal_to_noise_dev[genedf$williams_category == "Developmental_Non_Pleiotropic"] > SNR_thresh_dev)
sum(genedf$signal_to_noise_dev[genedf$williams_category == "Developmental_Non_Pleiotropic"] < SNR_thresh_dev)

# We're doing faceted histograms instead of scatterplots
# so that the SNR=0 genes can be in the same plot, and the 50-50 cutoff
# can be more apparent
p_dev <- ggplot(plotdf, aes(x = signal_to_noise_dev)) +
  geom_histogram(aes(fill = williams_category)) +
  geom_text(data = tibble(williams_category = c("Developmental_Non_Pleiotropic",
                                                "Immune_Non_Pleiotropic",
                                                "Pleiotropic"),
                          SNRthreshold = paste("SNR threshold = ",
                                               round(SNR_thresh_dev, digits = 2))),
            aes(label = SNRthreshold, x = SNR_thresh_dev + 1, y = 600),
            hjust = 0) +
  geom_text(data = group_by(plotdf, williams_category) |> 
              summarise(nAbove = sum(signal_to_noise_dev > SNR_thresh_dev)),
            aes(label = paste("# above = ", nAbove), y = 450, x = 25)) +
  geom_text(data = group_by(plotdf, williams_category) |> 
              summarise(nBelow = sum(signal_to_noise_dev <= SNR_thresh_dev)),
            aes(label = paste("# below = ", nBelow), y = 350, x = 25)) +
  facet_wrap(~williams_category) +
  geom_vline(xintercept = SNR_thresh_dev) +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("Signal to noise ratio, Embryonic development") +
  ylab("number of genes")
p_immune <- ggplot(plotdf, aes(x = signal_to_noise_immune)) +
  geom_histogram(aes(fill = williams_category)) +
  geom_text(data = tibble(williams_category = c("Developmental_Non_Pleiotropic",
                                                "Immune_Non_Pleiotropic",
                                                "Pleiotropic"),
                          SNRthreshold = paste("SNR threshold = ",
                                               round(SNR_thresh_immune, digits = 2))),
            aes(label = SNRthreshold, x = SNR_thresh_immune + 1, y = 750),
            hjust = 0) +
  geom_text(data = group_by(plotdf, williams_category) |> 
              summarise(nAbove = sum(signal_to_noise_immune > SNR_thresh_immune)),
            aes(label = paste("# above = ", nAbove), y = 450, x = 10)) +
  geom_text(data = group_by(plotdf, williams_category) |> 
              summarise(nBelow = sum(signal_to_noise_immune <= SNR_thresh_immune)),
            aes(label = paste("# below = ", nBelow), y = 350, x = 10)) +
  facet_wrap(~williams_category) +
  geom_vline(xintercept = SNR_thresh_immune) +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("Signal to noise ratio, Imd challenge") +
  ylab("number of genes")
ggarrange(p_dev, p_immune, nrow = 2, ncol = 1)
ggsave(filename = "figures/SNR_histograms.png", width = 7, height = 5)

# Accompanying stacked barplot of counts per williams category
plotdf_var <- plotdf |> 
  mutate(above_dev = signal_to_noise_dev > SNR_thresh_dev,
         above_immune = signal_to_noise_immune > SNR_thresh_immune) |> 
  mutate(var_group = if_else(above_dev,
                             true = if_else(above_immune,
                                            true = "var_both",
                                            false = "var_dev"),
                             false = if_else(above_immune,
                                             true = "var_immune",
                                             false = "var_neither")))
ggplot(plotdf_var, aes(y = factor(williams_category,
                                  levels = c("Developmental_Non_Pleiotropic",
                                             "Pleiotropic",
                                             "Immune_Non_Pleiotropic")))) +
  geom_bar(aes(fill = var_group), position = "fill", stat = "count") +
  scale_fill_brewer(palette = "Spectral") +
  theme_classic() +
  ylab("") +
  xlab("proportion of genes")
ggsave(filename = "figures/SNR_stackedBar.png", width = 7, height = 1.5)

### Conservative/embryonic (not sure if we need this b/c all the points are still included in the non-conservative version)
plotdf <- genedf |> 
  left_join(select(williams, FlyBaseID, Plei_allImmune_embDev),
            by = c("gene_name"="FlyBaseID")) |> 
  filter(!is.na(Plei_allImmune_embDev))
plotdf$williams_category <- plotdf$Plei_allImmune_embDev

# picking variability thresholds
threshdf <- plotdf |> select(gene_name, williams_category, 
                             signal_to_noise_immune, signal_to_noise_dev) |> 
  pivot_longer(cols = c("signal_to_noise_immune", "signal_to_noise_dev"),
               names_to = "environment", values_to = "signal_to_noise", 
               names_prefix = "signal_to_noise_") |>
  group_by(environment, williams_category) |> 
  summarise(quant50 = quantile(signal_to_noise, probs = 0.5))
SNR_thresh_immune <- threshdf |> 
  filter(environment == "immune" & 
           williams_category == "Immune_Non_Pleiotropic") |> 
  select(quant50) |> pull()
SNR_thresh_dev <- threshdf |> 
  filter(environment == "dev" & 
           williams_category == "Developmental_Non_Pleiotropic") |> 
  select(quant50) |> pull()
# QC: should be same number of genes less than and greater than threshold
sum(plotdf$signal_to_noise_immune[plotdf$williams_category == "Immune_Non_Pleiotropic"] > SNR_thresh_immune)
sum(plotdf$signal_to_noise_immune[plotdf$williams_category == "Immune_Non_Pleiotropic"] <= SNR_thresh_immune)
sum(plotdf$signal_to_noise_dev[plotdf$williams_category == "Developmental_Non_Pleiotropic"] > SNR_thresh_dev)
sum(plotdf$signal_to_noise_dev[plotdf$williams_category == "Developmental_Non_Pleiotropic"] <= SNR_thresh_dev)

# Scatterplots look too messy. We're doing faceted histograms instead
# so that the SNR=0 genes can be in the same plot and the 50-50 cutoff
# can be more apparent
p_dev <- ggplot(plotdf, aes(x = signal_to_noise_dev)) +
  geom_histogram(aes(fill = williams_category)) +
  geom_text(data = tibble(williams_category = c("Developmental_Non_Pleiotropic",
                                                "Immune_Non_Pleiotropic",
                                                "Pleiotropic"),
                          SNRthreshold = paste("SNR threshold = ",
                                               round(SNR_thresh_dev, digits = 2))),
            aes(label = SNRthreshold, x = SNR_thresh_dev + 1, y = 200),
            hjust = 0) +
  geom_text(data = group_by(plotdf, williams_category) |> 
              summarise(nAbove = sum(signal_to_noise_dev > SNR_thresh_dev)),
            aes(label = paste("# above = ", nAbove), y = 150, x = 13)) +
  geom_text(data = group_by(plotdf, williams_category) |> 
              summarise(nBelow = sum(signal_to_noise_dev <= SNR_thresh_dev)),
            aes(label = paste("# below = ", nBelow), y = 125, x = 13)) +
  facet_wrap(~williams_category) +
  geom_vline(xintercept = SNR_thresh_dev) +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("Signal to noise ratio, Embryonic development") +
  ylab("number of genes")
p_immune <- ggplot(plotdf, aes(x = signal_to_noise_immune)) +
  geom_histogram(aes(fill = williams_category)) +
  geom_text(data = tibble(williams_category = c("Developmental_Non_Pleiotropic",
                                                "Immune_Non_Pleiotropic",
                                                "Pleiotropic"),
                          SNRthreshold = paste("SNR threshold = ",
                                               round(SNR_thresh_immune, digits = 2))),
            aes(label = SNRthreshold, x = SNR_thresh_immune + 1, y = 200),
            hjust = 0) +
  geom_text(data = group_by(plotdf, williams_category) |> 
              summarise(nAbove = sum(signal_to_noise_immune > SNR_thresh_immune)),
            aes(label = paste("# above = ", nAbove), y = 150, x = 6)) +
  geom_text(data = group_by(plotdf, williams_category) |> 
              summarise(nBelow = sum(signal_to_noise_immune <= SNR_thresh_immune)),
            aes(label = paste("# below = ", nBelow), y = 125, x = 6)) +
  facet_wrap(~williams_category) +
  geom_vline(xintercept = SNR_thresh_immune) +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("Signal to noise ratio, Imd challenge") +
  ylab("number of genes")
ggarrange(p_dev, p_immune, nrow = 2, ncol = 1)
ggsave(filename = "figures/SNR_histograms_embryonic.png", width = 7, height = 5)

# Accompanying stacked barplot of counts per williams category
plotdf_var_emb <- plotdf |> 
  mutate(above_dev = signal_to_noise_dev > SNR_thresh_dev,
         above_immune = signal_to_noise_immune > SNR_thresh_immune) |> 
  mutate(var_group = if_else(above_dev,
                             true = if_else(above_immune,
                                            true = "var_both",
                                            false = "var_dev"),
                             false = if_else(above_immune,
                                             true = "var_immune",
                                             false = "var_neither")))
ggplot(plotdf_var_emb, aes(y = factor(williams_category,
                                      levels = c("Developmental_Non_Pleiotropic",
                                                 "Pleiotropic",
                                                 "Immune_Non_Pleiotropic")))) +
  geom_bar(aes(fill = var_group), position = "fill", stat = "count") +
  scale_fill_brewer(palette = "Spectral") +
  theme_classic() +
  ylab("") +
  xlab("proportion of genes")
ggsave(filename = "figures/SNR_stackedBar_embryonic.png", width = 7, height = 1.5)

# embryonic and all-development plots have the most dramatically different proportions
# for the pleiotropic genes:
plotdf_var$isEmbryonic <- FALSE
plotdf_var_emb$isEmbryonic <- TRUE
plotdf_var_full <- bind_rows(plotdf_var, plotdf_var_emb)
ggplot(plotdf_var_full, aes(y = isEmbryonic)) +
  geom_bar(aes(fill = var_group), position = "fill", stat = "count") +
  scale_fill_brewer(palette = "Spectral") +
  theme_classic() +
  facet_wrap(~williams_category) +
  ylab("is Embryonic") +
  xlab("proportion of genes")
