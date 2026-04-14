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

#### Are immune + pleiotropic genes more variable than dev genes in Imd challenge? ####
# In the adult, there aren't as many widespread changes in expression,
# so only the genes that need to change are changing
plotdf <- genedf |> select(signal_to_noise_immune, signal_to_noise_dev, williams_category) |> 
  filter(williams_category != "Other") |> unique()
p_imd <- ggplot(plotdf, aes(x = williams_category, y = signal_to_noise_immune)) +
  geom_boxplot(aes(fill = williams_category)) +
  stat_compare_means(comparisons = list(c("Developmental_Non_Pleiotropic", "Pleiotropic"),
                                        c("Immune_Non_Pleiotropic", "Pleiotropic"),
                                        c("Developmental_Non_Pleiotropic", "Immune_Non_Pleiotropic"))) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Signal to noise - Imd challenge")
p_imd
ggsave(p_imd, filename = "figures/immune_variability_box.png",
       width = 3, height = 7)

# this is in contrast to development, where most genes are quite variable
p_dev <- ggplot(plotdf, aes(x = williams_category, y = signal_to_noise_dev)) +
  geom_boxplot(aes(fill = williams_category)) +
  stat_compare_means(comparisons = list(c("Developmental_Non_Pleiotropic", "Pleiotropic"),
                                        c("Immune_Non_Pleiotropic", "Pleiotropic"),
                                        c("Developmental_Non_Pleiotropic", "Immune_Non_Pleiotropic"))) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Signal to noise - Embryonic Development")
p_dev
ggsave(p_dev, filename = "figures/dev_variability_box.png",
       width = 3, height = 7)

#### Example genes from Williams et al. ####

### Expression of Imd components from Williams et al. figure in 
# Imd challenge, Embyonic Dev, and Oogenesis
immune_color <- "#F04D23"
plei_color <- "#0083C3"
none_color <- "grey60"

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

#### Can we assign pleiotropic genes to Development or Immune based on enhancers/expression? ####
# Our enhancer exploration has 

