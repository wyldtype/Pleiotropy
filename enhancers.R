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
drosophila_lookup <- read_delim("data/immune_FlyBaseLookup.txt", delim = "\t")
colnames(drosophila_lookup) <- c("systematic_name", "full_common_name", "common_name")

# loading snr values and coefficients
load("data/SignalToNoise.RData")
SNR_thresh # decided in signal_to_noise.R

# Cohen et al. 2026, Imd enhancers
# loading enhancer assignments
enhancerdf_Cohen <- read_csv("data/EG_TFBS_cat_01_19_activityclass.csv") |> 
  filter(Gene != ".") # . are enhancers not assigned to any gene
table(enhancerdf_Cohen$Group)
#enhancerdf_allAccess <- read_csv("data/EG_all_access_07_28.csv")

# Kvon et al. 2014, Developmental enhancers
enhancerdf_Kvon <- read_csv("data/Kvon2014/2014-01-00083C-Supplementary Table 4.csv")
names(enhancerdf_Kvon) <- c("VTID", "Chrom_enhancer", "Start_enhancer", "End_enhancer",
                            "FlyBaseID", "Symbol", "Chrom_gene", "Start_gene", "End_gene",
                            "Orientation", "Match")

# creating dataframe of enhancers per gene (all dev enhancers)
nEnhancerdf_dev <- enhancerdf_Kvon |> 
  select(VTID, FlyBaseID) |> 
  group_by(FlyBaseID) |> 
  summarise(nEnhancers = n())

# creating dataframe of enhancers per gene (HKSM_only-enhancers)
nEnhancerdf_hksm <- enhancerdf_Cohen |> 
  filter(Group %in% c("HKSM_only") &
           Accessibility != "Always closed") |> 
  select(Enhancer, Gene) |> 
  separate_longer_delim(col = "Gene", delim = ";") |> 
  group_by(Gene) |> 
  summarise(nEnhancers = n())
# HKSM or 20E+HKSM enhancers
nEnhancerdf_20E_hksm <- enhancerdf_Cohen |> 
  filter(Group %in% c("HKSM_only", "HKSMn20E") &
           Accessibility != "Always closed") |> 
  select(Enhancer, Gene) |> 
  separate_longer_delim(col = "Gene", delim = ";") |> 
  group_by(Gene) |> 
  summarise(nEnhancers = n())

#### Positive control 1: Do developmental genes have more enhancers? ####

# all developmental enhancers
plotdf <- left_join(genedf, nEnhancerdf_dev,
                    by = c("gene_name"="FlyBaseID"))  |>  
  filter(williams_category != "Other") |> 
  mutate(nEnhancers = coalesce(nEnhancers, 0)) # genes missing from nEnhancerdf have 0 known enhancers

# % of genes in each category with at least one enhancer:
plotdf_pcts <- plotdf |> 
  group_by(williams_category) |> 
  summarise(pct_with_enhancer = (sum(nEnhancers > 0)/n())*100)

# gene type defined in Williams et al.
ggplot(plotdf,
       aes(x = factor(williams_category,
                      levels = c("Developmental_Non_Pleiotropic",
                                 "Pleiotropic",
                                 "Immune_Non_Pleiotropic")), 
           y = nEnhancers)) +
  geom_jitter(aes(color = williams_category), height = 0) +
  geom_text(data = plotdf_pcts,
            aes(x = williams_category, y = 16, 
                label = paste0("% of genes with at least\none enhancer: ",
                               round(pct_with_enhancer, digits = 2),
                               "%"))) +
  geom_violin() +
  theme_classic() +
  theme(legend.position = "none") +
  ylab("number of enhancers per gene") +
  xlab("gene category") +
  ggtitle("Developmental enhancers")

#### Pos-control 2: Do immune genes have more HKSM-only enhancers? ####

# HKSM-only enhancers
plotdf <- left_join(genedf, nEnhancerdf_hksm,
                    by = c("gene_name"="Gene")) |>  
  filter(williams_category != "Other") |> 
  mutate(nEnhancers = coalesce(nEnhancers, 0)) # genes missing from nEnhancerdf have 0 known enhancers


plotdf_pcts <- plotdf |> 
  group_by(williams_category) |> 
  summarise(pct_with_enhancer = (sum(nEnhancers > 0)/n())*100)

# gene type defined in Williams et al.
ggplot(plotdf,
       aes(x = factor(williams_category,
                      levels = c("Developmental_Non_Pleiotropic",
                                 "Pleiotropic",
                                 "Immune_Non_Pleiotropic")), 
           y = nEnhancers)) +
  geom_jitter(aes(color = williams_category), height = 0) +
  geom_text(data = plotdf_pcts,
            aes(x = williams_category, y = 5, 
                label = paste0("% of genes with at least\none enhancer: ",
                               round(pct_with_enhancer, digits = 2),
                               "%"))) +
  geom_violin() +
  theme_classic() +
  theme(legend.position = "none") +
  ylab("number of enhancers per gene") +
  xlab("gene category") +
  ggtitle("HKSM-only enhancers")

# Bar plot of just percents
plotdf_pcts_dev <- left_join(genedf, nEnhancerdf_dev,
                              by = c("gene_name"="FlyBaseID")) |>  
  filter(williams_category != "Other") |> 
  mutate(nEnhancers = coalesce(nEnhancers, 0)) |> 
  group_by(williams_category) |> 
  summarise(pct_with_enhancer = (sum(nEnhancers > 0)/n())*100) |> 
  mutate(enhancer_type = "Developmental Enhancers")
plotdf_pcts_hksm <- left_join(genedf, nEnhancerdf_hksm, # (nEnhancerdf_20E_hksm produces similar, there's just a lot more of them)
                                  by = c("gene_name"="Gene")) |>  
  filter(williams_category != "Other") |> 
  mutate(nEnhancers = coalesce(nEnhancers, 0)) |> 
  group_by(williams_category) |> 
  summarise(pct_with_enhancer = (sum(nEnhancers > 0)/n())*100) |> 
  mutate(enhancer_type = "HKSM-only Enhancers")

plotdf_pcts <- bind_rows(plotdf_pcts_dev, plotdf_pcts_hksm)

p <- ggplot(plotdf_pcts, aes(x = williams_category, y = pct_with_enhancer/100)) +
  geom_bar(aes(fill = williams_category), stat = "identity") +
  scale_fill_brewer(palette = "Set1") +
  theme_classic() +
  facet_wrap(~enhancer_type) +
  ylab("fraction of genes with at least one enhancer") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p
# figure pdf:
pdf("figures/enhancer_pcts_bar.pdf",
    width = 7, height = 5)
p
dev.off()

#### How many pleiotropic genes have developmental and immune enhancers? ####
# About 7% of pleiotropic genes have connections to developmental or immune enhancers.
# But are these different sets of 7?
plotdf <- genedf |> filter(williams_category != "Other") |> 
  left_join(nEnhancerdf_dev,
            by = c("gene_name"="FlyBaseID")) |> 
  dplyr::rename(c("nEnhancers_dev" = "nEnhancers")) |> 
  mutate(nEnhancers_dev = coalesce(nEnhancers_dev, 0)) |> 
  left_join(nEnhancerdf_hksm,
            by = c("gene_name"="Gene")) |> 
  dplyr::rename(c("nEnhancers_hksm" = "nEnhancers")) |> 
  mutate(nEnhancers_hksm = coalesce(nEnhancers_hksm, 0)) 

library(ggVennDiagram)
plotlist <- list(dev_enhancers = plotdf |> 
                   filter(williams_category == "Pleiotropic") |> 
                   filter(nEnhancers_dev > 0) |> 
                   select(gene_name) |> pull(),
                 immune_enhancers = plotdf |> 
                   filter(williams_category == "Pleiotropic") |> 
                   filter(nEnhancers_hksm > 0) |> 
                   select(gene_name) |> pull())
p <- ggVennDiagram(plotlist) + scale_fill_gradient(low="grey90",high = "red")
p # one gene
plotdf |> filter(williams_category == "Pleiotropic" &
                   nEnhancers_dev > 0 &
                   nEnhancers_hksm > 0) # FBgn0260635 = Diap1
plotExpressionProfile(.counts = counts$immune, .info = infodf$immune,
                      .gene_idxs = "FBgn0260635", .gene_names = "Diap1")
plotExpressionProfile(.counts = counts$dev, .info = infodf$dev,
                      .gene_idxs = "FBgn0260635", .gene_names = "Diap1")

#### Do more variable genes have more enhancers? ####
# for each type of gene: Developmental, Immune, or Pleiotropic, separate each by whether
# it has at least one enhancer in the environment where it is known to play a role
plotdf <- genedf |> filter(williams_category != "Other") |> 
  left_join(nEnhancerdf_dev,
            by = c("gene_name"="FlyBaseID")) |> 
  dplyr::rename(c("nEnhancers_dev" = "nEnhancers")) |> 
  mutate(nEnhancers_dev = coalesce(nEnhancers_dev, 0)) |> 
  left_join(nEnhancerdf_hksm,
            by = c("gene_name"="Gene")) |> 
  dplyr::rename(c("nEnhancers_hksm" = "nEnhancers")) |> 
  mutate(nEnhancers_hksm = coalesce(nEnhancers_hksm, 0)) 

# plotdf$onePlusEnhancer <- map(c(1:nrow(plotdf)), \(i) {
#   category <- plotdf$williams_category[i]
#   onePlusEnhancer <- NA
#   if (category == "Developmental_Non_Pleiotropic") {
#     onePlusEnhancer <- plotdf$nEnhancers_dev[i] > 0
#   }
#   if (category == "Immune_Non_Pleiotropic") {
#     onePlusEnhancer <- plotdf$nEnhancers_hksm[i] > 0
#   }
#   if (category == "Pleiotropic") {
#     onePlusEnhancer <- plotdf$nEnhancers_dev[i] > 0 |
#       plotdf$nEnhancers_hksm[i] > 0
#   }
#   return(onePlusEnhancer)
# }) |> unlist()

# nEnhancers in dev, SNR in Embryonic Development
plotdf_counts <- plotdf |> mutate(hasEnhancers = nEnhancers_dev > 0) |> 
  group_by(hasEnhancers, williams_category) |> 
  summarise(n = n())
ggplot(data = mutate(plotdf, hasEnhancers = nEnhancers_dev > 0),
       aes(x = hasEnhancers, y = signal_to_noise_dev)) +
  geom_boxplot(aes(group = hasEnhancers,
                   fill = williams_category)) +
  geom_text(data = plotdf_counts,
            aes(x = hasEnhancers, y = 30, label = n)) +
  facet_wrap(~williams_category) +
  theme_classic() +
  xlab("Has at least one enhancer in Embryonic Development")

# nEnhancers in immune, SNR in immune
plotdf_counts <- plotdf |> mutate(hasEnhancers = nEnhancers_hksm > 0) |> 
  group_by(hasEnhancers, williams_category) |> 
  summarise(n = n())
ggplot(data = mutate(plotdf, hasEnhancers = nEnhancers_hksm > 0),
       aes(x = hasEnhancers, y = signal_to_noise_immune)) +
  geom_boxplot(aes(group = hasEnhancers,
                   fill = williams_category)) +
  geom_text(data = plotdf_counts,
            aes(x = hasEnhancers, y = 10, label = n)) +
  facet_wrap(~williams_category) +
  theme_classic() +
  xlab("Has at least one enhancer in HKSM")

# what are those immune genes that seem so dev active?
plotgenes <- plotdf |> filter(nEnhancers_dev > 0 & 
                   williams_category == "Immune_Non_Pleiotropic") |> 
  left_join(drosophila_lookup, by = c("gene_name"="systematic_name")) |> 
  select(gene_name, common_name, signal_to_noise_dev, signal_to_noise_immune)
plotExpressionProfile(.counts = counts$dev, .info = infodf$dev,
                      .gene_idxs = plotgenes$gene_name, 
                      .gene_names = plotgenes$common_name)
plotExpressionProfile(.counts = counts$immune, .info = infodf$immune,
                      .gene_idxs = plotgenes$gene_name, 
                      .gene_names = plotgenes$common_name)

# tables of nEnhancers 0 or 1+ and SNR <= 2 or > 2
# dev
plotdf |> select(nEnhancers_dev, signal_to_noise_dev) |> 
  mutate(onePlusEnhancer = nEnhancers_dev > 0,
         SNRup = signal_to_noise_dev > SNR_thresh) |> 
  group_by(onePlusEnhancer, SNRup) |> 
  summarise(n = n())
(123/50)/(1453/948) # if you have an enhancer, you are 1.6 times more likely to have high SNR
fisher.test(cbind(c(123, 50), c(1453, 948)))

# immune
plotdf |> select(nEnhancers_hksm, signal_to_noise_immune) |> 
  mutate(onePlusEnhancer = nEnhancers_hksm > 0,
         SNRup = signal_to_noise_immune > 1) |> 
  group_by(onePlusEnhancer, SNRup) |> 
  summarise(n = n())
(101/30)/(1945/498)
fisher.test(cbind(c(101, 30), c(1945, 498))) # no direction

# The genes with more enhancers are more variable in dev, this makes sense
# But not so in immune. Why not? There are so few williams-classified immune genes 
# with enhancers, we can probably just look at each one manually

# Looking manually below suggested that immune genes with enhancers might have
# higher mean expression, rather than higher snr
# nEnhancers in immune, mean in immune
plotdf_counts <- plotdf |> mutate(hasEnhancers = nEnhancers_hksm > 0) |> 
  group_by(hasEnhancers, williams_category) |> 
  summarise(n = n())
ggplot(data = mutate(plotdf, hasEnhancers = nEnhancers_hksm > 0),
       aes(x = hasEnhancers, y = log2(mean_immune))) +
  geom_boxplot(aes(group = hasEnhancers,
                   fill = williams_category)) +
  geom_text(data = plotdf_counts,
            aes(x = hasEnhancers, y = 10, label = n)) +
  facet_wrap(~williams_category) +
  theme_classic() +
  xlab("Has at least one enhancer in HKSM") # yes for dev and plei, maybe for immune

# is that true in development too?
plotdf_counts <- plotdf |> mutate(hasEnhancers = nEnhancers_dev > 0) |> 
  group_by(hasEnhancers, williams_category) |> 
  summarise(n = n())
ggplot(data = mutate(plotdf, hasEnhancers = nEnhancers_dev > 0),
       aes(x = hasEnhancers, y = log2(mean_dev))) +
  geom_boxplot(aes(group = hasEnhancers,
                   fill = williams_category)) +
  geom_text(data = plotdf_counts,
            aes(x = hasEnhancers, y = 14, label = n)) +
  facet_wrap(~williams_category) +
  theme_classic() +
  xlab("Has at least one enhancer in Embryonic Development") # not for dev genes!

#### Inspecting expression of immune genes with enhancers ####
### Expression of Imd components from Williams et al. figure in 
# Imd challenge, Embyonic Dev, and Oogenesis
immune_color <- "#F04D23"
plei_color <- "#0083C3"

# Extracellular
p <- ggarrange(plotlist = plotExpressionProfileAllEnvironments(.gene_idxs = c("FBgn0035806", "FBgn0030695"), 
                                                          .gene_names = c("PGRP-SD", "PGRP-LE"),
                                                          .environments = c("oog", "dev", "immune"),
                                                          .color_vec = c(immune_color, plei_color)),
          nrow = 3, ncol = 1, legend = "none")
p
ggsave("figures/Imd_expression_extracellular.png",
    width = 3, height = 6, plot = p)

# Membrane bound
p <- ggarrange(plotlist = plotExpressionProfileAllEnvironments(.gene_idxs = c("FBgn0035976", "FBgn0035975", "FBgn0035977"), 
                                                          .gene_names = c("PGRP-LC", "PGRP-LA", "PGRP-LF"),
                                                          .environments = c("oog", "dev", "immune"),
                                                          .color_vec = c(immune_color, immune_color, plei_color)),
          nrow = 3, ncol = 1, legend = "none")
p
ggsave("figures/Imd_expression_membrane.png",
       width = 5, height = 6, plot = p)
# Intracellular pinwheel complex:
ggarrange(plotlist = plotExpressionProfileAllEnvironments(.gene_idxs = c("FBgn0035601", "FBgn0000173", "FBgn0011217", "FBgn0015247"), 
                                                          .gene_names = c("Uev1A", "bendless", "effete", "Diap2"),
                                                          .environments = c("oog", "dev", "immune"),
                                                          .color_vec = c(immune_color, plei_color, plei_color, plei_color)),
          nrow = 1, ncol = 3, legend = "none")

# Diap2 vs Diap1 (cause Diap1 has the double enhancer situation FBgn0260635 = Diap1)
ggarrange(plotlist = plotExpressionProfileAllEnvironments(.gene_idxs = c("FBgn0260635", "FBgn0015247"), 
                                                          .gene_names = c("Diap1", "Diap2"),
                                                          .environments = c("oog", "dev", "immune"),
                                                          .color_vec = c(plei_color, plei_color)),
          nrow = 1, ncol = 3, legend = "none")

# Intracellular pentagram complex (all plei so not that helpful):
ggarrange(plotlist = plotExpressionProfileAllEnvironments(.gene_idxs = c("FBgn0013983", "FBgn0020381", "FBgn0038928"), 
                                                          .gene_names = c("Imd", "Dredd", "Fadd"),
                                                          .environments = c("oog", "dev", "immune"),
                                                          .color_vec = c(plei_color, plei_color, plei_color, plei_color,
                                                                         plei_color, plei_color, plei_color, plei_color)),
          nrow = 1, ncol = 3, legend = "none")

# Tak1 and Tab2
ggarrange(plotlist = plotExpressionProfileAllEnvironments(.gene_idxs = c("FBgn0086358","FBgn0026323"), 
                                                          .gene_names = c("Tab2", "Tak1"),
                                                          .environments = c("oog", "dev", "immune"),
                                                          .color_vec = c(immune_color, plei_color)),
          nrow = 1, ncol = 3, legend = "none")

# Casper, kenny, basket, kayak, and Relish:
ggarrange(plotlist = plotExpressionProfileAllEnvironments(.gene_idxs = c("FBgn0034068", "FBgn0041205", "FBgn0014018", "FBgn0001297", "FBgn0000229"), 
                                                          .gene_names = c("Caspar", "kenny", "Relish", "kayak", "basket/JNK"),
                                                          .environments = c("oog", "dev", "immune"),
                                                          .color_vec = c(immune_color, immune_color, plei_color, plei_color, plei_color)),
          nrow = 1, ncol = 3, legend = "none")




#### visualizing genes by SNR ####
# immune genes with enhancers, SNR > 1:
gene_idxs <- plotdf |> filter(nEnhancers_hksm > 0 & 
                                signal_to_noise_immune > 1 &
                                williams_category == "Immune_Non_Pleiotropic") |> 
  select(gene_name) |> pull()
gene_names <- map(gene_idxs, \(g) {drosophila_lookup$common_name[which(drosophila_lookup$systematic_name == g)]}) |> 
  unlist()
plotExpressionProfile(.counts = counts$immune, .info = infodf$immune,
                      .gene_idxs = gene_idxs, .gene_names = gene_names)
# SNR <= 1:
gene_idxs <- plotdf |> filter(nEnhancers_hksm > 0 & 
                                signal_to_noise_immune <= 1 &
                                williams_category == "Immune_Non_Pleiotropic") |> 
  select(gene_name) |> pull()
gene_names <- map(gene_idxs, \(g) {drosophila_lookup$common_name[which(drosophila_lookup$systematic_name == g)]}) |> 
  unlist()
plotExpressionProfile(.counts = counts$immune, .info = infodf$immune,
                      .gene_idxs = gene_idxs, .gene_names = gene_names)
# that checks out

# conclusion: these don't have particularly high SNR or any uniform behavior. 
# The PGRPs are noisy between replicates. The most interesting genes are 3/4 of the
# un-named ones: CG16985 and CG16986 both increase steadily in both replicates and 
# CG17271 has a robust spike fairly early

#### Inspecting expression of immune genes with enhancers ####
# Pleiotropic genes with dev enhancers, SNR > 1
gene_idxs <- plotdf |> filter(nEnhancers_dev > 0 & 
                                signal_to_noise_dev > 1 &
                                williams_category == "Pleiotropic") |> 
  select(gene_name) |> pull()
plotExpressionProfile(.counts = counts$dev, .info = infodf$dev,
                      .gene_idxs = gene_idxs, .gene_names = gene_idxs)

# Pleiotropic genes with dev enhancers, SNR <= 1
gene_idxs <- plotdf |> filter(nEnhancers_dev > 0 & 
                                signal_to_noise_dev <= 1 &
                                williams_category == "Pleiotropic") |> 
  select(gene_name) |> pull()
plotExpressionProfile(.counts = counts$dev, .info = infodf$dev,
                      .gene_idxs = gene_idxs, .gene_names = gene_idxs)

# Pleiotropic genes with immune enhancers, SNR > 1
gene_idxs <- plotdf |> filter(nEnhancers_hksm > 0 & 
                                signal_to_noise_immune > 1 &
                                williams_category == "Pleiotropic") |> 
  select(gene_name) |> pull()
plotExpressionProfile(.counts = counts$immune, .info = infodf$immune,
                      .gene_idxs = gene_idxs, .gene_names = gene_idxs) # FBgn0039562 = Glycoprotein 93

# Pleiotropic genes with immune enhancers, SNR <= 1
plotExpressionProfile(.counts = counts$immune, .info = infodf$immune,
                      .gene_idxs = "FBgn0087035", .gene_names = "FBgn0087035") # highly expressed
gene_idxs <- plotdf |> filter(nEnhancers_hksm > 0 & 
                                signal_to_noise_immune <= 1 &
                                williams_category == "Pleiotropic") |> 
  select(gene_name) |> pull() |> setdiff("FBgn0087035")
plotExpressionProfile(.counts = counts$immune, .info = infodf$immune,
                      .gene_idxs = gene_idxs, .gene_names = gene_idxs)
# two of these do look fairly variable FBgn0004396 (CrebA) and FBgn0011274 (Dif):
plotdf |> filter(nEnhancers_hksm > 0 & 
                   signal_to_noise_immune <= 1 &
                   williams_category == "Pleiotropic") |> 
  arrange(desc(signal_to_noise_immune))

