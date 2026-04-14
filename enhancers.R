sapply(c("tidyr", "dplyr", "ggplot2", "ggpubr", "ggalluvial", "purrr", "readxl", "readr", "stringr", "plotly", "lmtest", "pheatmap"), FUN = require,
       character.only = TRUE)

# loading counts and functions
source("utils.r")
load("data/CountsList.RData")

# williams dataset of gene functions
williams <- read_xlsx(path = "data/17Jan2023_Supplemental_File_1.xlsx",
                      sheet = 2, na = "NA") |> 
  select(FlyBaseID, Plei_allImmune_allDev) 

# FlyBase gene info
drosophila_lookup <- read_delim("data/best_gene_summary_fb_2026_01.tsv", delim = "\t",
                                comment = "##")
colnames(drosophila_lookup) <- c("FlyBaseID", "Gene_Symbol", "Summary_Source", "Summary")

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

# Cohen et al. 2026, Imd enhancers
# loading enhancer assignments
enhancerdf_Cohen <- read_csv("data/EG_TFBS_cat_01_19_activityclass.csv") |> 
  filter(Gene != ".")  # . are enhancers not assigned to any gene
# If we need to parse the Enhancer Name:
#   mutate(Enhancer_Start = str_extract(Enhancer, "(?<=:)\\d+(?=-)"),
#          Enhancer_End = str_extract(Enhancer, "(?<=-)\\d+$"),
#          Chromosome = str_extract(Enhancer, "^.*(?=:)"))
table(enhancerdf_Cohen$Group)

#enhancerdf_allAccess <- read_csv("data/EG_all_access_07_28.csv")

# Kvon et al. 2014, Developmental enhancers
enhancerdf_Kvon <- read_csv("data/Kvon2014/2014-01-00083C-Supplementary Table 4.csv")
names(enhancerdf_Kvon) <- c("VTID", "Chrom_enhancer", "Start_enhancer", "End_enhancer",
                            "FlyBaseID", "Symbol", "Chrom_gene", "Start_gene", "End_gene",
                            "Orientation", "Match")

# creating dataframe of enhancers per gene (all conditions)

# creating dataframe of enhancers per gene (all dev enhancers)
nEnhancerdf_dev <- enhancerdf_Kvon |> 
  select(VTID, FlyBaseID) |> 
  group_by(FlyBaseID) |> 
  summarise(nEnhancers = n()) |> 
  mutate(condition = "dev")

# creating dataframe of enhancers per gene (HKSM_only-enhancers)
nEnhancerdf_hksm <- enhancerdf_Cohen |> 
  filter(Group %in% c("HKSM_only") &
           Accessibility != "Always closed") |> 
  select(Enhancer, Gene) |>
  dplyr::rename(c("FlyBaseID"="Gene")) |> 
  separate_longer_delim(col = "FlyBaseID", delim = ";") |> 
  group_by(FlyBaseID) |> 
  summarise(nEnhancers = n()) |> 
  mutate(condition = "HKSM_only")

# Enhancers that come on when 20E is added, regardless of whether HKSM is added
nEnhancerdf_20E <- enhancerdf_Cohen |> 
  filter(Group %in% c("HKSMn20E", "20E_only") &
           Accessibility != "Always closed") |> 
  select(Enhancer, Gene) |> 
  dplyr::rename(c("FlyBaseID"="Gene")) |> 
  separate_longer_delim(col = "FlyBaseID", delim = ";") |> 
  group_by(FlyBaseID) |> 
  summarise(nEnhancers = n()) |> 
  mutate(condition = "20E")

# Constitutive enhancers (control for Plei having more hksm enhancers)
nEnhancerdf_constit <- enhancerdf_Cohen |> 
  filter(Group %in% c("Constit") &
           Accessibility != "Always closed") |> 
  select(Enhancer, Gene) |> 
  dplyr::rename(c("FlyBaseID"="Gene")) |> 
  separate_longer_delim(col = "FlyBaseID", delim = ";") |> 
  group_by(FlyBaseID) |> 
  summarise(nEnhancers = n()) |>
  mutate(condition = "s2_constit")

# Control enhancers that aren't constitutive (i.e. off in at least one of HKSM or 20E)
nEnhancerdf_control <- enhancerdf_Cohen |> 
  filter(Group %in% c("Con_only", "Controln20E", "HKSMnCon") &
           Accessibility != "Always closed") |> 
  select(Enhancer, Gene) |> 
  dplyr::rename(c("FlyBaseID"="Gene")) |> 
  separate_longer_delim(col = "FlyBaseID", delim = ";") |> 
  group_by(FlyBaseID) |> 
  summarise(nEnhancers = n()) |>
  mutate(condition = "s2_control")

# nEnhancerdf for all conditions
nEnhancerdf <- bind_rows(nEnhancerdf_dev, nEnhancerdf_hksm,
                         nEnhancerdf_20E, nEnhancerdf_constit,
                         nEnhancerdf_control)

#### Saving enhancer datasets ####
save(nEnhancerdf, enhancerdf_Cohen, enhancerdf_Kvon,
     file = "data/Enhancers.RData")

#### Visualizing enhancer stats ####
# TODO: this plot can probably be cleaned up with the combined nEnhancerdf
# number of enhancers in each dataset
p_en <- ggplot(tibble(condition = c("develompent", "constit",
                                    "HKSM", "20E and HKSM"),
                      nEnhancers = c(length(unique(enhancerdf_Kvon$VTID)),
                                     length(unique(pull(select(filter(enhancerdf_Cohen,
                                                                      Group %in% c("Constit") &
                                                                        Accessibility != "Always closed"),
                                                               Enhancer)))),
                                     length(unique(pull(select(filter(enhancerdf_Cohen,
                                                                      Group %in% c("HKSM_only") &
                                                                        Accessibility != "Always closed"),
                                                               Enhancer)))),
                                     length(unique(pull(select(filter(enhancerdf_Cohen,
                                                                      Group %in% c("HKSM_only", "HKSMn20E") &
                                                                        Accessibility != "Always closed"),
                                                               Enhancer)))))),
               aes(x = factor(condition, levels = c("develompent", "constit", "HKSM", "20E and HKSM"),
                              labels = c("Develompent", "Constitutive (S2*)", "HKSM", "20E and HKSM")), 
                   y = nEnhancers)) +
  geom_bar(stat = "identity", color = "black", fill = "grey80") +
  geom_text(aes(label = nEnhancers, y = nEnhancers + 100)) +
  ylim(c(0, 1600)) +
  xlab("") +
  ylab("# unique enhancers") +
  theme_classic()
p_en

# number of genes targeted by enhancers in each dataset
p_g <- ggplot(tibble(condition = c("develompent", "constit", "HKSM", "20E and HKSM"),
              nEnhancers = c(nrow(nEnhancerdf_dev),
                             nrow(nEnhancerdf_constit),
                             nrow(nEnhancerdf_hksm),
                             nrow(nEnhancerdf_20E_hksm))),
       aes(x = factor(condition, levels = c("develompent", "constit", "HKSM", "20E and HKSM"),
                      labels = c("Develompent", "Constitutive (S2*)", "HKSM", "20E and HKSM")), 
           y = nEnhancers)) +
  geom_bar(stat = "identity", color = "black", fill = "grey80") +
  geom_text(aes(label = nEnhancers, y = nEnhancers + 100)) +
  ylim(c(0, 2000)) +
  xlab("") +
  ylab("# genes targeted by enhancers") +
  theme_classic()
p_g
p <- ggarrange(p_en, p_g, nrow = 1, ncol = 2)
ggsave(plot = p, filename = "figures/nEnhancers.png",
       width = 6, height = 3)


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
p <- ggplot(plotdf,
       aes(x = factor(williams_category,
                      levels = c("Developmental_Non_Pleiotropic",
                                 "Pleiotropic",
                                 "Immune_Non_Pleiotropic")), 
           y = nEnhancers)) +
  geom_jitter(aes(color = williams_category), height = 0) +
  geom_text(data = plotdf_pcts,
            aes(x = williams_category, y = 15, 
                label = paste0("% of genes with at least\none enhancer: ",
                               round(pct_with_enhancer, digits = 2),
                               "%"))) +
  geom_violin() +
  theme_classic() +
  theme(legend.position = "none") +
  ylab("number of enhancers per gene") +
  xlab("gene category") +
  ggtitle("Developmental enhancers")

p
ggsave("figures/nEnhancers_dev_scatter.png", plot = p,
       width = 6, height = 4)

#### Violin plots of number of enhancers per gene ####

# HKSM-only enhancers
plotdf <- left_join(genedf, nEnhancerdf_hksm,
                    by = c("gene_name"="Gene")) |>  
  filter(williams_category != "Other") |> 
  mutate(nEnhancers = coalesce(nEnhancers, 0)) # genes missing from nEnhancerdf have 0 known enhancers


plotdf_pcts <- plotdf |> 
  group_by(williams_category) |> 
  summarise(pct_with_enhancer = (sum(nEnhancers > 0)/n())*100)

# gene type defined in Williams et al.
p <- ggplot(plotdf,
       aes(x = factor(williams_category,
                      levels = c("Developmental_Non_Pleiotropic",
                                 "Pleiotropic",
                                 "Immune_Non_Pleiotropic")), 
           y = nEnhancers)) +
  geom_jitter(aes(color = williams_category), height = 0) +
  geom_text(data = plotdf_pcts,
            aes(x = williams_category, y = 3, 
                label = paste0("% of genes with at least\none enhancer: ",
                               round(pct_with_enhancer, digits = 2),
                               "%"))) +
  geom_violin() +
  theme_classic() +
  theme(legend.position = "none") +
  ylab("number of enhancers per gene") +
  xlab("gene category") +
  ggtitle("HKSM-only enhancers")

p
ggsave("figures/nEnhancers_hksm_scatter.png", plot = p,
       width = 6, height = 4)

# Constitutive S2* enhancers
plotdf <- left_join(genedf, nEnhancerdf_constit,
                    by = c("gene_name"="Gene")) |>  
  filter(williams_category != "Other") |> 
  mutate(nEnhancers = coalesce(nEnhancers, 0)) # genes missing from nEnhancerdf have 0 known enhancers


plotdf_pcts <- plotdf |> 
  group_by(williams_category) |> 
  summarise(pct_with_enhancer = (sum(nEnhancers > 0)/n())*100)

# gene type defined in Williams et al.
p <- ggplot(plotdf,
            aes(x = factor(williams_category,
                           levels = c("Developmental_Non_Pleiotropic",
                                      "Pleiotropic",
                                      "Immune_Non_Pleiotropic")), 
                y = nEnhancers)) +
  geom_jitter(aes(color = williams_category), height = 0) +
  geom_text(data = plotdf_pcts,
            aes(x = williams_category, y = 3, 
                label = paste0("% of genes with at least\none enhancer: ",
                               round(pct_with_enhancer, digits = 2),
                               "%"))) +
  geom_violin() +
  theme_classic() +
  theme(legend.position = "none") +
  ylab("number of enhancers per gene") +
  xlab("gene category") +
  ggtitle("Constitutive S2* enhancers")

p

#### Bar plots of pct genes with 1+ enhancer ####
plotdf_pcts <- genedf |> select(gene_name, williams_category) |> 
  expand_grid(condition = unique(nEnhancerdf$condition)) |> 
  left_join(y = filter(nEnhancerdf, FlyBaseID %in% genedf$gene_name),
            by = c("gene_name"="FlyBaseID", "condition")) |> 
  mutate(nEnhancers = coalesce(nEnhancers, 0)) |> 
  group_by(williams_category, condition) |> 
  summarise(n_with_enhancer = sum(nEnhancers > 0),
            n = n(),
            pct_with_enhancer = (n_with_enhancer/n)*100)

ggplot(plotdf_pcts, aes(x = williams_category, y = pct_with_enhancer/100)) +
  geom_bar(aes(fill = williams_category), stat = "identity") +
  theme_classic() +
  scale_fill_manual(values = c("Developmental_Non_Pleiotropic"="#F8766D",
                               "Pleiotropic"="#619CFF",
                               "Immune_Non_Pleiotropic"="#00BA38",
                               "Other"="grey80")) +
  facet_wrap(~factor(condition, levels = c("dev", "s2_constit", "s2_control",
                                           "20E", "HKSM_only"))) +
  ylab("fraction of genes with at least one enhancer") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Just dev, hksm-only, and constit (most useful comparisons)
p <- ggplot(filter(plotdf_pcts, condition %in% c("s2_constit", "HKSM_only", "dev")),
       aes(x = williams_category, y = pct_with_enhancer/100)) +
  geom_bar(aes(fill = factor(williams_category,
                             levels = c("Other", "Developmental_Non_Pleiotropic",
                                        "Immune_Non_Pleiotropic", "Pleiotropic"))), 
           stat = "identity") +
  theme_classic() +
  scale_fill_manual(values = c("Developmental_Non_Pleiotropic"="#F8766D",
                               "Pleiotropic"="#619CFF",
                               "Immune_Non_Pleiotropic"="#00BA38",
                               "Other"="grey80")) +
  facet_wrap(~factor(condition, levels = c("s2_constit", "HKSM_only", "dev"))) +
  ylab("fraction of genes with at least one enhancer") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none")
p
ggsave(plot = p,
       filename = "figures/enhancer_pcts_bar.png",
       width = 7, height = 5)

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
plotlist <- list(development = plotdf |> 
                   filter(williams_category == "Pleiotropic") |> 
                   filter(nEnhancers_dev > 0) |> 
                   select(gene_name) |> pull(),
                 immune = plotdf |> 
                   filter(williams_category == "Pleiotropic") |> 
                   filter(nEnhancers_hksm > 0) |> 
                   select(gene_name) |> pull())
p <- ggVennDiagram(plotlist, label_color = "black", 
                   label_alpha = 0) + 
  coord_flip() +
  scale_fill_gradient(low="grey90",high = "red")
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
  mutate(nEnhancers_hksm = coalesce(nEnhancers_hksm, 0)) |> 
  mutate(hasEnhancers_dev = if_else(nEnhancers_dev > 0,
                                    true = "1+",
                                    false = "none"),
         hasEnhancers_hksm = if_else(nEnhancers_hksm > 0,
                                     true = "1+",
                                     false = "none")) |> 
  mutate(hasEnhancers_dev = factor(hasEnhancers_dev, levels = c("none", "1+")),
         hasEnhancers_hksm = factor(hasEnhancers_hksm, levels = c("none", "1+")))

# wilcoxon tests of mean differences per group
wilcox.test(signal_to_noise_dev ~ hasEnhancers_dev, 
            data = filter(plotdf, williams_category == "Developmental_Non_Pleiotropic"))
wilcox.test(signal_to_noise_dev ~ hasEnhancers_dev, 
            data = filter(plotdf, williams_category == "Immune_Non_Pleiotropic"))
wilcox.test(signal_to_noise_dev ~ hasEnhancers_dev, 
            data = filter(plotdf, williams_category == "Pleiotropic"))
wilcox.test(signal_to_noise_immune ~ hasEnhancers_hksm, 
            data = filter(plotdf, williams_category == "Developmental_Non_Pleiotropic"))
wilcox.test(signal_to_noise_immune ~ hasEnhancers_hksm, 
            data = filter(plotdf, williams_category == "Immune_Non_Pleiotropic"))
wilcox.test(signal_to_noise_immune ~ hasEnhancers_hksm, 
            data = filter(plotdf, williams_category == "Pleiotropic"))

# nEnhancers in dev, SNR in Embryonic Development
plotdf_counts <- plotdf |>
  group_by(hasEnhancers_dev, williams_category) |> 
  summarise(n = n())
p_dev <- ggplot(plotdf, aes(x = hasEnhancers_dev, 
           y = signal_to_noise_dev)) +
  geom_boxplot(aes(group = hasEnhancers_dev,
                   fill = williams_category)) +
  geom_text(data = plotdf_counts,
            aes(x = hasEnhancers_dev, y = 30, label = n)) +
  facet_wrap(~williams_category) +
  theme_classic() +
  xlab("Gene has at least one enhancer") +
  ggtitle("Emybronic development enhancers")

# nEnhancers in immune, SNR in immune
plotdf_counts <- plotdf |> 
  group_by(hasEnhancers_hksm, williams_category) |> 
  summarise(n = n())
p_hksm <- ggplot(plotdf,
       aes(x = hasEnhancers_hksm, y = signal_to_noise_immune)) +
  geom_boxplot(aes(group = hasEnhancers_hksm,
                   fill = williams_category)) +
  geom_text(data = plotdf_counts,
            aes(x = hasEnhancers_hksm, y = 10, label = n)) +
  facet_wrap(~williams_category) +
  theme_classic() +
  xlab("Gene has at least one enhancer") +
  ggtitle("HKSM-only enhancers")

ggarrange(p_dev, p_hksm, nrow = 2, ncol = 1, common.legend = TRUE)
ggsave(plot = ggarrange(p_dev, p_hksm, nrow = 2, ncol = 1, common.legend = TRUE, legend = "right"),
       filename = "figures/nEnhancers_SNR_box.png",
       width = 10, height = 8)

# what are those immune genes that seem so dev active?
plotgenes <- plotdf |> filter(nEnhancers_dev > 0 & 
                   williams_category == "Immune_Non_Pleiotropic") |> 
  left_join(drosophila_lookup, by = c("gene_name"="FlyBaseID")) |> 
  select(gene_name, common_name, signal_to_noise_dev, signal_to_noise_immune)
scales::hue_pal()(n = 3)[2]
ggarrange(plotlist = plotExpressionProfileAllEnvironments(.gene_idxs = plotgenes$gene_name, 
                                               .gene_names = plotgenes$common_name,
                                               .environments = c("oog", "dev", "immune"),
                                               .color_vec = c("#00BA38", "#00BA38", "#00BA38")),
          nrow = 2, ncol = 2, legend = "none")

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
plotdf_counts <- plotdf |> 
  group_by(hasEnhancers_hksm, williams_category) |> 
  summarise(n = n())
p_hksm <- ggplot(plotdf,
       aes(x = hasEnhancers_hksm, y = log2(mean_immune))) +
  geom_boxplot(aes(group = hasEnhancers,
                   fill = williams_category)) +
  geom_text(data = plotdf_counts,
            aes(x = hasEnhancers, y = 10, label = n)) +
  facet_wrap(~williams_category) +
  theme_classic() +
  xlab("Gene has at least one enhancer") +
  ggtitle("HKSM-only enhancers")
p_hksm # yes

# is that true in development too?
plotdf_counts <- plotdf |> mutate(hasEnhancers = nEnhancers_dev > 0) |> 
  group_by(hasEnhancers, williams_category) |> 
  summarise(n = n())
p_dev <- ggplot(data = mutate(plotdf, hasEnhancers = nEnhancers_dev > 0),
       aes(x = hasEnhancers, y = log2(mean_dev))) +
  geom_boxplot(aes(group = hasEnhancers,
                   fill = williams_category)) +
  geom_text(data = plotdf_counts,
            aes(x = hasEnhancers, y = 14, label = n)) +
  facet_wrap(~williams_category) +
  theme_classic() +
  ggtitle("Embryonic development enhancers") +
  xlab("Gene has at least one enhancer") 
p_dev # not for dev genes!

# wilcoxon tests of mean differences per group
wilcox.test(mean_dev ~ hasEnhancers_dev, 
            data = filter(plotdf, williams_category == "Developmental_Non_Pleiotropic"))
wilcox.test(mean_dev ~ hasEnhancers_dev, 
            data = filter(plotdf, williams_category == "Immune_Non_Pleiotropic"))
wilcox.test(mean_dev ~ hasEnhancers_dev, 
            data = filter(plotdf, williams_category == "Pleiotropic"))
wilcox.test(mean_immune ~ hasEnhancers_hksm, 
            data = filter(plotdf, williams_category == "Developmental_Non_Pleiotropic"))
wilcox.test(mean_immune ~ hasEnhancers_hksm, 
            data = filter(plotdf, williams_category == "Immune_Non_Pleiotropic"))
wilcox.test(mean_immune ~ hasEnhancers_hksm, 
            data = filter(plotdf, williams_category == "Pleiotropic"))

ggsave(plot = ggarrange(p_dev, p_hksm, nrow = 2, ncol = 1, common.legend = TRUE, legend = "right"),
       filename = "figures/nEnhancers_mean_box.png",
       width = 10, height = 8)

#### visualizing genes by SNR ####
# immune genes with enhancers, SNR > 1:
gene_idxs <- plotdf |> filter(nEnhancers_hksm > 0 & 
                                signal_to_noise_immune > 1 &
                                williams_category == "Immune_Non_Pleiotropic") |> 
  select(gene_name) |> pull()
gene_names <- map(gene_idxs, \(g) {drosophila_lookup$Gene_Symbol[which(drosophila_lookup$FlyBaseID == g)]}) |> 
  unlist()
plotExpressionProfile(.counts = counts$immune, .info = infodf$immune,
                      .gene_idxs = gene_idxs, .gene_names = gene_names)
# SNR <= 1:
gene_idxs <- plotdf |> filter(nEnhancers_hksm > 0 & 
                                signal_to_noise_immune <= 1 &
                                williams_category == "Immune_Non_Pleiotropic") |> 
  select(gene_name) |> pull()
gene_names <- map(gene_idxs, \(g) {drosophila_lookup$Gene_Symbol[which(drosophila_lookup$FlyBaseID == g)]}) |> 
  unlist()
plotExpressionProfile(.counts = counts$immune, .info = infodf$immune,
                      .gene_idxs = gene_idxs, .gene_names = gene_names)
# that checks out

# conclusion: these don't have particularly high SNR or any uniform behavior. 
# The PGRPs are noisy between replicates. The most interesting genes are 3/4 of the
# un-named ones: CG16985 and CG16986 both increase steadily in both replicates and 
# CG17271 has a robust spike fairly early

#### Inspecting expression of immune genes with enhancers ####
# Immune genes with > 1 hksm enhancer
plotdf |> filter(nEnhancers_hksm > 1) |> 
  left_join(drosophila_lookup, by = c("gene_name"="FlyBaseID")) |> 
  select(gene_name, Gene_Symbol, nEnhancers_hksm, nEnhancers_dev)
# Eip75B
plotExpressionProfileAllEnvironments(.gene_idxs = "FBgn0000568",
                                     .gene_names = "Eip75B",
                                     .environments = "immune")

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

#### Pure proximity enhancer-gene assignments ####
# all Cohen enhancers, regardless of ATACseq accessibility or gene assignment
all_enhancers_cohen <- read_csv("data/EG_TFBS_cat_01_19_activityclass.csv") |> 
  select(Enhancer) |> pull()

# quick function to remove enhancer strand information, as this is the thing
# that might be different between closest enhancer assignments and Cohen enhancer assignments
destrandEnhancer <- function(.enhancer_name) {
  .enhancer_name |> 
  gsub(pattern = "\\(\\+\\)", replacement = "") |> 
    gsub(pattern = "\\(\\-\\)", replacement = "")
}
# tests for destrandEnhancer
destrandEnhancer("2L:119700-121346(+)")
destrandEnhancer("X:119700-121346(-)")

# function to load bedtools closest output into a tibble
# .col_select is used to specify which columns in bed file correspond to which
# feature. There MUST be these 9 features in colnames:
# c("Enhancer_Chromosome","Enhancer_Start", "Enhancer_End",
# "Gene_Chromosome", "Gene_Start", "Gene_End",
# "Gene_Strand", "FlyBaseID", "Distance"
bed2tibble <- function(.file_name, .file_dir = "data/", 
                       .col_select, .col_names = c("Enhancer_Chromosome",
                                                   "Enhancer_Start",
                                                   "Enhancer_End",
                                                   "Gene_Chromosome",
                                                   "Gene_Start",
                                                   "Gene_End",
                                                   "Gene_Strand",
                                                   "FlyBaseID",
                                                   "Distance")) {
  stopifnot(length(.col_select) == length(.col_names))
  stopifnot(all(c("Enhancer_Chromosome",
                 "Enhancer_Start",
                 "Enhancer_End",
                 "Gene_Chromosome",
                 "Gene_Start",
                 "Gene_End",
                 "Gene_Strand",
                 "FlyBaseID",
                 "Distance") %in% .col_names))
  outdf <- read_delim(file = paste0(.file_dir, .file_name),
                      delim = "\t", col_names = FALSE,
                      col_select = all_of(.col_select),
                      show_col_types = FALSE)
  colnames(outdf) <- .col_names
  # bedtools closest output occasionally outputs the same enhancer-gene assignment in
  # multiple rows, possibly b/c of input dmel_genes.bed file having duplicate gene entries
  duplicated_rows <- duplicated(outdf)
  cat("removing", sum(duplicated_rows), "duplicated rows")
  outdf <- outdf[!duplicated_rows,]
  # naming enhancers like Cohen et al.
  outdf$Enhancer <- map(c(1:nrow(outdf)), \(i) {
    min_pos <- min(outdf$Enhancer_Start[i],
                   outdf$Enhancer_End[i])
    max_pos <- max(outdf$Enhancer_Start[i],
                   outdf$Enhancer_End[i])
    out_string <- paste0(gsub("chr", "", outdf$Enhancer_Chromosome[i]),
                         ":", min_pos, "-", max_pos, 
                         "(+)") # Note: this sets enhancer strand to always be +, same as Cohen et al. enhancer-gene assignments
    return(out_string)
  }) |> unlist()
  # checking a few things before outputting
  stopifnot(all(outdf$Enhancer_Chromosome == outdf$Gene_Chromosome)) # genes can't be on different chromosomes
  # if there are 2+ entries for an enhancer, all the gene distances must be equal, 
  # that is true ties (all genes overlapping the feature will have a Distance of 0, 
  # which counts as a true tie)
  is_fake_tie <- map2(c(1:nrow(outdf)),
                      duplicated(outdf$Enhancer), 
                      \(idx, is_dupe) {
                        if (is_dupe) {
                          enhancer_name <- outdf$Enhancer[idx]
                          distances <- outdf |> filter(Enhancer == enhancer_name) |> 
                            select(Distance) |> pull()
                          return(length(unique(distances)) != 1)
                        }
                        else {
                          return(FALSE)
                        }
  }) |> unlist()
  stopifnot(all(!is_fake_tie))
  return(select(outdf, Enhancer, Enhancer_Chromosome,
                Enhancer_Start, Enhancer_End, FlyBaseID,
                Gene_Start, Gene_End, Gene_Strand))
}
# loading alternative enhancer-gene assignments to compare
# to Cohen et al.

# closest single gene using bedtools closest
enhancerdf_closestHKSM <- bed2tibble("HKSM_ONLY_closest.bed",
                                     .col_select = c(1:6, 9, 7, 14))
enhancerdf_closest20E <- bed2tibble("20E_closest.bed",
                                    .col_select = c(1:6, 9, 7, 14))
enhancerdf_closest_s2control <- bed2tibble("s2_control_closest.bed",
                                           .col_select = c(1:6, 9, 7, 14))
enhancerdf_closest_s2constit <- bed2tibble("s2_constit_closest.bed",
                                           .col_select = c(1:6, 9, 7, 14))

# # Kvon et al. 2014 enhancers
# # writing to BED file for use in bedtools
# out_bed_Kvon <- select(enhancerdf_Kvon, Chrom_enhancer, Start_enhancer, End_enhancer)
# names(out_bed_Kvon) <- NULL
# write_delim(out_bed_Kvon, col_names = FALSE, file = "data/kvon_strong_dev_enhancers_dm3.bed")
# # Kvon enhancers are dm3, so using USCS liftOver tool to get dm6 coordinates
enhancerdf_closest_dev <- bed2tibble("dev_closest.bed", .col_select = c(1:6, 9, 7, 14))
all_enhancers_kvon <- map(c(1:nrow(enhancerdf_Kvon)), \(i) {
  min_pos <- min(enhancerdf_Kvon$Start_enhancer[i],
                 enhancerdf_Kvon$End_enhancer[i])
  max_pos <- max(enhancerdf_Kvon$Start_enhancer[i],
                 enhancerdf_Kvon$End_enhancer[i])
  out_string <- paste0(gsub("chr", "", enhancerdf_Kvon$Chrom_enhancer[i]),
                       ":", min_pos, "-", max_pos, 
                       "(+)") # Note: this sets enhancer strand to always be +, same as Cohen et al. enhancer-gene assignments
  return(out_string)
}) |> unlist() |> unique()

# how many enhancers in closest are in Cohen et al.?
sum(enhancerdf_closestHKSM$Enhancer %in% all_enhancers_cohen) # should be nrow(enhancerdf_closestHKSM)
sum(!(enhancerdf_closestHKSM$Enhancer %in% all_enhancers_cohen)) # should be 0
sum(!(enhancerdf_closest20E$Enhancer %in% all_enhancers_cohen)) # should be 0
sum(!(enhancerdf_closest_s2control$Enhancer %in% all_enhancers_cohen)) # should be 0
sum(!(enhancerdf_closest_s2constit$Enhancer %in% all_enhancers_cohen)) # should be 0

# TODO: 201 in Kvon 2014, and 575 not, are these dm3 -> dm6 issues?
sum((enhancerdf_closest_dev$Enhancer %in% all_enhancers_kvon))
sum(!(enhancerdf_closest_dev$Enhancer %in% all_enhancers_kvon)) # should be 0
setdiff(all_enhancers_kvon, enhancerdf_closest_dev$Enhancer)[1:10]

# counting nEnhancers per gene
# HKSM-only
nEnhancerdf_closestHKSM <- enhancerdf_closestHKSM |> 
  group_by(FlyBaseID) |> 
  summarise(nEnhancers = n()) |>
  mutate(condition = "closest_HKSM_only")
table(nEnhancerdf_closestHKSM$nEnhancers)
# 20E
nEnhancerdf_closest20E <- enhancerdf_closest20E |> 
  group_by(FlyBaseID) |> 
  summarise(nEnhancers = n()) |>
  mutate(condition = "closest_20E")
table(nEnhancerdf_closest20E$nEnhancers)
# Control
nEnhancerdf_closest_s2_control <- enhancerdf_closest_s2control |> 
  group_by(FlyBaseID) |> 
  summarise(nEnhancers = n()) |>
  mutate(condition = "closest_s2_control")
table(nEnhancerdf_closest_s2_control$nEnhancers)
# Constit
nEnhancerdf_closest_s2_constit <- enhancerdf_closest_s2constit |> 
  group_by(FlyBaseID) |> 
  summarise(nEnhancers = n()) |>
  mutate(condition = "closest_s2_constit")
table(nEnhancerdf_closest_s2_constit$nEnhancers)
# Dev
nEnhancerdf_closest_dev <- enhancerdf_closest_dev |> 
  group_by(FlyBaseID) |> 
  summarise(nEnhancers = n()) |>
  mutate(condition = "closest_dev")
table(nEnhancerdf_closest_dev$nEnhancers)

# QC_nEnhancerdf for all conditions, comparing published assignments to closest
QC_nEnhancerdf <- bind_rows(nEnhancerdf_closestHKSM,
                         nEnhancerdf_closest20E,
                         nEnhancerdf_closest_s2_control,
                         nEnhancerdf_closest_s2_constit,
                         nEnhancerdf_closest_dev)
# plotting
plotdf_pcts <- genedf |> select(gene_name, williams_category) |> 
  expand_grid(condition = unique(QC_nEnhancerdf$condition)) |> 
  left_join(y = filter(QC_nEnhancerdf, FlyBaseID %in% genedf$gene_name),
            by = c("gene_name"="FlyBaseID", "condition")) |>
  mutate(nEnhancers = coalesce(nEnhancers, 0)) |> 
  group_by(williams_category, condition) |> 
  summarise(n_with_enhancer = sum(nEnhancers > 0),
            n = n(),
            pct_with_enhancer = (n_with_enhancer/n)*100)

ggplot(filter(plotdf_pcts, condition %in% c("closest_s2_constit", "closest_HKSM_only", 
                                            "closest_dev")), 
       aes(x = williams_category, y = pct_with_enhancer/100)) +
  geom_bar(aes(fill = williams_category), stat = "identity") +
  theme_classic() +
  facet_wrap(~factor(condition, levels = c("closest_s2_constit", 
                                           "closest_HKSM_only",
                                           "closest_dev"))) +
  ylab("fraction of genes with at least one enhancer") +
  xlab("") +
  scale_fill_manual(values = c("Developmental_Non_Pleiotropic"="#F8766D",
                               "Pleiotropic"="#619CFF",
                               "Immune_Non_Pleiotropic"="#00BA38",
                               "Other"="grey80")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none") +
  ggtitle("Bedtools closest Enhancer-Gene assignments")

#### Are pleiotropic genes more likely to be near the middle of each chromosome leg? ####
# We've seen that plei genes have more enhancers across the board, also when done
# using "pure proximity" assignments (bedtools closest). Is this driven by their
# positions on each chromosome?

# The goal: Geom_density where x axis is chrom position ~ facet_wrap by chromosome/leg
# To properly integrate density and count each gene equally, we'll take each gene's midpoint
# as it's location (we could also check if Plei genes are longer)

# These might be useful regex to know:
# enhancer start
str_extract("2L:16836-18924(+)", pattern = "(?<=:)\\d+(?=-)") # matches any # consecutive digits between : and -
# enhancer end
str_extract("2L:16836-18924(+)", pattern = "(?<=\\-)\\d+") # matches any # consecutive digits after -
# chromsome
str_extract("2L:16836-18924(+)", pattern = "^.*(?=:)") # matches any characters from start until :

# But we might not need those b/c dmel genes already have Chrom, start, and end parsed
dmel_genes <- read_delim("data/dmel-all-r6.67.gtf", delim = "\t", col_names = FALSE)
colnames(dmel_genes) <- c("Chromosome", "Database", "Feature", "Start", "End", "Who_knows",
                          "Strand", "Who_knows2", "Information")
dmel_genes <- dmel_genes |> filter(Feature == "gene") |> 
  mutate(FlyBaseID = str_remove_all(str_extract(Information,
                                                pattern = '(?<=gene_id\\s).+?(?=;)'),
                                    pattern = '\"')) |> 
  select(FlyBaseID, Chromosome, Start, End, Strand) |> 
  arrange(Chromosome, Start)

chromdf <- tibble(Chromosome = c("2L", "2R", "3L", "3R", "4", "X", "Y"),
                  Chrom_Length = c(23513712, 25286936, 28110227, 32079331, 
                                   1348131, 23542271, 3667352))
plotdf <- dmel_genes |> filter(Chromosome %in% chromdf$Chromosome) |> 
  left_join(chromdf, by = "Chromosome") |> 
  inner_join(genedf, by = c("FlyBaseID"="gene_name")) |> 
  mutate(midpoint = (Start + End)/2,
         midpoint_norm = midpoint/Chrom_Length)

# counts
p <- ggplot(plotdf, aes(x = midpoint_norm)) +
  geom_histogram(aes(fill = williams_category)) +
  facet_wrap(~Chromosome) +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("Percent along chromosome (midpoint of gene)") +
  scale_fill_manual(values = c("Developmental_Non_Pleiotropic"="#F8766D",
                               "Pleiotropic"="#619CFF",
                               "Immune_Non_Pleiotropic"="#00BA38",
                               "Other"="grey80"))
p
ggsave(plot = p, filename = "figures/chrom_position_hist.pdf",
       width = 7, height = 5)

# density
p <- ggplot(plotdf, aes(x = midpoint_norm)) +
  geom_density(aes(fill = williams_category), alpha = 0.75) +
  facet_wrap(~Chromosome, scales = "free_y") +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("Percent along chromosome (midpoint of gene)") +
  scale_fill_manual(values = c("Developmental_Non_Pleiotropic"="#F8766D",
                               "Pleiotropic"="#619CFF",
                               "Immune_Non_Pleiotropic"="#00BA38",
                               "Other"="grey80"))
p
ggsave(plot = p, filename = "figures/chrom_position_density.pdf",
       width = 7, height = 5)

# Nothing obvious

### Are plei genes longer?
plotdf$gene_length <- abs(plotdf$Start - plotdf$End)

# length distribution per chromosome
p <- ggplot(plotdf, aes(x = log2(gene_length))) +
  geom_density(aes(fill = williams_category), alpha = 0.75) +
  facet_wrap(~Chromosome, scales = "free_y") +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("Gene length (log2(bp))") +
  scale_fill_manual(values = c("Developmental_Non_Pleiotropic"="#F8766D",
                               "Pleiotropic"="#619CFF",
                               "Immune_Non_Pleiotropic"="#00BA38",
                               "Other"="grey80"))
p
ggsave(plot = p, filename = "figures/gene_length.pdf",
       width = 7, height = 5)

# boxplot comparison of lengths
p <- ggplot(plotdf, aes(y = log2(gene_length), 
                        x = factor(williams_category,
                                   levels = c("Other",
                                              "Immune_Non_Pleiotropic",
                                              "Developmental_Non_Pleiotropic",
                                              "Pleiotropic")))) +
  geom_boxplot(aes(fill = williams_category)) +
  #facet_wrap(~Chromosome, scales = "free_y") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("Gene length (log2(bp))") +
  stat_compare_means(comparisons = list(c("Other", "Immune_Non_Pleiotropic"),
                                        c("Other", "Developmental_Non_Pleiotropic"),
                                        c("Other", "Pleiotropic"))) +
  scale_fill_manual(values = c("Developmental_Non_Pleiotropic"="#F8766D",
                               "Pleiotropic"="#619CFF",
                               "Immune_Non_Pleiotropic"="#00BA38",
                               "Other"="grey80"))
p
ggsave(plot = p, filename = "figures/gene_length_box.pdf",
       width = 3, height = 6)

# Yup they are longer, so are dev genes

#### Enhancer-gene assignments using gene midpoints ####
# closest single gene using bedtools closest, genes are only 1bp long at their midpoint
enhancerdf_midpointsHKSM <- bed2tibble("HKSM_ONLY_midpoints.bed",
                                       .col_select = c(1:6, 8, 7, 9))
enhancerdf_midpoints20E <- bed2tibble("20E_midpoints.bed",
                                      .col_select = c(1:6, 8, 7, 9))
enhancerdf_midpoints_s2control <- bed2tibble("s2_control_midpoints.bed",
                                             .col_select = c(1:6, 8, 7, 9))
enhancerdf_midpoints_s2constit <- bed2tibble("s2_constit_midpoints.bed",
                                             .col_select = c(1:6, 8, 7, 9))

# # Kvon et al. 2014 enhancers
# # writing to BED file for use in bedtools
# out_bed_Kvon <- select(enhancerdf_Kvon, Chrom_enhancer, Start_enhancer, End_enhancer)
# names(out_bed_Kvon) <- NULL
# write_delim(out_bed_Kvon, col_names = FALSE, file = "data/kvon_strong_dev_enhancers_dm3.bed")
# # Kvon enhancers are dm3, so using USCS liftOver tool to get dm6 coordinates
enhancerdf_midpoints_dev <- bed2tibble("dev_midpoints.bed",
                                       .col_select = c(1:6, 8, 7, 9))
all_enhancers_kvon <- map(c(1:nrow(enhancerdf_Kvon)), \(i) {
  min_pos <- min(enhancerdf_Kvon$Start_enhancer[i],
                 enhancerdf_Kvon$End_enhancer[i])
  max_pos <- max(enhancerdf_Kvon$Start_enhancer[i],
                 enhancerdf_Kvon$End_enhancer[i])
  out_string <- paste0(gsub("chr", "", enhancerdf_Kvon$Chrom_enhancer[i]),
                       ":", min_pos, "-", max_pos, 
                       "(+)") # Note: this sets enhancer strand to always be +, same as Cohen et al. enhancer-gene assignments
  return(out_string)
}) |> unlist() |> unique()

# how many enhancers in midpoints are in Cohen et al.?
sum(enhancerdf_midpointsHKSM$Enhancer %in% all_enhancers_cohen) # should be nrow(enhancerdf_midpointsHKSM)
sum(!(enhancerdf_midpointsHKSM$Enhancer %in% all_enhancers_cohen)) # should be 0
sum(!(enhancerdf_midpoints20E$Enhancer %in% all_enhancers_cohen)) # should be 0
sum(!(enhancerdf_midpoints_s2control$Enhancer %in% all_enhancers_cohen)) # should be 0
sum(!(enhancerdf_midpoints_s2constit$Enhancer %in% all_enhancers_cohen)) # should be 0

# TODO: 201 in Kvon 2014, and 575 not, are these dm3 -> dm6 issues?
sum((enhancerdf_midpoints_dev$Enhancer %in% all_enhancers_kvon))
sum(!(enhancerdf_midpoints_dev$Enhancer %in% all_enhancers_kvon)) # should be 0
setdiff(all_enhancers_kvon, enhancerdf_midpoints_dev$Enhancer)[1:10]

# counting nEnhancers per gene
# HKSM-only
nEnhancerdf_midpointsHKSM <- enhancerdf_midpointsHKSM |> 
  group_by(FlyBaseID) |> 
  summarise(nEnhancers = n()) |>
  mutate(condition = "midpoints_HKSM_only")
table(nEnhancerdf_midpointsHKSM$nEnhancers)
# 20E
nEnhancerdf_midpoints20E <- enhancerdf_midpoints20E |> 
  group_by(FlyBaseID) |> 
  summarise(nEnhancers = n()) |>
  mutate(condition = "midpoints_20E")
table(nEnhancerdf_midpoints20E$nEnhancers)
# Control
nEnhancerdf_midpoints_s2_control <- enhancerdf_midpoints_s2control |> 
  group_by(FlyBaseID) |> 
  summarise(nEnhancers = n()) |>
  mutate(condition = "midpoints_s2_control")
table(nEnhancerdf_midpoints_s2_control$nEnhancers)
# Constit
nEnhancerdf_midpoints_s2_constit <- enhancerdf_midpoints_s2constit |> 
  group_by(FlyBaseID) |> 
  summarise(nEnhancers = n()) |>
  mutate(condition = "midpoints_s2_constit")
table(nEnhancerdf_midpoints_s2_constit$nEnhancers)
# Dev
nEnhancerdf_midpoints_dev <- enhancerdf_midpoints_dev |> 
  group_by(FlyBaseID) |> 
  summarise(nEnhancers = n()) |>
  mutate(condition = "midpoints_dev")
table(nEnhancerdf_midpoints_dev$nEnhancers)

# QC_nEnhancerdf for all conditions, comparing published assignments to midpoints
QC_nEnhancerdf <- bind_rows(nEnhancerdf_midpointsHKSM,
                            nEnhancerdf_midpoints20E,
                            nEnhancerdf_midpoints_s2_control,
                            nEnhancerdf_midpoints_s2_constit,
                            nEnhancerdf_midpoints_dev)
# plotting
plotdf_pcts <- genedf |> select(gene_name, williams_category) |> 
  expand_grid(condition = unique(QC_nEnhancerdf$condition)) |> 
  left_join(y = filter(QC_nEnhancerdf, FlyBaseID %in% genedf$gene_name),
            by = c("gene_name"="FlyBaseID", "condition")) |>
  mutate(nEnhancers = coalesce(nEnhancers, 0)) |> 
  group_by(williams_category, condition) |> 
  summarise(n_with_enhancer = sum(nEnhancers > 0),
            n = n(),
            pct_with_enhancer = (n_with_enhancer/n)*100)

ggplot(filter(plotdf_pcts, condition %in% c("midpoints_s2_constit",
                                            "midpoints_HKSM_only",
                                            "midpoints_dev")), 
       aes(x = williams_category, y = pct_with_enhancer/100)) +
  geom_bar(aes(fill = williams_category), stat = "identity") +
  theme_classic() +
  facet_wrap(~factor(condition, levels = c("midpoints_s2_constit",
                                           "midpoints_HKSM_only",
                                           "midpoints_dev"))) +
  ylab("fraction of genes with at least one enhancer") +
  xlab("") +
  scale_fill_manual(values = c("Developmental_Non_Pleiotropic"="#F8766D",
                               "Pleiotropic"="#619CFF",
                               "Immune_Non_Pleiotropic"="#00BA38",
                               "Other"="grey80")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none") +
  ggtitle("Bedtools midpoints Enhancer-Gene assignments")

### Are any chromosomes enriched for Plei gene enhancer assignments?
enhancerdf_closest_s2constit |> left_join(y = select(genedf, gene_name, williams_category),
                                          by = c("FlyBaseID"="gene_name")) |> 
  select(williams_category, Enhancer_Chromosome) |> table()
# 3R perhaps?
enhancerdf_closest_s2constit |> left_join(y = select(genedf, gene_name, williams_category),
                                          by = c("FlyBaseID"="gene_name")) |> 
  filter(williams_category == "Pleiotropic" & Enhancer_Chromosome == "chr3R") |> 
  left_join(y = drosophila_lookup, by = "FlyBaseID") |> 
  select(Gene_Symbol, Summary) |> View() # puckered, pointed, kayak, foxo, p53

#### Are the pleiotropic genes with s2 constitutivly-active enhancers hematopoesis genes? ####

# TODO: Both of these should be filtering for the pleiotropic genes with 
# s2* constitutive enhancers, but taking from nEnhancers rather than enhancer_cohen
# results in more genes
more_s2constit_plei_genes <- nEnhancerdf_constit |> 
  left_join(genedf, by = c("FlyBaseID"="gene_name")) |> 
  filter(williams_category == "Pleiotropic") |> 
  left_join(drosophila_lookup, by = "FlyBaseID") |> 
  select(Gene_Symbol) |> pull() |> unique()

fewer_s2constit_plei_genes <- enhancerdf_Cohen |> 
  filter(Group == "Constit" & Accessibility != "Always closed") |> 
  left_join(genedf, by = c("Gene"="gene_name")) |> 
  filter(williams_category == "Pleiotropic") |> 
  left_join(drosophila_lookup, by = c("Gene"="FlyBaseID")) |> 
  mutate(Enhancer_Chrom = str_extract(Enhancer, pattern = "^.*(?=:)")) |> 
  group_by(Gene_Symbol) |> 
  mutate(nEnhancers = n()) |> 
  select(Gene_Symbol) |> pull() |> unique()

setdiff(more_s2constit_plei_genes, fewer_s2constit_plei_genes)

