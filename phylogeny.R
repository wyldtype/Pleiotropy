sapply(c("tidyr", "dplyr", "ggplot2", "ggpubr", "ggalluvial", "purrr", "readxl", "readr", "stringr", "plotly", "lmtest", "pheatmap"), FUN = require,
       character.only = TRUE)

# loading counts and functions
source("utils.r")
load("data/CountsList.RData")
load("data/SignalToNoise.RData")
load("data/Enhancers.RData")

# loading HOG data
dmel <- read_tsv("../phylogeny/215Fly_Orthogroups_tables/Dmel_HOG_association.tsv")
HOG_to_OG <- read_tsv("../phylogeny/215Fly_Orthogroups_tables/HOG_to_OG_table.tsv")
ortho <- read_tsv("../phylogeny/215Fly_Orthogroups_tables/Orthogroups.tsv")
ortho_GeneCounts <- read_tsv("../phylogeny/215Fly_Orthogroups_tables/Orthogroups.GeneCount.tsv")

#### Validating/exploring orthogroup datasets ####
# duplications that aren't 
# found across all species I think are rows in ortho
# where there are NA values for certain species
# Why do I think that? 93% of NAs in Dmel or Dsim are NAs in both:
sum(is.na(ortho$DROSOPHILA_MELANOGASTER) & 
      is.na(ortho$DROSOPHILA_SIMULANS))/sum(is.na(ortho$DROSOPHILA_MELANOGASTER) |
                                              is.na(ortho$DROSOPHILA_SIMULANS))

### Polarizing duplication history by Dmel
# Which duplications ocurred before Drosophilidae split?
OGs_in_all <- ortho$Orthogroup[apply(ortho[,-1], MARGIN=1, \(x) sum(!is.na(x))/length(x) >= 0.99)] # "universal" defined by Dhakad et al. 2026 as in >=99% species
length(OGs_in_all)
length(ortho$Orthogroup)

# Which duplications are specific to Dmel? 
OGs_in_dmel <- ortho$Orthogroup[sapply(ortho$DROSOPHILA_MELANOGASTER,\(x) !is.na(x))] |> 
  setdiff(y = OGs_in_all)

# Which are duplicated in other lineages but missing from the Dmel lineage?
OGs_not_in_dmel <- ortho$Orthogroup[sapply(ortho$DROSOPHILA_MELANOGASTER,\(x) is.na(x))] |> 
  setdiff(y = OGs_in_all)

# These should now be 3 mutually exclusive sets:
sum(length(OGs_in_all), length(OGs_in_dmel), length(OGs_not_in_dmel)) # should be same as next line
length(ortho$Orthogroup) # should be same as previous line
setdiff(x = unique(ortho$Orthogroup), y = unique(c(OGs_in_all, OGs_in_dmel, OGs_not_in_dmel))) # should be empty
setdiff(y = unique(ortho$Orthogroup), x = unique(c(OGs_in_all, OGs_in_dmel, OGs_not_in_dmel))) # should be empty

#### Assigning Dmel gene names to HOGs ####
# Helper function for creating HOGdf, summarises williams et al. annotations 
detectWilliams <- function(.williams_categories) {
  filtered_categories <- .williams_categories[.williams_categories %in% c("Immune_Non_Pleiotropic",
                                                                          "Developmental_Non_Pleiotropic",
                                                                          "Pleiotropic",
                                                                          "Mixed")] |> unique()
  if (length(filtered_categories) == 1) {
    return(filtered_categories)
  }
  if (length(filtered_categories) == 0) {
    return("None")
  }
  if (length(filtered_categories) > 1) {
    return("Mixed")
  }
}
# tests for detectWilliams
test_cats <- c("Other", "Other", "Other") # all "Other"
detectWilliams(test_cats)
test_cats <- c("Other", "Immune_Non_Pleiotropic", "Other", NA,
               "Immune_Non_Pleiotropic", "Other") # some "Immune"
detectWilliams(test_cats)
test_cats <- c("Other", "Immune_Non_Pleiotropic", "Other", NA,
               "Mixed", "Other") # one "Mixed"
detectWilliams(test_cats)

### creating HOGdf, a dataframe with one row per HOG
HOGdf <- HOG_to_OG |> select(HOG, OG) |> # Note: the `Gene Tree Parent Clade` column seemed initially promising, but the documentation is lacking and the nXXX number does not appear to reflect the branching order (there are 215 species and 1400+ clades), so we'll probably use the species annotation files to figure out exactly where the duplication occurred
  mutate(OG_cons = if_else(OG %in% OGs_in_all,
                           true = "all", false = if_else(OG %in% OGs_in_dmel,
                                                         true = "dmel", false = "notDmel"))) |> 
  left_join(y = select(dmel, HOG, Dmel_FBgnID, Dmel_CommonNames),
            by = "HOG") |> 
  mutate(Dmel_FBgnID = coalesce(Dmel_FBgnID, "none")) |> 
  separate_longer_delim(col = c("Dmel_FBgnID", "Dmel_CommonNames"), delim = ",") |> 
  left_join(select(genedf, gene_name, williams_category), # We can add SNR etc later if that is useful
            by = c("Dmel_FBgnID"="gene_name")) |> 
  group_by(HOG, OG, OG_cons) |> 
  summarise(williams_HOG = detectWilliams(williams_category),
            Dmel_FBgnID = paste(unique(Dmel_FBgnID), collapse = ","), # Note: sep = "," does not collapse concatenated (c()) strings into one!
            Dmel_CommonNames = paste(unique(Dmel_CommonNames), collapse = ",")) |> 
  ungroup() |> 
  group_by(OG) |> 
  mutate(williams_OG = detectWilliams(williams_HOG),
         nOG = n(),
         nOG_dmel = sum(Dmel_FBgnID != "none")) |> 
  ungroup()

# QC: We shouldn't have lost any HOGs:
nrow(HOG_to_OG)
length(unique(HOG_to_OG$HOG))
nrow(HOGdf)
length(unique(HOGdf$HOG))

#### Saving HOGdf ####
save(HOGdf, file = "data/HOG.RData")

### How many Dmel genes have duplications in other lineages?

# 12,145 HOGs in the 215 Drosophila dataset should have a dmel annotation:
sum(is.na(dmel$HOG))
dmel$HOG |> unique() |> length()

# 26,046 total HOGs:
sum(is.na(HOG_to_OG$HOG))
HOG_to_OG$HOG |> unique() |> length()
HOGdf$HOG |> unique() |>  length()

# 26,046 - 12,145 = 13,901 HOGs that don't contain a Dmel gene
setdiff(HOG_to_OG$HOG, dmel$HOG) |> length()
HOGdf |> filter(Dmel_FBgnID == "none") |> nrow()

# Mmp1 example:
HOGdf |> filter(OG == "OG0003666")
# While the orthogroup is conserved, there are multiple HOGs. So this counts as a plei gene that has one
# lineage-specific duplicate
# Jonah example:
HOGdf |> filter(OG == "OG0000000") |> arrange(HOG)
# most of these HOGs are in Dmel, so they aren't lineage-specific duplicates, just a big family
HOGdf |> filter(OG == "OG0000000") |> mutate(inDmel = !is.na(Dmel_FBgnID)) |> 
  select(HOG, inDmel) |> unique() |> select(inDmel) |> table() # 15 HOGs duplicated elsewhere
# note how if we just took every OG entry we would over-count HOGs  that
# correspond to duplicates with similar sequence (ex Jon65Ai and Jon65Aii) 
# because OGs are only based on all-against-all blast, not genomic position
HOGdf |> filter(OG == "OG0000000") |> mutate(inDmel = !is.na(Dmel_FBgnID)) |>
  select(inDmel) |> table()
# Note that the HOG assignments aren't always accurate to genomic location:
HOGdf |> filter(HOG == "N1.HOG0000019") # Jon65Ai/Jon65Aii are next to Jon65Aii on chrom 3L, Jon25Bii is on a separate chromosome (2L)
dmel |> filter(HOG == "N1.HOG0000019")

# so for fairness, we'll collapse these cases in dmel b/c we don't have the gene-level resolution in non-mel
# but b/c we want accurate Immune/Plei/Dev categories, we'll take a HOG to be "Immune" if any of the 
# dmel genes in it are "Immune_Non_Pleiotropic". If any are Plei or Dev, we'll call it "Mixed". That shouldn't happen very often
plotdf <- HOGdf |> group_by(OG) |> 
  mutate(nHOGs = n(),
         pct_inDmel = sum(Dmel_FBgnID != "none")/nHOGs)

# Orthogroups where D.mel has at least one but not 100% of HOGs
plotdf_not100pct <- filter(plotdf, williams_HOG %in% c("Developmental_Non_Pleiotropic", 
                                                       "Immune_Non_Pleiotropic",
                                                       "Pleiotropic") &
                             pct_inDmel != 1) |> 
  mutate(williams_HOG = factor(williams_HOG, levels = c("Immune_Non_Pleiotropic",
                                                        "Developmental_Non_Pleiotropic", 
                                                        "Pleiotropic")))
p <- ggplot(plotdf_not100pct,
       aes(x = williams_HOG,
           y = pct_inDmel)) +
  geom_violin(aes(fill = williams_HOG)) +
  geom_point(data = summarise(group_by(plotdf_not100pct, williams_HOG), mean_pct = mean(pct_inDmel)),
             aes(x = williams_HOG, y = mean_pct)) +
  geom_text(data = summarise(group_by(plotdf_not100pct, williams_HOG), n = n()),
             aes(x = williams_HOG, y = 1.1, label = paste0("n = ", n, "\northogroups"))) +
  scale_fill_manual(values = c("Developmental_Non_Pleiotropic"="#F8766D",
                               "Pleiotropic"="#619CFF",
                               "Immune_Non_Pleiotropic"="#00BA38")) +
  stat_compare_means(comparisons = list(c("Immune_Non_Pleiotropic", "Developmental_Non_Pleiotropic"),
                                        c("Developmental_Non_Pleiotropic", "Pleiotropic"),
                                        c("Immune_Non_Pleiotropic", "Pleiotropic")),
                     label.y = c(1.2, 1.3, 1.4)) +
  scale_y_continuous(labels = \(x) ifelse(x > 1, "", x), limits = c(0, 1.5)) +
  #geom_jitter(aes(color = williams_HOG, size = (n_inDmel + n_not_inDmel))) +
  ylab("percent of HOGs in Dmel") +
  xlab("") +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Gene family expansions outside of D.mel:\nOrthogroups with partial D. mel membership")
p
ggsave(filename = "figures/violin_ortho_partial.pdf", width = 7, height = 5)

# Orthogroups where D.mel has 100% of HOGs and other lineages have lost the Orthogroup
plotdf_100pct <- filter(plotdf, williams_HOG %in% c("Developmental_Non_Pleiotropic", 
                                                       "Immune_Non_Pleiotropic",
                                                       "Pleiotropic") &
                             pct_inDmel == 1) |> 
  mutate(williams_HOG = factor(williams_HOG, levels = c("Immune_Non_Pleiotropic",
                                                        "Developmental_Non_Pleiotropic", 
                                                        "Pleiotropic")),
         universality = if_else(OG_cons == "all",
                                true = "universal",
                                false = "non-universal")) |> 
  group_by(williams_HOG, universality) |> 
  summarise(n = n())

p <- ggplot(plotdf_100pct,
            aes(x = williams_HOG, y = n)) +
  geom_bar(aes(fill = williams_HOG), stat = "identity") +
  geom_text(aes(x = williams_HOG, y = 1.1, label = paste0("n = ", n, "\northogroups"))) +
  scale_fill_manual(values = c("Developmental_Non_Pleiotropic"="#F8766D",
                               "Pleiotropic"="#619CFF",
                               "Immune_Non_Pleiotropic"="#00BA38")) +
  facet_wrap(~universality) +
  ylab("number of orthogroups") +
  xlab("") +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Gene family losses outside of D.mel:\nOrthogroups with full D. mel membership")
p

#### Exploring specific HOG family expansions/losses ####