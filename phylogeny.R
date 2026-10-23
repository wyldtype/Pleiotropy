sapply(c("tidyr", "dplyr", "ggplot2", "ggpubr", "ggalluvial", "purrr", "readxl", "readr", "stringr", "plotly", "lmtest", "pheatmap"), FUN = require,
       character.only = TRUE)

# loading counts and functions
source("utils.r")
source("global.r")
load("data/CountsList.RData")
load("data/SignalToNoise.RData")
load("data/Enhancers.RData")

### Loading files

# Dmel HOG associations
dmel <- read_tsv("../phylogeny/215Fly_Orthogroups_tables_v5/Dmel_HOG_association.tsv")

# representation of each HOG in each species
HOG_cons <- read_tsv("../phylogeny/215Fly_Orthogroups_tables_v5/HOG_stat_filtered.tsv")
dim(HOG_cons) # Should have nFlies + 1 columns

#### Adding tree id column to species info df ####
if (!file.exists(file = "data/Species215.RData")) {
  # species names and stats
  species <- read_xlsx("../phylogeny/215Fly_Orthogroups_tables_v5/Species_summary_215Fly.xlsx",
                       col_names = TRUE)
  
  # Issue: some spids in the species xlsx info sheet are shortenings or other
  # abbreviations of the species ID in the gene tree files,
  # which is a problem for treerecs
  # Treerecs will be ok if we have a very simple deliminated text file dictionary
  
  # But first we have to add a column to species that has all those
  # shortened tree ids in it
  is_conserved_single_copy_HOG <- apply(HOG_cons[,-1], 1, FUN = \(x) {
    all(x == 1)
  })
  HOG_cons$HOG[is_conserved_single_copy_HOG] |> sample(1)
  full_HOG <- "N1.HOG0013221"
  # reading names in from an alignment file
  full_tree_ids <- read_delim(delim = "\n", col_names = FALSE,
                              file = paste0("../phylogeny/215Fly_Orthogroups_tables_v5/HOG_Alignments/", 
                                            full_HOG, "/", full_HOG, "_final_align_AA.aln"),
                              col_types = "c")
  full_tree_ids <- full_tree_ids$X1[grepl(pattern = "^>", 
                                          x = full_tree_ids$X1)] |>
    str_split_i(pattern = "_", i = 1) |> 
    gsub(pattern = "^>", replacement = "")
  full_tree_ids[1:3]
  table(full_tree_ids) |> table() # 215 should all be 1
  # The Dmel tree_id is the transcript name in dmel (rna-NM_##.##### or something):
  dmel |> filter(HOG == full_HOG) |> select(Dmel_transcripts) |> 
    pull()
  # So splitting it by _ removes the transcript-unique part:
  dmel |> filter(HOG == full_HOG) |> select(Dmel_transcripts) |> 
    pull() |> str_split_i(pattern = "_", i = 1)
  grep(full_tree_ids, pattern = "rna-NM", value = TRUE)
  # adding tree_ids to species
  # first manually check that the names match up (removing dmel):
  spids_minusDmel <- setdiff(species$spid, "Dmelanogas")
  tree_ids_minusDmel <- setdiff(full_tree_ids, "rna-NM")
  spid_tree_lookup <- tibble(spid = sort(spids_minusDmel), 
                             tree_id = sort(tree_ids_minusDmel)) |> 
    bind_rows(tibble(spid = "Dmelanogas", 
                     tree_id = "rna-NM"))
  spid_tree_lookup |> head()
  spid_tree_lookup |> tail() # any mis-matches in the middle will result in the species at the end being mismatched
  # Annoyingly, the tree id isn't even always a prefix of spid, but sorting does seem to work
  
  # finally adding the tree_id column to species
  species <- left_join(species, spid_tree_lookup, by = "spid")
  # saving new version of species dataframe
  save(species, file = "data/Species215.RData")
}

# species data plus the above column
load(file = "data/Species215.RData")

#### Creating a dictionary of species names, gene tree names, and spids for treerecs

# Each HOG has different gene names for each species gene copies, 
# so we'll have to make a dictionary giving the gene name and corresponding
# species name for each gene in the HOG
writeTreeIDDictionary <- function(.HOG, 
                                  .folder = "../phylogeny/215Fly_Orthogroups_tables_v5/HOG_Dictionaries/") {
  # Getting all the tree_ids for this HOG,
  # With and without the _M####### identifier that is unique to each gene copy
  tree_names <- tryCatch({
    read_delim(delim = "\n", col_names = FALSE,
               file = paste0("../phylogeny/215Fly_Orthogroups_tables_v5/HOG_Alignments/", 
                             .HOG, "/", .HOG, "_final_align_AA.aln"),
               col_types = "c")
  }, error = function(e) {
    cat(.HOG, "is not in gene tree directory, returning NULL\n")
    return(NULL)
  })
  if (!is.null(tree_names)) {
    tree_names <- tree_names$X1[grepl(pattern = "^>", x = tree_names$X1)] |>
      gsub(pattern = "^>", replacement = "")
    tree_names_minusM <- str_split_i(tree_names, pattern = "_", i = 1)
    tree_lookup <- tibble(tree_name = tree_names,
                          tree_name_minusM = tree_names_minusM)
    speciesdf <- species |> select(species_name, tree_id) |> 
      filter(tree_id %in% tree_names_minusM) |> 
      mutate(species_name = toupper(species_name))
    dict <- left_join(tree_lookup, speciesdf, by = c("tree_name_minusM"="tree_id")) |> 
      select(tree_name, species_name)
    cat("writing dictionary to file\n")
    write_delim(dict, file = paste0(.folder, .HOG, ".txt"),
                delim = " ", col_names = FALSE)
    return(.HOG)
  }
}
# Tests for writeTreeIDDictionary
missing_HOG <- "N1.HOG0000425"
test_HOG <- "N1.HOG0000454" # Send1, spermatothecal enzyme
Diap1_HOG <- "N1.HOG0007077"
Diap2_HOG <- "N1.HOG0010308"
writeTreeIDDictionary(.HOG = missing_HOG)
# test <- writeTreeIDDictionary(.HOG = test_HOG)
# test <- writeTreeIDDictionary(.HOG = Diap1_HOG)
# test <- writeTreeIDDictionary(.HOG = Diap2_HOG)

### making all the dictionaries (only need to run this once)
# dmel_HOGs <- dmel |> drop_na(Dmel_FBgnID) |> select(HOG)
# write_delim(dmel_HOGs[1,], delim = "\n", col_names = FALSE,
#             file = "../phylogeny/treerecs_data/dmel_HOGs_test.txt")
# good_HOGs <- map(dmel_HOGs$HOG, writeTreeIDDictionary) |> unlist()
# write_delim(tibble(HOG = good_HOGs), delim = "\n", col_names = FALSE,
#             file = "../phylogeny/treerecs_data/dmel_HOGs.txt")

#### Testing treerecs parameters ####
# Making a HOGlist for testing treerecs contraction threshold and dup/loss costs
write_delim(tibble(HOG = c(test_HOG, Diap1_HOG, Diap2_HOG)), delim = "\n", col_names = FALSE,
                         file = "../phylogeny/treerecs_data/test_HOGs.txt")



#### Loading in the full Treerecs output ####
treerecs <- read_delim(file = "../phylogeny/treerecs_data/treerecs_output/combined_output.txt",
                       delim = "\n", skip = 2)
names(treerecs) <- "first_line"
treerecs$HOG <- str_split_i(treerecs$first_line, pattern = "_", i = 1)
treerecs$duplications <- str_extract(treerecs$first_line, "(?<=duplications =\\s)[0-9]+") |> 
  as.numeric()
treerecs$losses <- str_extract(treerecs$first_line, "(?<=losses =\\s)[0-9]+") |> 
  as.numeric()

#### Visualizing treerecs inferred duplications and losses ####
# Without annotation, b/c Dmel dupes are dropped with my current joining method
plotdf <- treerecs |> select(HOG, duplications, losses) |> 
  left_join(dmel, by = "HOG")
ggplot(plotdf, aes(x = duplications, y = losses)) +
  geom_point(alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_classic()

ggplot(plotdf, aes(x = duplications, y = losses)) +
  geom_point(alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_classic() +
  xlim(c(0, 200)) +
  ylim(c(0, 200))
ggplot(plotdf, aes(x = losses/duplications)) +
  geom_density() +
  theme_classic()
ggplot(plotdf, aes(x = duplications)) +
  geom_histogram(binwidth = 1) +
  theme_classic()
ggplot(plotdf, aes(x = losses)) +
  geom_histogram(binwidth = 1) +
  theme_classic()

# Why are there so many losses around 75?
loss_table <- treerecs$losses |> table()
loss_table[order(as.numeric(names(loss_table)))]
# Do the n losses correspond to the number of 0s in the phylogeny?
plotdf <- left_join(plotdf,
                    tibble(HOG = HOG_cons$HOG, 
                           n0 = apply(HOG_cons[,-1], 1, \(x) {sum(x == 0)})),
                    by = "HOG")
ggplot(plotdf, aes(x = losses, y = n0)) + geom_point()
# Not really? Although one loss can result in any number of 0s

# TODO: Are all 2000ish of those HOGs showing a similar loss pattern?
# save the problem HOGs into a new list and re-run the treerecs with svg output

plotdf <- tibble(mean_family_size = rowMeans(HOG_cons[HOG_cons$HOG %in% treerecs$HOG, -1]),
                 sum_family_size = rowSums(HOG_cons[HOG_cons$HOG %in% treerecs$HOG, -1]),
                 sd_family_size = apply(HOG_cons[HOG_cons$HOG %in% treerecs$HOG, -1], 1, FUN = sd),
                 HOG = HOG_cons$HOG[HOG_cons$HOG %in% treerecs$HOG]) |> 
  left_join(treerecs, by = "HOG")
ggplot(plotdf) +
  geom_histogram(aes(x = sum_family_size),
                 binwidth = 20) +
  geom_vline(xintercept = 215, color = "red") +
  theme_classic() +
  xlab("number of copies, summed for all 215 species")
plotdf$is_about_215 <- if_else(plotdf$sum_family_size >= 205 & 
                                 plotdf$sum_family_size <= 225, 
                               true = "Family size ~215 +/- 10",
                               false = "Family size rather different from 215")
ggplot(plotdf, aes(x = duplications)) +
  geom_histogram(binwidth = 1) +
  theme_classic() +
  facet_wrap(~is_about_215) +
  xlim(c(0, 500))
ggplot(plotdf, aes(x = losses)) +
  geom_histogram(binwidth = 1) +
  theme_classic() +
  facet_wrap(~is_about_215) +
  xlim(c(0, 500))

plotdf <- treerecs |> select(HOG, duplications, losses) |> 
  left_join(separate_longer_delim(select(dmel, HOG, Dmel_FBgnID, Dmel_CommonNames), 
                                  delim = ",", cols = c("Dmel_FBgnID", "Dmel_CommonNames")), 
            by = "HOG", relationship = "one-to-many") |> 
  left_join(williams, by = c("Dmel_FBgnID"="FlyBaseID")) |> 
  select(HOG, duplications, losses, Dmel_CommonNames, Plei_allImmune_allDev)
plotdf$Plei_allImmune_allDev[is.na(plotdf$Plei_allImmune_allDev)] <- "Other"

library(ggExtra)
p <- ggplot(filter(plotdf, Plei_allImmune_allDev != "Other"), 
            aes(x = duplications, y = losses)) +
  geom_point(aes(color = Plei_allImmune_allDev)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_classic() +
  theme(legend.position = "left")
ggMarginal(p, groupFill = TRUE, groupColour = TRUE, margins = "both")

ggplot(plotdf, aes(x = losses/duplications)) +
  geom_density(aes(fill = Plei_allImmune_allDev)) +
  facet_wrap(~Plei_allImmune_allDev)
ggplot(filter(plotdf, Plei_allImmune_allDev != "Other"),
       aes(x = duplications/losses)) +
  geom_density(aes(fill = Plei_allImmune_allDev), alpha = 0.5) +
  xlim(c(0, 1))
p <- ggplot(filter(plotdf, Plei_allImmune_allDev != "Other"),
       aes(x = duplications/losses, 
           y = factor(Plei_allImmune_allDev,
                      levels = c("Immune_Non_Pleiotropic",
                                 "Pleiotropic",
                                 "Developmental_Non_Pleiotropic")))) +
  geom_point(aes(color = Plei_allImmune_allDev), alpha = 0) +
  geom_boxplot(aes(fill = Plei_allImmune_allDev), alpha = 0.5) +
  xlim(c(0, 1)) +
  theme_classic() +
  ylab("")
ggMarginal(p, margins = "x", groupColor = TRUE, groupFill = TRUE) 

p <- ggplot(filter(plotdf, Plei_allImmune_allDev != "Other"),
            aes(x = duplications, 
                y = factor(Plei_allImmune_allDev,
                           levels = c("Immune_Non_Pleiotropic",
                                      "Pleiotropic",
                                      "Developmental_Non_Pleiotropic")))) +
  geom_point(aes(color = Plei_allImmune_allDev), alpha = 0) +
  geom_boxplot(aes(fill = Plei_allImmune_allDev), alpha = 0.5) +
  xlim(c(0, 1)) +
  theme_classic() +
  ylab("")
ggMarginal(p, margins = "x", groupColor = TRUE, groupFill = TRUE) 

plotdf |> filter(duplications > losses/4) |> View()

ggplot(plotdf, aes(y = losses, x = Plei_allImmune_allDev == "Immune_Non_Pleiotropic")) +
  geom_violin(aes(fill = Plei_allImmune_allDev)) +
  theme_classic() +
  stat_compare_means()
  
plotdf <- tibble(mean_family_size = rowMeans(HOG_cons[HOG_cons$HOG %in% treerecs$HOG, -1]),
                 sum_family_size = rowSums(HOG_cons[HOG_cons$HOG %in% treerecs$HOG, -1]),
                 sd_family_size = apply(HOG_cons[HOG_cons$HOG %in% treerecs$HOG, -1], 1, FUN = sd),
                 HOG = HOG_cons$HOG[HOG_cons$HOG %in% treerecs$HOG]) |> 
  left_join(dmel, by = "HOG") |> left_join(williams, by = c("Dmel_FBgnID"="FlyBaseID")) |> 
  select(HOG, sum_family_size, mean_family_size, sd_family_size, 
         Dmel_CommonNames, Plei_allImmune_allDev) |> 
  drop_na()

library(forcats)
ggplot(plotdf) +
  geom_histogram(aes(x = sum_family_size, fill = Plei_allImmune_allDev),
                 binwidth = 20) +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  facet_wrap(~Plei_allImmune_allDev, nrow = 3) +
  xlab("number of copies, summed for all 215 species")

ggplot(plotdf, aes(y = log(sd_family_size), x = Plei_allImmune_allDev)) +
  geom_violin(aes(fill = Plei_allImmune_allDev)) +
  theme_classic() + stat_compare_means()

-### Test: Checking random HOG alignments to verify the species the same as HOG_cons
# TODO: Once the alingment tables finish downloaded and I can delete all the HOGs
# I don't want to test, I can run some tests on selected HOGs:

# Test 1: N1.HOG0015892, randomly selected
# test_HOG <- "N1.HOG0015892"
# Test 2: N1.HOG0000454, some duplications present in in HOG_cons
# test_HOG <- "N1.HOG0000454"
# Test3: "N1.HOG0013221", one copy for all 301 species, to create a dictionary of spids

# testing
test <- read_delim(delim = "\n", col_names = FALSE,
                    file = paste0("../phylogeny/215Fly_Orthogroups_tables_v5/HOG_Alignments/", 
                                  test_HOG, "/", test_HOG, "_final_align_AA.aln"),
                    col_types = "c")
test <- test$X1[grepl(pattern = "^>", x = test$X1)] |> 
  str_split_i(pattern = "_", i = 1) |> 
  gsub(pattern = "^>", replacement = "")
test[1:3] # spid column in species
table(test) |> table() # shows duplications, if there are any

# For some reason, spids in the xlsx table can be shorter than their alignment name
test2 <- sapply(test, \(nm) {
  as.character(species$spid[startsWith(prefix = species$spid, x = nm)][1])
})
tibble(pre_shortening = test, 
       post_shortening = test2)
test <- test2
rm(test2)

# Verifying in HOG_cons
# should be zero
no_species <- species |> filter(!(spid %in% test)) |> 
  select(species_name) |> pull() |> toupper()
HOG_cons |> filter(HOG == test_HOG) |> select(all_of(no_species)) |> 
  as.numeric() |> table() # current issue: not 0
# TODO eventually: the above name repair is not perfect yet because the spid
# names aren't just prefixes of the alignment name---at least 8 of them
# are some other abbreviation. I don't currently think it's worth it to fix
# this, but it could be considered. For now, I trust that HOG_cons gives
# the correct copy number for each species

# should be length test
yes_species <- species$spid[startsWith(x = test, prefix = species$spid)] |> 
  select(species_name) |> pull() |> toupper()
HOG_cons |> filter(HOG == test_HOG) |> select(all_of(yes_species)) |> 
  sum()
length(test)

#### Quantifying copy number variation for each HOG ####
# Do pleiotropic genes have more copy number variation?
plotdf <- tibble(HOG = HOG_cons$HOG,
                 copy_number_sd = apply(HOG_cons[,-1], FUN = sd, MARGIN = 1),
                 copy_number_mean = apply(HOG_cons[,-1], FUN = mean, MARGIN = 1)) |> 
  left_join(select(dmel, -Dmel_transcripts, -Dmel_geneIDs), by = "HOG") |> 
  drop_na() |> 
  left_join(williams, by = c("Dmel_FBgnID"="FlyBaseID")) |> 
  drop_na(Plei_immuneResponse_embDev)
ggplot(filter(plotdf, copy_number_sd > 0 & copy_number_mean > 0),
       aes(y = log(copy_number_sd/copy_number_mean), 
           x = Plei_immuneResponse_embDev)) +
  geom_violin(aes(fill = Plei_immuneResponse_embDev)) +
  stat_compare_means() # immune do, but these include losses

# Are pleiotropic genes more likely to be retained in duplicate?
# TODO: to do this right, we'll need to infer duplication events
# Treerecs seems like a good tool for this. The trick will be getting
# The species tree in the right format

#### Finding duplicates that occurred outside of Dmel ####
SpeciesNames <- setdiff(colnames(HOG_cons), "HOG")

# this will likely be more useful for phylogenetic inference on where
# the duplication occurred
dupdf <- HOG_cons |> filter(DROSOPHILA_MELANOGASTER == 1) |> 
  filter(if_any(all_of(setdiff(SpeciesNames, 
                               "DROSOPHILA_MELANOGASTER")),
                ~ .x > 1)) |> 
  pivot_longer(cols = all_of(SpeciesNames),
               names_to = "species",
               values_to = "nCopies")

### Quantile-based criteria
# finding HOGs where Dmel copy number is below quantileThresh (but not 0)
quantileThresh <- 0.25
HOG_q_low <- apply(HOG_cons[, -1], FUN = quantile, MARGIN = 1, 
               probs = quantileThresh)
hist(HOG_q_low)
HOG_q_high <- apply(HOG_cons[, -1], FUN = quantile, MARGIN = 1, 
                   probs = 1-quantileThresh)
hist(HOG_q_high)

HOGdf <- tibble(HOG = HOG_cons$HOG,
                Dmel_High = HOG_q_high < HOG_cons$DROSOPHILA_MELANOGASTER,
                Dmel_Low = HOG_q_low > HOG_cons$DROSOPHILA_MELANOGASTER,
                Dmel_0 = HOG_cons$DROSOPHILA_MELANOGASTER == 0) |> 
  left_join(select(dmel, -Dmel_transcripts, -Dmel_geneIDs), by = "HOG")

HOGdf |> filter(!Dmel_0) |> select(Dmel_High, Dmel_Low) |> table()
# 17 isn't a large number to work with, discrete criteria will probably be better

### Discrete criteria
minSpecies <- 10
maxDmel <- 1
minDmel <- 1

# finding HOGs where Dmel copy number is smaller than the maximum seen
# in at least minSpecies, provided Dmel has no more than maxDmel

checkIfDmelHasLowCopyNumber <- function(.HOGs, .minSpecies,
                                        .minDmel, .maxDmel) {
  HOG_is_low_Dmel <- apply(HOG_cons[, -1], MARGIN = 1, FUN = \(x) {
    dmel <- x["DROSOPHILA_MELANOGASTER"]
    max_x <- max(x)
    n_max_x <- sum(x == max_x)
    isLowDmel <- (dmel <= .maxDmel) & 
      (dmel >= .minDmel) &
      (n_max_x >= .minSpecies) &
      (max_x > dmel)
    return(isLowDmel)
  })
  return(HOG_is_low_Dmel)
}

# checking a few values of minSpecies
nSpeciesdf <- map(c(2:20), .f = checkIfDmelHasLowCopyNumber,
              .HOGs = HOG_cons, .minDmel = minDmel, 
              .maxDmel = maxDmel) |> reduce(.f = bind_cols)
colnames(nSpeciesdf) <- paste0("nSpecies_", 2:20)
plot(x = c(2:20), y = colSums(nSpeciesdf),
     xlab = "Min nSpecies with higher HOG number",
     ylab = "Number of low copy HOGs in Dmel")
# Conclusion: probably worth having a non-conservative (2),
# med-conservative (5), and quite conservative (10) cutoff

HOGdf <- tibble(HOG = HOG_cons$HOG,
                low_copy_2 = nSpeciesdf$nSpecies_2,
                low_copy_5 = nSpeciesdf$nSpecies_5,
                low_copy_10 = nSpeciesdf$nSpecies_10,
                nDmel = HOG_cons$DROSOPHILA_MELANOGASTER) |> 
  filter(nDmel > 0) |> 
  left_join(select(dmel, -Dmel_transcripts, -Dmel_geneIDs), by = "HOG")

HOGdf |> select(nDmel, low_copy_5) |> table()

# Saving HOGdf
save(HOGdf, file = "data/HOG.RData")

#### Previous analysis that uses orthogroup presence and absence to assess duplication ####
# orthogroup info that I don't think we need anymore, b/c HOGs include all the same info
#HOG_to_OG <- read_tsv("../phylogeny/215Fly_Orthogroups_tables/HOG_to_OG_table.tsv")
#ortho <- read_tsv("../phylogeny/215Fly_Orthogroups_tables/Orthogroups.tsv")
#ortho_GeneCounts <- read_tsv("../phylogeny/215Fly_Orthogroups_tables/Orthogroups.GeneCount.tsv")

#### Validating/exploring orthogroup datasets ####
# duplications that aren't 
# found across all species I think are rows in ortho
# where there are NA values for certain species
# Why do I think that? 93% of NAs in Dmel or Dsim are NAs in both:
sum(is.na(ortho$DROSOPHILA_MELANOGASTER) & 
      is.na(ortho$DROSOPHILA_SIMULANS))/sum(is.na(ortho$DROSOPHILA_MELANOGASTER) |
                                              is.na(ortho$DROSOPHILA_SIMULANS))

### Polarizing duplication history by Dmel
# Which duplications occurred before Drosophilidae split?
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