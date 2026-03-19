sapply(c("tidyr", "dplyr", "ggplot2", "ggpubr", "ggalluvial", "purrr", "readxl", "readr", "stringr", "plotly", "lmtest", "pheatmap"), FUN = require,
       character.only = TRUE)

source("utils.r")

williams <- read_xlsx(path = "data/17Jan2023_Supplemental_File_1.xlsx",
                      sheet = 2, na = "NA") |> 
  select(FlyBaseID, Plei_allImmune_allDev) 
# using the least-strict definition of pleiotropy,
# because their strict threshold was based on 
# expression, which we can make our own judgments on

load("data/ImmuneCounts.RData")
counts_immune <- counts
clustdf_immune <- tibble("gene_name" = counts_immune$gene_name,
                         "environment" = "immune")
infodf_immune <- tibble("sample_name" = colnames(counts_immune)[-1],
                        "replicate" = if_else(grepl("A", colnames(counts_immune)[-1]),
                                                    true = "A", false = "B"),
                        "hour" = parse_number(colnames(counts_immune)[-1]))
rm(counts)

load("data/DevCounts.RData")
counts_dev <- counts
clustdf_dev <- tibble("gene_name" = counts_dev$gene_name,
                      "environment" = "dev")
infodf_dev <- infodf |> select(-sample_name) |> dplyr::rename("sample_name"="colname")
rm(counts, infodf)

load("data/OogCounts.RData")
counts_oog <- counts
clustdf_oog <- tibble("gene_name" = counts_oog$gene_name,
                      "environment" = "oog")
infodf_oog <- infodf
rm(counts, infodf)

#### Limiting Oogenesis to WholeOvary ####
samples_wholeOvary <- grep("^WholeOvary", colnames(counts_oog), value = TRUE)
counts_oog <- counts_oog[, c("gene_name", samples_wholeOvary)]
infodf_oog <- filter(infodf_oog, sample_name %in% samples_wholeOvary)

#### Combining datasets ####
# combining clustering results
clustdf <- bind_rows(clustdf_immune, clustdf_dev, clustdf_oog) |> 
  left_join(y = select(williams, FlyBaseID, 
                       Plei_allImmune_allDev),
            by = c("gene_name"="FlyBaseID"))
clustdf$williams_category <- if_else(is.na(clustdf$Plei_allImmune_allDev),
                                     true = "Other",
                                     false = clustdf$Plei_allImmune_allDev)

# combining counts and infodfs
counts <- list("immune" = counts_immune,
               "dev" = counts_dev,
               "oog" = counts_oog)
infodf <- list("immune" = infodf_immune,
               "dev" = infodf_dev,
               "oog" = infodf_oog)

#### Saving counts and info lists ####
save(counts, infodf, file = "data/CountsList.RData")

#### Parsing sample metadata ####
# getter function for any column in the infodf that is needed
# obtained by sample name
getInfo <- function(.sample_vec, .environment, .info_colname) {
  info_env <- infodf[[.environment]]
  mask <- sapply(.sample_vec, \(s) which(info_env$sample_name == s))
  out_vec <- info_env[mask, .info_colname] |> unlist() |> as.numeric()
  return(out_vec)
}
# test for getInfo
test_col <- colnames(counts_immune)[sample(2:ncol(counts_immune), size = 1)]
test_col
getInfo(test_col, .environment = "immune", .info_colname = "hour")
test_cols <- colnames(counts_oog)[sample(2:ncol(counts_oog), size = 5)]
test_cols
getInfo(test_cols, .environment = "oog", .info_colname = "hour")

#### Functions for Signal to Noise ####
# We'll define "signal" as the proportion of the expression variance (ESS,
# explained sum of squares) captured
# in a linear model with time as the predictor with a third degree polynomial fit
# The third degree polynomial should capture expression tending to increase
# or decrease over time, but not necessarily with one single peak or dip
# like it would with a second degree polynomial
# "noise" is therefore the residual sum of squares (RSS) variance
# Very old school simple linear regression: Total SS = ESS + RSS
# So our signal to noise ratio will be ESS/RSS. If we want it to be from 0
# to 1, we can change to R^2: ESS/TSS. But I think we don't really want this
# because we want the genes with more variability to still have that variability

### Helper for modeling, get the counts and time data for the model from
# the gene name and environment (for looping through clustdf)
getData <- function(.gene, .environment) {
  cts <- counts[[.environment]]
  counts_vec <- cts[cts$gene_name == .gene, -1] |> as.numeric()
  hour_vec <- getInfo(colnames(cts)[-1], 
                      .environment = .environment,
                      .info_colname = "hour")
  return(list(counts = counts_vec, hour = hour_vec))
}

### Fit our polynomial linear model
fitModel <- function(.counts_vec, .hour_vec, .degree = 3) {
  x <- data.frame(count = .counts_vec,
                  hour = .hour_vec)
  mod <- lm(count ~ poly(hour, degree = .degree), data = x)
  return(mod)
}
# tests for fitModel and getData
test_data <- getData(.gene = "FBgn0034329", .environment = "immune") # Bomanin Short 1
test_mod <- fitModel(.counts_vec = test_data$counts, 
                     .hour_vec = test_data$hour)
test_anova <- anova(test_mod)
test_anova[1,"Sum Sq"]/test_anova["Residuals", "Sum Sq"]

### fits a null model where time is not a predictor, for likelihood ratio test
fitNullModel <- function(.counts_vec, .hour_vec) {
  x <- data.frame(count = .counts_vec,
                  hour = .hour_vec)
  mod <- lm(count ~ 1, data = x)
  return(mod)
}
# tests for fitNullModel
test_data <- getData(.gene = "FBgn0034329", .environment = "immune") # Bomanin Short 1
test_null_mod <- fitNullModel(.counts_vec = test_data$counts, 
                              .hour_vec = test_data$hour)
test_full_lik <- logLik(test_mod)
test_null_lik <- logLik(test_null_mod)
test_lrt_stat <- -2 * (as.numeric(test_null_lik)-as.numeric(test_full_lik)) # why scale by -2? No idea but it makes it a chisq distributed stat: https://api.rpubs.com/tomanderson_34/lrt
pchisq(test_lrt_stat, df = 1, lower.tail = FALSE)

### get the anova info
getSNR <- function(.mod) {
  mod_anova <- anova(.mod)
  snr <- mod_anova[1,"Sum Sq"]/mod_anova["Residuals", "Sum Sq"]
  return(snr)
}
# tests for getSNR
test_anova <- anova(test_mod)
test_anova[1,"Sum Sq"]/test_anova["Residuals", "Sum Sq"]
getSNR(test_mod) # should be the same as above

### get model coefficients and model test results
getPolyCoefficients <- function(.mod, .null_mod) {
  coefs <- coefficients(.mod)
  coef_names <- c("Intercept", paste0("coef", 1:(length(coefs) - 1)))
  # likelihood ratio test of model fit
  full_lik <- logLik(.mod)
  null_lik <- logLik(.null_mod)
  lrt_stat <- -2 * (as.numeric(null_lik)-as.numeric(full_lik))
  lrt_pval <- pchisq(lrt_stat, df = 1, lower.tail = FALSE)
  # formatting data for output as one row of tibble
  outdf <- tibble(name = coef_names,
                  value = as.numeric(coefs)) |> 
    pivot_wider(names_from = "name", values_from = "value")
  outdf$lrt_pval <- lrt_pval
  return(outdf)
}
# tests for getPolyCoefficients
test_data <- getData(.gene = "FBgn0034329", .environment = "immune") # Bomanin Short 1
test_mod <- fitModel(.counts_vec = test_data$counts, 
                     .hour_vec = test_data$hour,
                     .degree = 3)
test_null_mod <- fitNullModel(.counts_vec = test_data$counts, 
                              .hour_vec = test_data$hour)
getPolyCoefficients(.mod = test_mod, .null_mod = test_null_mod) |> t()

#### Calculating Signal to Noise ####
# dataframe of modeling results for each gene in each environment
modeldf <- map2(clustdf$gene_name, clustdf$environment,
                \(g, e) {
                  cat(which(paste(g, e) == paste(clustdf$gene_name, clustdf$environment)),
                      "/", nrow(clustdf), "\n")
                  e_degree <- if_else(e == "oog",
                                      true = 2, false = 3) # Oog only has 3 time categories, polynomial can only be degree 2
                  g_data <- getData(.gene = g, .environment = e)
                  mod <- fitModel(.counts_vec = g_data$counts,
                                  .hour_vec = g_data$hour,
                                  .degree = e_degree)
                  null_mod <- fitNullModel(.counts_vec = g_data$counts,
                                           .hour_vec = g_data$hour)
                  snr <- getSNR(mod)
                  coefdf <- getPolyCoefficients(.mod = mod, .null_mod = null_mod)
                  outdf <- tibble("gene_name" = g, "environment" = e,
                                  "signal_to_noise" = snr) |> 
                    bind_cols(y = coefdf)
                  return(outdf)
                }) |> purrr::reduce(.f = bind_rows)

# checking for missing values
sum(is.na(modeldf$signal_to_noise))
# adding modeling results to clustdf ()
clustdf <- left_join(clustdf, modeldf, by = c("gene_name", "environment"))

sort(clustdf$signal_to_noise, decreasing = TRUE)[1:10] # some major outliers with huge SNR
sort(clustdf$signal_to_noise, decreasing = FALSE)[1:20] # some major outliers, all in "Other" category
plotdf <- filter(clustdf, signal_to_noise < 1e30)
plotdf$williams_category <- factor(plotdf$williams_category,
                                   levels = c("Developmental_Non_Pleiotropic", "Immune_Non_Pleiotropic",
                                              "Pleiotropic", "Other"),
                                   labels = c("Developmental", "Immune",
                                              "Pleiotropic", "Neither"))
group_colors <- c("Developmental"="#E41A1CFF",
                  "Immune"="#377EB8FF",
                  "Pleiotropic"="#984EA3FF",
                  "Neither"="#4DAF4AFF")
p <- ggplot(plotdf, aes(x = williams_category, y = log2(signal_to_noise + 1e-10))) +
  #geom_jitter(aes(color = williams_category)) +
  geom_violin(aes(fill = williams_category, group = williams_category)) +
  stat_summary(fun = "mean", geom = "point") +
  geom_hline(data = summarise(group_by(filter(plotdf, williams_category == "Neither"),
                                       environment), meanSNR = mean(log2(signal_to_noise + 1e-10))),
             aes(yintercept = meanSNR)) +
  stat_compare_means(comparisons = list(c("Developmental", "Neither"),
                                        c("Immune", "Neither"),
                                        c("Pleiotropic", "Neither"))) +
  #scale_fill_brewer(type = "qual", palette = "Set1") +
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  facet_wrap(~factor(environment, levels = c("dev", "immune", "oog"), 
                     labels = c("Embryonic Development", "Imd Challenge", "Oogenesis"))) +
  xlab("gene category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none") +
  ylim(c(-10, log2(max(plotdf$signal_to_noise)) + 14))
p
# pdf("figures/violin.pdf", width = 5, height = 3.5)
# p
# dev.off()

#### Inspecting individual genes ####
# looking at our massive outliers
# high outliers
test_highSNR_df <- clustdf |> arrange(desc(signal_to_noise)) |> 
  slice_head(n = 10) |> 
  select(all_of(c("gene_name", "environment")))
test_highSNR_df$environment |> table() # oog
plotExpressionProfile(.counts = counts$oog, .info = infodf$oog,
                      .gene_idxs = test_highSNR_df$gene_name,
                      .gene_names = c(1:10))
# These look like genes in two categories:
# 1) genes with high variability over time (what we want)
# Or 2) genes with little to no variability over time,
#       but absolutely no variation between replicates (what we don't want)

# genes 1, 2, 4, and 7 are the main culprits of the second category
test_highSNR_df$gene_name[c(1,2,4,7)] # "FBgn0034808" "FBgn0003254" "FBgn0032213" "FBgn0053056"
counts$oog[counts$oog$gene_name == test_highSNR_df$gene_name[1],]
counts$oog[counts$oog$gene_name == test_highSNR_df$gene_name[2],]
counts$oog[counts$oog$gene_name == test_highSNR_df$gene_name[4],]
counts$oog[counts$oog$gene_name == test_highSNR_df$gene_name[7],]

# low outliers
test_lowSNR_df <- clustdf |> arrange(signal_to_noise) |> 
  slice_head(n = 10) |> 
  select(all_of(c("gene_name", "environment")))
test_lowSNR_df$environment |> table() # oog
plotExpressionProfile(.counts = counts$oog, .info = infodf$oog,
                      .gene_idxs = test_lowSNR_df$gene_name,
                      .gene_names = c(1:10)) # checks out. these are good ones to have near-0 SNR

# 10 random genes to compare replicate variation
test_idxs <- sample(counts$oog$gene_name, size = 10)
plotExpressionProfile(.counts = counts$oog, .info = infodf$oog,
                      .gene_idxs = test_idxs,
                      .gene_names = c(1:10))
counts$oog[counts$oog$gene_name == test_idxs[9],] # pick any one that looks like a flat line
# these do tend to have more var between replicates

# conclusions: SNR appears to be quantifying what we're interested in. But it
# is worth incorporating info on mean expression of each gene b/c the low expression
# genes that just happen to have low variation between replicates can have a very high SNR

#### Calculating mean expression ####
clustdf$mean <- map2(clustdf$gene_name, clustdf$environment,
                     \(g, e) {
                       g_data <- getData(.gene = g, .environment = e)
                       return(mean(g_data$counts))
                     }) |> unlist() # takes a moment

# makes SNR be between 0 and 1
normalizeSNR <- function(.snr_vec) {
  min_snr <- min(.snr_vec)
  out_vec <- .snr_vec + (0 - min_snr)
  new_max_snr <- max(out_vec)
  out_vec <- out_vec/new_max_snr
  return(out_vec)
}
# tests for normalizeSNR
test_devdf <- clustdf |> filter(mean >= 10 & 
                                  signal_to_noise < 100 &
                                  environment == "dev")
test_normSNR <- normalizeSNR(test_devdf$signal_to_noise)
quantile(test_normSNR, probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 1))
test_immunedf <- clustdf |> filter(mean >= 10 & 
                                     signal_to_noise < 100 &
                                     environment == "immune")
test_normSNR <- normalizeSNR(test_immunedf$signal_to_noise)
quantile(test_normSNR, probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 1))
# TODO: if we do this, probably not good to do a linear scale, everything gets piled
# up at the 99th percentile
# (currently we're doing dense_rank instead of normalizing)

#### creating data frame where each gene has 1 row ####
genedf <- clustdf |> 
  select(gene_name, environment, signal_to_noise, williams_category, mean) |> 
  pivot_wider(id_cols = c("gene_name", "williams_category"),
              names_from = "environment",
              values_from = c("signal_to_noise", "mean")) |> 
  # if a gene is missing/below expression threshold for a dataset, 
  # we set it to mean = 0 and SNR = 0 instead of NA
  mutate(signal_to_noise_immune = coalesce(signal_to_noise_immune, 0), 
         signal_to_noise_dev = coalesce(signal_to_noise_dev, 0),
         signal_to_noise_oog = coalesce(signal_to_noise_oog, 0),
         mean_immune = coalesce(mean_immune, 0),
         mean_dev = coalesce(mean_dev, 0),
         mean_oog = coalesce(mean_oog, 0))
# ranking SNR so it is comparable between datasets
nGenes <- nrow(genedf)
genedf$immune_rank <- dplyr::dense_rank(desc(genedf$signal_to_noise_immune)) # highest SNR = 1st rank
genedf$dev_rank <- dplyr::dense_rank(desc(genedf$signal_to_noise_dev))
genedf$oog_rank <- dplyr::dense_rank(desc(genedf$signal_to_noise_oog))
# rank threshold with uniform SNR cutoff
SNR_thresh <- 2
immune_rank_thresh <- max(genedf$immune_rank[which(genedf$signal_to_noise_immune > SNR_thresh)])
dev_rank_thresh <- max(genedf$dev_rank[which(genedf$signal_to_noise_dev > SNR_thresh)])
oog_rank_thresh <- max(genedf$oog_rank[which(genedf$signal_to_noise_oog > SNR_thresh)])
# # rank threshold with uniform nGenes
# immune_rank_thresh <- nGenes*0.25 # top 25% of genes
# dev_rank_thresh <- nGenes*0.25
# oog_rank_thresh <- nGenes*0.25

# assigning gene its variability category by SNR rank
genedf$max_snr <- map(c(1:nrow(genedf)), \(i) {
  rank_imm <- genedf$immune_rank[i]
  rank_dev <- genedf$dev_rank[i]
  rank_oog <- genedf$oog_rank[i]
  if (rank_imm > immune_rank_thresh) {
    rank_imm <- nGenes
  }
  if (rank_dev > dev_rank_thresh) {
    rank_dev <- nGenes
  }
  if (rank_oog > oog_rank_thresh) {
    rank_oog <- nGenes
  }
  if (rank_imm > immune_rank_thresh &
      rank_dev > dev_rank_thresh & 
      rank_oog > oog_rank_thresh) {
    return("low_SNR")
  }
  rank_vec <- c(rank_imm,
                rank_dev,
                rank_oog)
  return(c("immune", "dev", "oog")[which.min(rank_vec)])
}) |> unlist()

# does the gene have variable expression in the conditions(s) where it plays a role?
# Developmental_Non_Pleiotropic: Embryonic Development or Oogenesis
# Immune_Non_Pleiotropic: Imd Challenge
# PLeiotropic: (Embryonic Development or Oogenesis) AND Imd Challenge

#### saving ####
save(clustdf, genedf, SNR_thresh, file = "data/SignalToNoise.RData")

###################### Exploration ###################### 
#### 3D scatterplot of SNR ranks in each dataset ####
plotdf <- clustdf |> # filter(mean >= 10 & signal_to_noise < 100) |> 
  filter(williams_category != "Other") |>
  select(all_of(c("gene_name", "environment", "williams_category", "signal_to_noise"))) |> 
  pivot_wider(id_cols = c("gene_name", "williams_category"),
              names_from = "environment", values_from = "signal_to_noise") |> 
  drop_na() |> # removes genes not expressed in all three environments
  mutate(dev_norm = rank(dev, ties.method = "first"),
         immune_norm = rank(immune, ties.method = "first"),
         oog_norm = rank(oog, ties.method = "first"))
  # mutate(dev_norm = as.numeric(scale(dev)),
  #        immune_norm = as.numeric(scale(immune)),
  #        oog_norm = as.numeric(scale(oog)))


# Did scaling produce similar distributions of SNR?
ggplot(pivot_longer(plotdf, cols = c("immune_norm", "dev_norm", "oog_norm")),
       aes(x = value, fill = name)) + 
  geom_density(aes(fill = name), alpha = 0.5) # dev is the weird one

# 3D plot
group_colors <- c("Developmental_Non_Pleiotropic"="#E41A1CFF",
                  "Immune_Non_Pleiotropic"="#377EB8FF",
                  "Pleiotropic"="#984EA3FF",
                  "Other"="#4DAF4AFF")
exampledf_dev <- filter(plotdf, gene_name == "FBgn0003612") # Pleiotropic gene that is high variability in dev and low in both other environments
exampledf_immune <- filter(plotdf, gene_name == "FBgn0001125") # Dev gene that is high variability in immune and low in both other environments
exampledf_immune2 <- filter(plotdf, gene_name == "FBgn0029990") # Immune gene that is high variability in immune and low in both other environments
exampledf_immune3 <- filter(plotdf, gene_name == "FBgn0061356") # Immune gene that is high variability in immune and low in both other environments

fig <- plot_ly(plotdf, x = ~dev_norm, y = ~immune_norm, 
               z = ~oog_norm, color = ~williams_category,
               hoverinfo = 'text',
               text = ~gene_name,
               colors = group_colors, alpha = 0.5, size = 0.25)
fig <- fig |> add_markers()
fig <- fig |> layout(scene = list(xaxis = list(title = 'Dev'),
                                   yaxis = list(title = 'Immune'),
                                   zaxis = list(title = 'Oog')))
fig

### following up on examples from figure
# Plei genes wiht high var in Oog and Dev, low in Immune:
gene_idxs <- "FBgn0010389" # heartless
gene_idxs <- "FBgn0033438" # MMP2
# genes with high var in all three:
gene_idxs <- "FBgn0028540" # immune, not well characteriezd Predicted to enable glucose-6-phosphate 1-epimerase activity. Involved in defense response to virus. Predicted to be active in cytoplasm
gene_idxs <- "FBgn0261547" # plei, Ephexin a rhoGEF
gene_idxs <- "FBgn0027111" # dev, mple1

# plotting
p_immune <- plotExpressionProfile(.counts = counts$immune, 
                      .info = infodf$immune, 
                      .gene_idxs = gene_idxs,
                      .gene_names = gene_idxs)
p_dev <- plotExpressionProfile(.counts = counts$dev, 
                      .info = infodf$dev, 
                      .gene_idxs = gene_idxs,
                      .gene_names = gene_idxs)
p_oog <- plotExpressionProfile(.counts = counts$oog, 
                      .info = infodf$oog, 
                      .gene_idxs = gene_idxs,
                      .gene_names = gene_idxs)
ggarrange(p_immune, p_dev, p_oog, nrow = 1, ncol = 3,
          common.legend = TRUE)

#### Determining a signal-to-noise threshold for expression variability ####

ggplot(filter(clustdf, signal_to_noise < 1e30 & signal_to_noise > 1e-30), 
       aes(x = log2(signal_to_noise), y = -log10(lrt_pval))) +
  geom_point(aes(color = environment)) +
  geom_vline(xintercept = 1) + # (log2(2) = 1)
  facet_wrap(~environment)

# By an SNR of 2, all the pvalues clearly "lift off" of the x-axis, so we can
# use that as our cutoff for variability. Intuitively it makes sense, 
# twice as much signal than noise

# There has to be a mathematical reason why the likelihood ratio results and the 
# sum of squares ratio are monotonically related, they're essentially both ratios
# of complex vs simple models, but slightly different ones

#### Which dataset does each gene have its more variable expression in? ####
# max_snr assigns each gene to the dataset where it has the highest SNR 
# (by rank among all genes, rank of 1 is the gene with the highest SNR in that dataset)
plotdf <- group_by(genedf, max_snr, williams_category) |> 
  summarise(n = n()) |> ungroup()
plotdf_counts <- group_by(genedf, williams_category) |> 
  summarise(n = n()) |> ungroup()
plotdf$williams_category <- factor(plotdf$williams_category, 
                                   levels = c("Developmental_Non_Pleiotropic",
                                              "Immune_Non_Pleiotropic",
                                              "Pleiotropic",
                                              "Other"))
plotdf$max_snr <- factor(plotdf$max_snr, 
                              levels = c("dev", "immune", "oog", "low_SNR"),
                              labels = c("Embryonic Development",
                                         "Imd Challenge",
                                         "Oogenesis",
                                         "None (lowly-varying)"))
# stacked bar proportions
ggplot(plotdf, aes(x = williams_category, y = n)) +
  geom_bar(aes(fill = max_snr), stat = "identity", position = "fill") +
  geom_text(data = plotdf_counts, aes(label = paste0("n = ", n),
                                      x = williams_category, y = 1.1)) +
  xlab("function of gene") +
  scale_fill_brewer(name = "environment of most\nvariable expression",
                    palette = "Set1") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# But some genes may have plenty of variability in multiple datasets, so we should
# also check a Venn diagram visualization
library(ggVennDiagram)
# this hurts my head every time
# instead of tidy data, we have a named list where each name is a dataset
# that each gene can or cannot have variable expression in
# the items in the list are genes that do have variable expression in that dataset
plotVenn <- function(.df, .williams_category) {
  plotlist <- list(Embryonic_Development = .df |> 
                         filter(williams_category == .williams_category) |> 
                         filter(dev_rank < dev_rank_thresh) |> 
                         select(gene_name) |> pull(),
                       Imd_Challenge = .df |> 
                         filter(williams_category == .williams_category) |> 
                         filter(immune_rank < immune_rank_thresh) |> 
                         select(gene_name) |> pull(),
                       Oogenesis = .df |> 
                         filter(williams_category == .williams_category) |> 
                         filter(oog_rank < oog_rank_thresh) |> 
                         select(gene_name) |> pull())
  p <- ggVennDiagram(plotlist) + scale_fill_gradient(low="grey90",high = "red")
  return(p)
}

plotVenn(genedf, .williams_category = "Developmental_Non_Pleiotropic") +
  ggtitle("Developmental, non-pleiotropic genes") # 4% immune
plotVenn(genedf, .williams_category = "Immune_Non_Pleiotropic") +
  ggtitle("Immune, non-pleiotropic genes") # 16% immune
plotVenn(genedf, .williams_category = "Pleiotropic") +
  ggtitle("Pleiotropic genes")# 8% immune

# Conclusion: There are differences in proportion for each class of gene, but
# modest ones. Especially for the developmental dataset, it is important to 
# distinguish between "active" and "passive" expression shifts. "Active" shifts
# are created by the genes that need to change expression because they
# have a role in the process. "Passive" are genes that change expression as a
# result of the "active" group changing. These genes theoretically would not need
# to change expression if there was infinite transcriptional machinery in the cell

# The way we will distinguish these classes is based on Williams et al. functional
# annotation. "Active" genes are the genes williams et al. annotates as having
# a role in that process (Development, Immune, or Both). "Passive" are the genes that
# are also changing expression but don't have a role in that environment
# (i.e. Immune_Non_Pleiotropic genes that
# have variable expression in Embryonic Development)

#### Clustering variable genes based on polynomial coefficients ####
# plot polynomial curve
fitCurve <- function(.coefs, .x_values) {
  nterms <- length(.coefs) - 1
  y <- map2(.x = c(0:nterms),
            .y = .coefs, \(term, cf) {
              return(cf*.x_values^term)
            }) |> Reduce(f = `+`)
  #y <- .coefs[1] + .coefs[2]*x + .coefs[3]*x^2 + .coefs[4]*x^3
  return(y)
}
test_x <- runif(min = -2, max = 2, n = 1000)
test_y <- fitCurve(.coefs = as.numeric(c(0, -2, 50, -40)),
                   .x_values = test_x)
plot(test_x, test_y)
abline(v = 0, col = "red")
abline(h = 0, col = "red")

### Heatmap of polynomial coefficients in each environment

# signal-to-noise <= 1
plotdf <- clustdf |> filter(signal_to_noise <= 1) |> 
  select(all_of(c("gene_name", "environment", "coef1", "coef2", "coef3"))) |> 
  pivot_wider(id_cols = "gene_name", names_from = "environment", 
              values_from = c("coef1", "coef2", "coef3")) |> 
  select(-all_of("coef3_oog")) |> 
  drop_na()
plot_mat <- as.matrix(plotdf[,-1])
rownames(plot_mat) <- plotdf$gene_name
plot_mat <- plot_mat[, c("coef1_immune", "coef2_immune", "coef3_immune",
                         "coef1_dev", "coef2_dev", "coef3_dev",
                         "coef1_oog", "coef2_oog")]
plot_mat <- apply(plot_mat, 1, scale)
pheatmap(t(plot_mat), cluster_cols = FALSE, show_rownames = FALSE)

# signal-to-noise > 1, immune
plotdf <- clustdf |> filter(environment == "immune" & signal_to_noise > 1) |> 
  select(all_of(c("gene_name", "coef1", "coef2", "coef3"))) |>
  drop_na()
plot_mat <- as.matrix(plotdf[,-1])
rownames(plot_mat) <- plotdf$gene_name
plot_mat <- apply(plot_mat, 1, scale, scale = TRUE, center = FALSE) |> t()
pheatmap(plot_mat, cluster_cols = FALSE, 
         show_rownames = FALSE,
         show_colnames = TRUE, na_col = "grey30",
         main = paste("Imd challenge,", nrow(plot_mat), "genes"))

# signal-to-noise > 1, dev
plotdf <- clustdf |> filter(environment == "dev" & signal_to_noise > 1) |> 
  select(all_of(c("gene_name", "coef1", "coef2", "coef3"))) |>
  drop_na()
plot_mat <- as.matrix(plotdf[,-1])
rownames(plot_mat) <- plotdf$gene_name
plot_mat <- apply(plot_mat, 1, scale, scale = TRUE, center = FALSE) |> t()
pheatmap(plot_mat, cluster_cols = FALSE, 
         show_rownames = FALSE,
         show_colnames = FALSE, na_col = "grey30",
         main = paste("Embryonic Development,", nrow(plot_mat), "genes"))

# signal-to-noise > 1, oog
plotdf <- clustdf |> filter(environment == "oog" & signal_to_noise > 1 &
                              signal_to_noise < 1e30) |> 
  select(all_of(c("gene_name", "coef1", "coef2"))) |>
  drop_na()
plot_mat <- as.matrix(plotdf[,-1])
rownames(plot_mat) <- plotdf$gene_name
plot_mat <- apply(plot_mat, 1, scale, scale = TRUE, center = FALSE) |> t()
pheatmap(plot_mat, cluster_cols = FALSE, 
         show_rownames = FALSE,
         show_colnames = TRUE, na_col = "grey30",
         main = paste("Oogenesis,", nrow(plot_mat), "genes"))

### Plot example genes from each cluster
# immune
test_gene_idxs <- clustdf |> 
  filter(environment == "immune" & signal_to_noise > 1.1 &
           coef1 > 0 & coef2 < -10 & coef3 > 10) |> 
  select(gene_name) |> pull()
test_gene_idxs # Syndecan
plotExpressionProfile(.counts = counts$immune, 
                      .info = infodf$immune, 
                      .gene_idxs = test_gene_idxs,
                      .gene_names = c(1:length(test_gene_idxs)))

# dev
test_gene_idxs <- clustdf |> 
  filter(environment == "dev" & signal_to_noise > 2 &
           coef1 < -100 & coef2 > 3000 & coef3 < -100) |> 
  select(gene_name) |> pull()
plotExpressionProfile(.counts = counts$dev, 
                      .info = infodf$dev, 
                      .gene_idxs = test_gene_idxs,
                      .gene_names = c(1:length(test_gene_idxs)))


# oog
test_gene_idxs <- clustdf |> 
  filter(environment == "oog" & signal_to_noise > 20 &
           coef1 < -100 & coef2 > 100) |> 
  select(gene_name) |> pull()
plotExpressionProfile(.counts = counts$oog, 
                      .info = infodf$oog, 
                      .gene_idxs = test_gene_idxs,
                      .gene_names = c(1:length(test_gene_idxs)))

# TODO maybe: This def does not do the right polynomials yet
plotExpressionProfileWithPolynomial <- function(.counts, .info, .gene_idxs, .gene_names, .polydf) {
  gene_name_lookup <- tibble(gene_idx = .gene_idxs,
                             gene = .gene_names)
  countdf <- map2(.gene_idxs, .gene_names, \(g, nm) {
    outdf <- tibble(sample = colnames(.counts)[-1],
                    gene_vec = as.numeric(.counts[.counts$gene_name == g, -1]))
    names(outdf) <- c("sample_name", nm)
    return(outdf)
  }) |> purrr::reduce(.f = left_join, by = "sample_name")
  plotdf <- countdf |> pivot_longer(cols = setdiff(names(countdf), "sample_name"),
                                    names_to = "gene", values_to = "expr") |> 
    left_join(y = .info, by = "sample_name") |> 
    left_join(y = gene_name_lookup, by = "gene")
  plotdf$model <- map2(plotdf$hour, plotdf$gene_idx, \(h, g) {
    coefs <- .polydf |> filter(gene_name == g) |>
      select(-all_of(c("gene_name", "environment", "signal_to_noise"))) |> 
      unlist() |> 
      as.numeric()
    return(fitCurve(.coefs = coefs, .x_values = h))
  }) |> unlist()
  ggplot(plotdf, aes(x = factor(hour), y = expr)) +
    geom_point(aes(color = factor(gene)),
               size = 0.25) +
    geom_line(aes(color = factor(gene),
                  group = interaction(gene, replicate),
                  linetype = factor(replicate))) +
    geom_line(aes(x = factor(hour), y = model, group = factor(gene))) +
    #scale_color_brewer(palette = "Set1") +
    labs(color = "Cluster") +
    theme_classic() +
    ylab("Expression (tpm)") +
    xlab("Time (hours)") +
    facet_wrap(~factor(gene, levels = .gene_names, 
                       labels = .gene_names))
}
# tests for plotExpressionProfileWithPolynomial
test_gene_idxs <- clustdf |> 
  filter(environment == "immune" & signal_to_noise > 1 &
           coef1 > 100 & coef2 < -10 & coef3 > 0) |> 
  select(gene_name) |> pull()
test_gene_idxs # Hrb27C, CG34166
plotExpressionProfileWithPolynomial(.counts = counts$immune, 
                                    .info = infodf$immune, 
                                    .gene_idxs = test_gene_idxs,
                                    .gene_names = c("Hrb27C", "CG34166"),
                                    .polydf = filter(modeldf, environment == "immune"))


################ Archive ################ 
#### Ternary Plot ####
# Archived b/c I didn't realize that by definition in 
# Ternary plots, you can't be high in 2/3 variables --- 
# If you are high in one variably you must be low in the other
# Two b/c their values must sum to 100% (i.e. % clay, sand, organic for soil)
# library(ggtern)
# ggtern(data = plotdf,
#        aes(x = immune_norm, y = dev_norm, z = oog_norm)) + # L = x, T = y, R = z
#   geom_point(aes(color = williams_category), alpha = 1, size = 0.5) +
#   geom_point(data = exampledf_dev,
#              color = group_colors[exampledf_dev$williams_category],
#              alpha = 1, size = 3) +
#   geom_point(data = exampledf_immune,
#              color = group_colors[exampledf_immune$williams_category],
#              alpha = 1, size = 3) +
#   geom_point(data = exampledf_immune2,
#              color = group_colors[exampledf_immune2$williams_category],
#              alpha = 1, size = 3) +
#   geom_point(data = exampledf_immune3,
#              color = "gold",
#              alpha = 1, size = 3) +
#   # facet_wrap(~factor(williams_category,
#   #                    levels = c("Developmental_Non_Pleiotropic",
#   #                               "Immune_Non_Pleiotropic",
#   #                               "Pleiotropic",
#   #                               "Other"),
#   #                    labels = c("Developmental",
#   #                               "Immune",
#   #                               "Pleiotropic",
#   #                               "Neither"))) +
#   scale_color_manual(values = group_colors) +
#   geom_Tline(Tintercept = 0.9) +
#   geom_Lline(Lintercept = 0.9) +
#   scale_L_continuous(limits = c(min(plotdf$immune_norm), max(plotdf$immune_norm))) +
#   scale_T_continuous(limits = c(min(plotdf$dev_norm), max(plotdf$dev_norm))) +
#   scale_R_continuous(limits = c(min(plotdf$oog_norm), max(plotdf$oog_norm))) +
#   tern_limit(T = 1, L = 1, R = 1)


