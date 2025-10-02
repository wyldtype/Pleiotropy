sapply(c("tidyr", "dplyr", "ggplot2", "purrr", "readxl", "readr", "stringr"), FUN = require,
       character.only = TRUE)
williams <- read_xlsx(path = "data/17Jan2023_Supplemental_File_1.xlsx",
                      sheet = 2, na = "NA") 
immune_lookup <- read_tsv(file = "data/immune_FlyBaseLookup.txt") 
colnames(immune_lookup) <- c("FlyBaseID", "gene_name", "symbol")
# any GO annotation in immune or dev
williams |> 
  select(Plei_allImmune_allDev) |> table(useNA = "always")
# responsive in immune, active all dev
williams |> 
  select(Plei_immuneResponse_allDev) |> table(useNA = "always")
# responsive in immune, active in embryonic dev
williams |> 
  select(Plei_immuneResponse_embDev) |> table(useNA = "always")

# separating pleiotropic and non-pleiotropic parts
pleidf <- williams |> select(FlyBaseID, Plei_immuneResponse_embDev) |> 
  filter(Plei_immuneResponse_embDev == "Pleiotropic")
immunedf <- williams |> select(FlyBaseID, Plei_immuneResponse_embDev) |> 
  filter(Plei_immuneResponse_embDev == "Immune_Non_Pleiotropic")
immunedf <- left_join(immunedf, immune_lookup,
                      by = "FlyBaseID")

# exporting for lit review
write_delim(pleidf, "data/pleiotropic_BLANK.csv", delim = ",")
write_delim(immunedf, "data/immune_BLANK.csv", delim = ",")

# importing for lit review
pleidf <- read_delim("data/pleiotropic.csv", delim = ",",
                     col_select = c(1:4))
colnames(pleidf) <- c("FlyBaseID", "gene_name", "pleioType", "mechanism")

# plotting
plotdf <- select(pleidf, gene_name, pleioType)
plotdf$pleioType <- gsub("\\?", "", plotdf$pleioType) # removing ?'s and /'s from pleioType
plotdf$pleioType1 <- sapply(plotdf$pleioType, \(x) {
  types <- str_split(x, pattern = "/") |> unlist()
  types[1]
})
plotdf$pleioType2 <- sapply(plotdf$pleioType, \(x) {
  types <- str_split(x, pattern = "/") |> unlist()
  types[2]
})
plotdf <- plotdf |> select(gene_name, pleioType1, pleioType2) |> 
  pivot_longer(cols = c("pleioType1", "pleioType2"),
               values_to = "pleioType") |> 
  drop_na()
plotdf$pleioType <- factor(plotdf$pleioType,
                           levels = c("not pleiotropic",
                                      "artefactual",
                                      "secondary",
                                      "adoptive",
                                      "parsimonious",
                                      "opportunistic",
                                      "combinatorial",
                                      "unifying"))
ggplot(plotdf, aes(x = pleioType)) +
  geom_bar(aes(fill = pleioType)) +
  theme_classic() +
  scale_x_discrete(drop = FALSE) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   vjust = 1),
        legend.position = "none")
  
