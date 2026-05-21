#### Williams et al. pleiotropic gene catalog ####
# pleiotropic gene classifications
williams <- read_xlsx(path = "data/17Jan2023_Supplemental_File_1.xlsx",
                      sheet = 2, na = "NA")
# Adding "most conservative option" which was not looked at in Williams et al. for some reason: 
# the conservative definition of development (embryonic development)
# plus the conservative definition of immune (immune functional annotation, not immune response)
williams$Plei_allImmune_embDev <- if_else(williams$Plei_allImmune_allDev == "Immune_Non_Pleiotropic",
                                          true = "Immune_Non_Pleiotropic",
                                          false = if_else(williams$GO_Dev_Emb == "embryonic_dev",
                                                          true = williams$Plei_immuneResponse_embDev,
                                                          false = NA))

#### Colors used for plotting ####
ColorValues <- c("Developmental_Non_Pleiotropic"="#F8766D",
                 "Pleiotropic"="#619CFF",
                 "Immune_Non_Pleiotropic"="#00BA38",
                 "Other"="grey80")

#### Drosophila gene summaries ####
# flybase gene summaries
drosophila_lookup <- read_delim("data/best_gene_summary_fb_2026_01.tsv", delim = "\t",
                                comment = "##")
colnames(drosophila_lookup) <- c("FlyBaseID", "Gene_Symbol", "Summary_Source", "Summary")
