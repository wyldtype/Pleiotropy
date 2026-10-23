sapply(c("tidyr", "dplyr", "ggplot2", "ggpubr", "ggalluvial", "purrr", "readxl", "readr", "stringr", "plotly", "lmtest", "pheatmap"), FUN = require,
       character.only = TRUE)

# Cohen et al. 2026, Imd enhancers
# loading enhancer assignments
enhancers <- read_csv("data/EG_TFBS_cat_01_19_activityclass.csv") |> 
  filter(Gene != ".") |>  # . are enhancers not assigned to any gene
# If we need to parse the Enhancer Name:
  mutate(Enhancer_Start = str_extract(Enhancer, "(?<=:)\\d+(?=-)"),
         Enhancer_End = str_extract(Enhancer, "(?<=-)\\d+"),
         Chromosome = str_extract(Enhancer, "^.*(?=:)")) |> 
  select(Enhancer, Enhancer_Start, Enhancer_End, Chromosome,
         Group, Gene, Rel_3, Accessibility, `Broad Immune Role`) |> 
  rename("ImmuneRole"="Broad Immune Role")

# which condition or combination of conditions the enhancer was active in:
table(enhancers$Group)

# number of enhancers with at least on relish motif
sum(enhancers$Rel_3 > 0)

# Version of enhancer dataset where there are separate
# rows for each gene assigned to an enhancer in a specific
# activity group
enhancers_longer <- enhancers |>
  separate_longer_delim(col = "Gene", delim = ";") 
# Note: even though it would be great to also be able to separate_longer_delim the ImmuneRole column, it doesn't work b/c if multiple genes were "nan", the column only contains one "nan" and you can't tell which genes correspond to the missing info 

