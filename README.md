# Pleiotropy
Investigating why drosophila haven't bothered to evolve separate signaling pathways for development and innate immunity

## Script descriptions
dataset\_processing: reads in the RNAseq (currently as count matrices rather than raw reads) datasets from Schlamp (immune), Becker (embryogenesis), and Tarikere (oogenesis), normalizes them to transcripts per million, filters out genes with a mean count of <1 tpm, and saves count files to .RData files

clustering/clustering\_dev/clustering\_oog: heirarchical clustering to group genes into clusters

pleiotropy/pleiotropy\_dev/pleiotropy\_oog: summarizing clustering results for pleiotropic vs non-pleiotropic genes according to Williams et al. classifications

compare: compare clustering results from each dataset

classification: short script to generate excel spreadsheets for lit reviews on specific pleiotropic genes

kmeans: potentially to-be-implemented alternative to heiarchical clustering

## Literature cited
Williams, Alissa M., Thi Minh Ngo, Veronica E. Figueroa, and Ann T. Tate. 2023. “The Effect of Developmental Pleiotropy on the Evolution of Insect Immune Genes.” Genome Biology and Evolution 15 (3): 1–16. https://doi.org/10.1093/gbe/evad044.

Becker, Kolja, Alina Bluhm, Nuria Casas-Vila, et al. 2018. “Quantifying Post-Transcriptional Regulation in the Development of Drosophila Melanogaster.” Nature Communications 9 (1): 4970. https://doi.org/10.1038/s41467-018-07455-9.

Tarikere, Shreeharsha, Guillem Ylla, and Cassandra G Extavour. 2022. “Distinct Gene Expression Dynamics in Germ Line and Somatic Tissue during Ovariole Morphogenesis in Drosophila Melanogaster.” G3 Genes|Genomes|Genetics 12 (2): jkab305. https://doi.org/10.1093/g3journal/jkab305.

Schlamp, Florencia, Sofie Y. N. Delbare, Angela M. Early, Martin T. Wells, Sumanta Basu, and Andrew G. Clark. 2021. “Dense Time-Course Gene Expression Profiling of the Drosophila Melanogaster Innate Immune Response.” BMC Genomics 22 (1): 304. https://doi.org/10.1186/s12864-021-07593-3.

