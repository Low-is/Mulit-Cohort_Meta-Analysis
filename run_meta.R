library(yaml)
library(jsonlite)

config <- yaml::read_yaml("Multi-Cohort_Meta-Analysis/config/config.yaml")

dna_studies <- jsonlite::fromJSON(config$analysis$input$dna_gse_file)
rna_studies <- jsonlite::fromJSON(config$analysis$input$rna_gse_file)

# Need to add code that filters matrices to match dimmensions of pData

# Find common genes across all studies being used for meta-analysis
common_genes <- find_common_genes(DNA = config$analysis$modalities$DNA,
                                  RNA = config$analysis$modalities$RNA,
                                  list_of_dna_mtx = matrices_dna[dna_studies],
                                  list_of_rna_mtx = matrices_RNA[rna_studies],
                                  use_DEG = config$analysis$use_DEG
                                 )
