library(yaml)
library(jsonlite)
library(COCONUT) # need to include installation of COCONUT in setup.R script

# Loading config file
message("Loading config file...")
config <- yaml::read_yaml("Multi-Cohort_Meta-Analysis/config/config.yaml")
message("Config file loaded!")

# Loading named list of studies
message("Loading list of studies...")
dna_studies <- jsonlite::fromJSON(config$analysis$input$dna_gse_file)
rna_studies <- jsonlite::fromJSON(config$analysis$input$rna_gse_file)
message("Studies loaded!")

# Loading expression matrices
message("Loading expression matrices...")
dna_matrices <-  readRDS("meta/matrices/dna_matrices.rds")
rna_matrices <- readRDS("meta/matrices/rna_matrices.rds")
message("Matrices loaded!")

# Need to add code that filters matrices to match dimmensions of pData

# Find common genes across all studies being used for meta-analysis
common_genes <- find_common_genes(DNA = config$analysis$modalities$DNA,
                                  RNA = config$analysis$modalities$RNA,
                                  list_of_dna_mtx = dna_matrices[dna_studies],
                                  list_of_rna_mtx = rna_matrices[rna_studies],
                                  use_DEG = config$analysis$use_DEG
                                 )

