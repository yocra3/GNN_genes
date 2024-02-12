#'#################################################################################
#'#################################################################################
#' Process MGI pheenotype data
#' Generate table with genotypes and phenotypes
#'#################################################################################
#'#################################################################################

docker run -v $PWD:$PWD -w $PWD -it gnn_rsession:1.1 bash

## Load libraries
library(tidyverse)
library(ontologyIndex)

## Load data
load("results/preprocess/MGI_phenotypes_genotypes_preprocess.Rdata")

mpo_ont <- get_ontology("data/MGI/MPheno_OBO.ontology")

## Select top level MPO terms
top_level <- mpo_ont$children[[1]]

## Select MPO terms present in the genotype tibble
all_phenos <- strsplit(phenotypes_summary_noeffect$PhenotypeID, ",") %>%
    unlist() %>% unique()

## Map phenos to top level MPO terms
getTopVec <- function(pheno){
    top_level %in% mpo_ont$ancestor[[which(mpo_ont$id == pheno)]]
}
phenos_top_map <- tibble(pheno = all_phenos) %>%
    mutate(top_pheno_vec = lapply(pheno, getTopVec))
phenos_top_map$top_pheno <- sapply(phenos_top_map$top_pheno_vec, function(x){
    paste(top_level[x], collapse = ",")
})

genotype_phenotype <- mutate(phenotypes_summary_noeffect,
    TopPhenotypeIDs = sapply(PhenotypeID, function(x){
        vec <- strsplit(PhenotypeID, ",")[[1]]
        all_top <- subset(phenos_top_map, pheno %in% vec)
        all_phenos <- unlist(strsplit(all_top$top_pheno, ","))
        un_phenos <- unique(all_phenos)
        out <- paste(un_phenos, collapse = ",")
    })
)