#'#################################################################################
#'#################################################################################
#' Map MGI ids of genes in models to Uniprot and Symbol
#'#################################################################################
#'#################################################################################

docker run -v $PWD:$PWD -w $PWD -it gnn_rsession:1.2 bash

## Load libraries
library(tidyverse)

## Load data
gene_mapping <- read_delim("data/MGI/MRK_Sequence.rpt", 
    comment = "#", 
    col_names = TRUE)
colnames(gene_mapping) <- gsub(" ", "", colnames(gene_mapping))

load("results/preprocess/GO_gene_map.Rdata")
load("results/preprocess/MGI_phenotypes_genotypes_preprocess.Rdata")


## Map genes to Uniprot IDs
ko_genes <- gsub("\\/[1-2]", "", phenotypes_summary_noeffect$MGIMarkerGenotype) %>%
    strsplit(",") %>% unlist() 
go_genes <- sing_leaves_go_gene$MGIMarkerID

all_genes <- unique(c(ko_genes, go_genes))


uniprot_map <- gene_mapping %>%
    filter(MGIMarkerAccessionID %in% all_genes) %>%
    select(MGIMarkerAccessionID, MarkerSymbol, UniProtIDs) %>%
    separate_rows(UniProtIDs, sep = "\\|")

write.table(uniprot_map, file = "results/preprocess/MGI_uniprot_map.tsv",
    sep = "\t", row.names = FALSE, quote = FALSE)