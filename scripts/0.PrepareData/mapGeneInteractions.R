#'#################################################################################
#'#################################################################################
#' Map OmniPath gene interactions to our data
#'#################################################################################
#'#################################################################################

docker run -v $PWD:$PWD -w $PWD -it gnn_rsession:1.2 bash

## Load libraries
library(tidyverse)
library(OmnipathR)

## Load interactions
post_transcrip <- import_post_translational_interactions(organism = "10090")
transcrip <- import_transcriptional_interactions(organism = "10090")

all_interactions <- rbind(post_transcrip, transcrip[, colnames(post_transcrip)])

## Load mice data
load("results/preprocess/GO_gene_map.Rdata")
gene_map <- read_delim("results/preprocess/MGI_uniprot_map.tsv", delim = "\t")


## Merge all datasets
all_comb <- left_join(sing_leaves_go_gene, select(gene_map, MGIMarkerAccessionID, MarkerSymbol), 
    by = join_by(MGIMarkerID == MGIMarkerAccessionID))

all_gos <- unique(all_comb$GO_ID)
gene_interaction_count <- sapply(all_gos[1:10], function(go){

    go_genes <- subset(all_comb, GO_ID == go)$MarkerSymbol
    go_inter <- subset(all_interactions, source_genesymbol %in% go_genes & target_genesymbol %in% go_genes)
    c(length(go_genes), sum(go_genes %in% c(go_inter$source_genesymbol, go_inter$target_genesymbol)))

})
## Very few connections