#'#################################################################################
#'#################################################################################
#' Prepare training data 
#' NN models 
#' Select all genes and all gene combination
#' Collapse gene features present in multiple genes
#' Allow gene combinations
#'#################################################################################
#'#################################################################################

docker run -v $PWD:$PWD -w $PWD -it gnn_rsession:1.2 bash

## Load libraries
library(tidyverse)
library(caret)

## Load data
load("results/preprocess/final_mice_genotype_phenotype.Rdata")
load("results/preprocess/Compartment_gene_map_proportions.Rdata")
load("results/preprocess/GO_MF_gene_map_proportions.Rdata")
load("results/preprocess/Tissue_gene_map_proportions.Rdata")
load("results/preprocess/GO_gene_map_proportions.Rdata")

## Get list of genes
genes_l <- gsub("\\/[1-2]", "", genotype_phenotype$MGIMarkerGenotype) %>%
    strsplit(",")

## Define functiosn to create matrices from feature data
createFeatureVector <- function(terms){
    out <- rep(0, length(terms))
    names(out) <- terms
    out
}

createFeatureMatrix <- function(gene_char, feature_df, feature_id){

    gene_vec <- strsplit(gene_char, ",")
    out <- createFeatureVector(unique(feature_df[[feature_id]]))

    ## Get genes and genos in a vector
    genes <- sapply(gene_vec, function(x) strsplit(x, split = "/") %>% sapply(FUN = `[`, 1))
    genos <- sapply(gene_vec, function(x) strsplit(x, split = "/") %>% sapply(FUN = `[`, 2))
    
    tab <- data.frame(MGIMarkerID = genes, geno = as.numeric(genos)) %>%
        merge(feature_df, by = "MGIMarkerID") %>%
        group_by(!!sym(feature_id)) %>%
        summarize(n = sum(geno))
    
    out[tab[[feature_id]]] <- tab$n
    out
}

## GO BP matrices
go_genes <- unique(genes_go_map_bp$MGIMarkerID)
go_bp_matrices <- lapply(sel_go_N, function(df) {

    ## Create a "new" GO indicating if the gene is present in any GO term
    gene_df <- data.frame(MGIMarkerID = go_genes, GO_ID = "Any")
    parallel::mclapply(genotype_phenotype$MGIMarkerGenotype, createFeatureMatrix, rbind(df, gene_df), "GO_ID", mc.cores = 20) %>%
        do.call(what = rbind)

}
)

## GO MF matrices
mf_genes <- unique(genes_go_map_mf$MGIMarkerID)
go_mf_matrices <- lapply(sel_mf_N, function(df) {

    ## Create a "new" GO indicating if the gene is present in any GO term
    gene_df <- data.frame(MGIMarkerID = mf_genes, GO_ID = "Any")
    parallel::mclapply(genotype_phenotype$MGIMarkerGenotype, createFeatureMatrix,  rbind(df, gene_df), "GO_ID", mc.cores = 20) %>%
        do.call(what = rbind)
}
)

## Tissue matrices
colnames(tissue_map)[1] <- "MGIMarkerID"
tissue_genes <- unique(tissue_map$MGIMarkerID)
tissue_matrices <- lapply(sel_tissue_N, function(df) {

    ## Create a "new" GO indicating if the gene is present in any GO term
    gene_df <- data.frame(MGIMarkerID = tissue_genes, BTO_ID = "Any")
    parallel::mclapply(genotype_phenotype$MGIMarkerGenotype, createFeatureMatrix,  rbind(df, gene_df), "BTO_ID", mc.cores = 20) %>%
        do.call(what = rbind)
}
)

## Compartment matrices
colnames(compartment_map)[1] <- "MGIMarkerID"
compartment_genes <- unique(compartment_map$MGIMarkerID)
compartment_matrices <- lapply(sel_compartment_N, function(df) {

    ## Create a "new" GO indicating if the gene is present in any GO term
    gene_df <- data.frame(MGIMarkerID = compartment_genes, GO_ID = "Any")
    parallel::mclapply(genotype_phenotype$MGIMarkerGenotype, createFeatureMatrix,  rbind(df, gene_df), "GO_ID", mc.cores = 20) %>%
        do.call(what = rbind)
}
)

## Feature matrices 
### Start by combining the matrices by proportion
combine_matrices <- function(comb_list, index) {
  # Extract matrices at the specified index from each list
  matrices <- lapply(comb_list, function(x) x[[index]])
  # Combine matrices by column
  combined_matrix <- do.call(cbind, matrices)
  return(combined_matrix)
}
comb_list <- list(go_bp_matrices, go_mf_matrices, tissue_matrices, compartment_matrices)
out_matrices <- lapply(seq_len(length(go_bp_matrices)), combine_matrices, comb_list = comb_list)
names(out_matrices) <- gsub("%", "", names(go_bp_matrices))

## Define phenotype column
### Only consider MP:0010768 (Mortality/Aging)
output <-  grepl("MP:0010768", genotype_phenotype$TopPhenotypeIDs) %>%
        as.numeric()

## Split the data keepìng unique combinations either on test or training
gene_comp_tab <- tibble(vec = genotype_phenotype$MGIMarkerGenotype, output = output) %>%
    distinct()

set.seed(27)
train_indices <- createDataPartition(gene_comp_tab$output, p = 0.78, list = FALSE)[,1]
train_combs <- unique(gene_comp_tab$vec[train_indices])

train_mask <- genotype_phenotype$MGIMarkerGenotype %in% train_combs
test_mask <- !genotype_phenotype$MGIMarkerGenotype %in% train_combs

## Create the output features file for each dataset
out_tibbles <- lapply(names(out_matrices), function(name){

    out_mat <- out_matrices[[name]]
    out_tibble <- cbind(out_mat, output) %>%
        as_tibble()
    
    train_tibble <- out_tibble[train_mask, ]
    write.table(train_tibble, 
        file = paste0("results/NNmodels/data/collapse_feats_", name, "_percen_train.tsv"),
        sep = "\t", row.names = FALSE, quote = FALSE)

    test_tibble <- out_tibble[test_mask, ]
    write.table(test_tibble, 
        file = paste0("results/NNmodels/data/collapse_feats_", name, "_percen_test.tsv"),
        sep = "\t", row.names = FALSE, quote = FALSE)
    out_tibble
})
names(out_tibbles) <- names(out_matrices)
save(out_tibbles, out_matrices, file = "results/NNmodels/data/collapse_feats.Rdata")

## Create train and test with the same combinations than full and sparse NN
gene_go_map <- sel_go_N[["0.5%"]]
gene_gos <- unique(gene_go_map$MGIMarkerID)

all_genos_l <- gsub("\\/[1-2]", "", genotype_phenotype$MGIMarkerGenotype) %>%
    strsplit(",")
genotype_phenotype_go <- genotype_phenotype[unlist(parallel::mclapply(all_genos_l, function(x) all(x %in% gene_gos), mc.cores = 20 )),]



load("results/NNmodels/data/processed_data.Rdata")
train_mask2 <- genotype_phenotype$MGIMarkerGenotype %in% train_combs

test_combs <- genotype_phenotype_go$MGIMarkerGenotype[!genotype_phenotype_go$MGIMarkerGenotype %in% train_combs]
test_mask2 <- genotype_phenotype$MGIMarkerGenotype %in% test_combs

out_tibbles_filt <- lapply(names(out_tibbles), function(name){

    out_tibble <- out_tibbles[[name]]
    
    train_tibble <- out_tibble[train_mask2, ]
    write.table(train_tibble, 
        file = paste0("results/NNmodels/data/collapse_feats_", name, "_percen_reduc_train.tsv"),
        sep = "\t", row.names = FALSE, quote = FALSE)

    test_tibble <- out_tibble[test_mask2, ]
    write.table(test_tibble, 
        file = paste0("results/NNmodels/data/collapse_feats_", name, "_percen_reduc_test.tsv"),
        sep = "\t", row.names = FALSE, quote = FALSE)
    out_tibble
})

## Explore a bit the features before training the model
pc_features <- lapply(out_matrices, prcomp)

df <- pc_features[[2]]$x %>% as_tibble() %>%
    mutate(output = output)
ggplot(df, aes(x = PC1, y = PC2, col = factor(output))) +
    geom_point()