#'#################################################################################
#'#################################################################################
#' Prepare training data 
#' NN models
#' Select genes present in GOs selected with 0.5%
#' No include gene features
#' Allow gene combinations
#'#################################################################################
#'#################################################################################

docker run -v $PWD:$PWD -w $PWD -it gnn_rsession:1.2 R

## Load libraries
library(tidyverse)
library(caret)

## Load data
load("results/preprocess/GO_gene_map_proportions.Rdata")
load("results/preprocess/final_mice_genotype_phenotype.Rdata")


## Subset alleles combinations to genes present in GO 0.5%
gene_go_map <- sel_go_N[["0.5%"]]
gene_gos <- unique(gene_go_map$MGIMarkerID)

all_genos_l <- gsub("\\/[1-2]", "", genotype_phenotype$MGIMarkerGenotype) %>%
    strsplit(",")
genotype_phenotype_go <- genotype_phenotype[unlist(parallel::mclapply(all_genos_l, function(x) all(x %in% gene_gos), mc.cores = 20 )),]


## Create gene matrix
all_genes <- gsub("\\/[1-2]", "", genotype_phenotype_go$MGIMarkerGenotype) %>%
    strsplit(",") %>% unlist() %>% unique()

createGenesVector <- function(){
    out <- rep(0, length(all_genes))
    names(out) <- all_genes
    out
}

gene_mat <- sapply(strsplit(genotype_phenotype_go$MGIMarkerGenotype, ","), function(l){
    
    out <- createGenesVector()

    ## Get genes and genos in a vector
    genes <- sapply(l, function(x) strsplit(x, split = "/") %>% sapply(FUN = `[`, 1))
    genos <- sapply(l, function(x) strsplit(x, split = "/") %>% sapply(FUN = `[`, 2))
    
    
    out[genes] <- as.numeric(genos)/2
    out

}) %>% t()
colnames(gene_mat) <- all_genes

## Define phenotype column
### Only consider MP:0010768 (Mortality/Aging)
pheno_tab <- genotype_phenotype_go %>%
    mutate(Output = grepl("MP:0010768", TopPhenotypeIDs) %>%
        as.numeric()
    ) %>%
    ungroup() %>%
    select(MGIMarkerGenotype, Output)

## Create output file
out_tibble <- cbind(as_tibble(gene_mat), Output = pheno_tab$Output) %>%
    as_tibble()

## Split the data keepìng unique combinations either on test or training
set.seed(27)
train_indices <- createDataPartition(pheno_tab$Output, p = 0.60, list = FALSE)[,1]
train_combs <- unique(pheno_tab$MGIMarkerGenotype[train_indices])

train_mask <- pheno_tab$MGIMarkerGenotype %in% train_combs
train_df <- out_tibble[train_mask, ]
test_df <-  out_tibble[!train_mask, ]

write.table(train_df, file = "results/NNmodels/data/train.tsv",
    sep = "\t", row.names = FALSE, quote = FALSE)

write.table(test_df, file = "results/NNmodels/data/test.tsv",
    sep = "\t", row.names = FALSE, quote = FALSE)

write.table(colnames(gene_mat), file = "results/NNmodels/data/gene_order.tsv",
    sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

write.table(gene_go_map, file = "results/NNmodels/data/gene_go_map.tsv",
    sep = "\t", row.names = FALSE, quote = FALSE)

save(out_tibble, gene_mat, train_df, test_df, train_combs, file = "results/NNmodels/data/processed_data.Rdata")

