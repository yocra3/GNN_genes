#'#################################################################################
#'#################################################################################
#' Process InterPRO domains data
#'#################################################################################
#'#################################################################################

docker run -v $PWD:$PWD -w $PWD -it gnn_rsession:1.2 bash

## Load libraries
library(tidyverse)

## Load data
domains <- read_delim("results/preprocess/UniprotPFAMDomains.tsv",
    delim = "\t", col_names = TRUE)
load("results/preprocess/final_mice_genotype_phenotype.Rdata")

## Get list of markers
all_genos <- unlist(strsplit(genotype_phenotype$MGIMarkerGenotype, ","))
all_genes <- unlist(strsplit(all_genos, "/"))
all_genes <- unique(all_genes[grep("MGI", all_genes)])

domain_map <- select(domains,  MGIMarkerAccessionID, Domain_ID, N)

## Apply the same approach than for GO BP (4.ProcessGOdata.R)
## Define function to get Domains present in a number of genes 
getDomainIDs <- function(N, domain_map){

    genes_per_domain <- domain_map %>%
        group_by(Domain_ID) %>%
        summarize(n = n())

    sel_domains <- subset(genes_per_domain, n >= N)$Domain_ID
    map_sel <- subset(domain_map, Domain_ID %in%  sel_domains)
    map_sel
}

genes_per_domain <- domain_map %>%
        group_by(Domain_ID) %>%
        summarize(n = n())

## Select GO terms based on proportion of genes covered
total_genes <- length(all_genes)
props <- c(0.005, 0.01, 0.02, 0.05, 0.10)
names(props) <- c("0.5%", "1%", "2%", "5%", "10%")
n_genes <- props*total_genes

sel_domain_N <- lapply(n_genes, getDomainIDs, 
    domain_map = domain_map)

invisible <- lapply(names(sel_domain_N), function(n){
    x <- sel_domain_N[[n]]
    message("Proportion ", n, " contained ", length(unique(x$Domain_ID)), " domains and ", 
    sum(all_genes %in% unique(x$MGIMarkerAccessionID)), " genes. ", round(mean(all_genes %in% unique(x$MGIMarkerAccessionID))*100, 2), 
    "% of original genes were present in any domain term.")
})
## Very few domains shared between proteins