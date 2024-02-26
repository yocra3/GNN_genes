#'#################################################################################
#'#################################################################################
#' Process Tissues data
#'#################################################################################
#'#################################################################################

docker run -v $PWD:$PWD -w $PWD -it gnn_rsession:1.2 bash

## Load libraries
library(tidyverse)
library(ontologyIndex)
library(ontologyPlot)

## Load data
tissue <- read_delim("data/JENSEN/mouse_tissue_knowledge_full.tsv",
    delim = "\t", col_names = FALSE)
colnames(tissue) <- c("EnsemblID", "MarkerSymbol", "BTO_ID", "BTO_Name", "Source", 
    "EvidenceCode", "EvidenceLevel")

gene_map <- read_delim("results/preprocess/MGI_uniprot_map.tsv", delim = "\t")
load("results/preprocess/final_mice_genotype_phenotype.Rdata")

basic_tissue_raw <- get_ontology("data/JENSEN/tissue.obo")

## Merge both dataset
comb_df <- left_join(gene_map, tissue, by = "MarkerSymbol")
## Select compartment(s) with evidence equal or larger than 3
tissue_top <- comb_df %>% group_by(MGIMarkerAccessionID) %>% 
    filter(EvidenceLevel >= 3 ) %>%
    select(MGIMarkerAccessionID, MarkerSymbol, BTO_ID, BTO_Name, EvidenceLevel) %>%
    distinct()

tissue_map <- select(tissue_top, MGIMarkerAccessionID, BTO_ID)

## Filter GO tree
basic_tissue <- lapply(basic_tissue_raw, function(x) x[!basic_tissue_raw$obsolete])
class(basic_tissue) <- "ontology_index"

## Get list of markers
all_genos <- unlist(strsplit(genotype_phenotype$MGIMarkerGenotype, ","))
all_genes <- unlist(strsplit(all_genos, "/"))
all_genes <- unique(all_genes[grep("MGI", all_genes)])

## Apply the same approach than for GO BP (4.ProcessGOdata.R)
## Define function to get BTOs present in a number of genes (make more general)
get_BTO_IDs <- function(N, gene_ont_map, ontology){

    genes_per_bto <- gene_ont_map %>%
        group_by(BTO_ID) %>%
        summarize(n = n())

    sel_btos <- subset(genes_per_bto, n >= N)$BTO_ID
    map_sel <- subset(gene_ont_map, BTO_ID %in%  sel_btos)

    ## Filter tree to go terms with the required number of genes
    geneNum_mask <- ontology$id %in% sel_btos
    ontology_filt <- lapply(ontology, function(x) x[geneNum_mask])
    class(ontology_filt) <- "ontology_index"

    ## Filter parents, ancestors and children slots to GO present in our first filter
    ontology_filt$parents <- parallel::mclapply(ontology_filt$parents, 
                                function(x) x[x %in% sel_btos], mc.cores = 20)
    ontology_filt$children <-  parallel::mclapply(ontology_filt$children, 
                                function(x) x[x %in% sel_btos], mc.cores = 20)
    ontology_filt$ancestors <- parallel::mclapply(ontology_filt$ancestors, 
                                function(x) x[x %in% sel_btos], mc.cores = 20)

    ## Select GOs without children
    leaves_geneNum_filt <- ontology_filt$id[lengths(ontology_filt$children) == 0]
    leaves_go_gene <- subset(map_sel, BTO_ID %in% leaves_geneNum_filt)
    leaves_go_gene
}

genes_per_bto <- tissue_map %>%
        group_by(BTO_ID) %>%
        summarize(n = n())

## Select GO terms based on proportion of genes covered
total_genes <- length(all_genes)
props <- c(0.005, 0.01, 0.02, 0.05, 0.10)
names(props) <- c("0.5%", "1%", "2%", "5%", "10%")
n_genes <- props*total_genes

sel_tissue_N <- lapply(n_genes, get_BTO_IDs, 
    gene_ont_map = tissue_map, ontology = basic_tissue)

invisible <- lapply(names(sel_tissue_N), function(n){
    x <- sel_tissue_N[[n]]
    message("Proportion ", n, " contained ", length(unique(x$BTO_ID)), " BTO ids and ", 
    sum(all_genes %in% unique(x$MGIMarkerAccessionID)), " genes. ", round(mean(all_genes %in% unique(x$MGIMarkerAccessionID))*100, 2), 
    "% of original genes were present in any BTO term.")
})
save(sel_tissue_N, tissue_map, file = "results/preprocess/Tissue_gene_map_proportions.Rdata")
