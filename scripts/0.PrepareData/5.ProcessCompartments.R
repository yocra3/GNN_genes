#'#################################################################################
#'#################################################################################
#' Process compartments data
#'#################################################################################
#'#################################################################################

docker run -v $PWD:$PWD -w $PWD -it gnn_rsession:1.2 bash

## Load libraries
library(tidyverse)
library(ontologyIndex)
library(ontologyPlot)

## Load data
compartments <- read_delim("data/JENSEN/mouse_compartment_knowledge_full.tsv",
    delim = "\t", col_names = FALSE)
colnames(compartments) <- c("EnsemblID", "MarkerSymbol", "GO_ID", "GO_Name", "Source", 
    "EvidenceCode", "EvidenceLevel")

gene_map <- read_delim("results/preprocess/MGI_uniprot_map.tsv", delim = "\t")
load("results/preprocess/final_mice_genotype_phenotype.Rdata")

basic_tissue_raw <- get_ontology("data/GO/go-basic.obo")

## Filter GO tree
basic_go <- lapply(basic_go_raw, function(x) x[!basic_go_raw$obsolete])
class(basic_go) <- "ontology_index"

## Merge both dataset
comb_df <- left_join(gene_map, compartments, by = "MarkerSymbol")
## Select compartment(s) with evidence equal or larger than 3
compartment_top <- comb_df %>% group_by(MGIMarkerAccessionID) %>% 
    filter(EvidenceLevel >= 3 ) %>%
    select(MGIMarkerAccessionID, MarkerSymbol, GO_ID, GO_Name, EvidenceLevel) %>%
    distinct()

compartment_map <- select(compartment_top, MGIMarkerAccessionID, GO_ID)

## Get list of markers
all_genos <- unlist(strsplit(genotype_phenotype$MGIMarkerGenotype, ","))
all_genes <- unlist(strsplit(all_genos, "/"))
all_genes <- unique(all_genes[grep("MGI", all_genes)])


## Apply the same approach than for GO BP (4.ProcessGOdata.R)
## Define function to get GOs present in a number of genes (make more general)
getGOIDs <- function(N, gene_ont_map, ontology){

    genes_per_go <- gene_ont_map %>%
        group_by(GO_ID) %>%
        summarize(n = n())

    sel_gos <- subset(genes_per_go, n >= N)$GO_ID
    map_sel <- subset(gene_ont_map, GO_ID %in%  sel_gos)

    ## Filter tree to go terms with the required number of genes
    geneNum_mask <- ontology$id %in% sel_gos
    ontology_filt <- lapply(ontology, function(x) x[geneNum_mask])
    class(ontology_filt) <- "ontology_index"

    ## Filter parents, ancestors and children slots to GO present in our first filter
    ontology_filt$parents <- parallel::mclapply(ontology_filt$parents, 
                                function(x) x[x %in% sel_gos], mc.cores = 20)
    ontology_filt$children <-  parallel::mclapply(ontology_filt$children, 
                                function(x) x[x %in% sel_gos], mc.cores = 20)
    ontology_filt$ancestors <- parallel::mclapply(ontology_filt$ancestors, 
                                function(x) x[x %in% sel_gos], mc.cores = 20)

    ## Select GOs without children
    leaves_geneNum_filt <- ontology_filt$id[lengths(ontology_filt$children) == 0]
    leaves_go_gene <- subset(map_sel, GO_ID %in% leaves_geneNum_filt)
    leaves_go_gene
}

genes_per_go <- compartment_map %>%
        group_by(GO_ID) %>%
        summarize(n = n())

## Select GO terms based on proportion of genes covered
total_genes <- length(all_genes)
props <- c(0.005, 0.01, 0.02, 0.05, 0.10)
names(props) <- c("0.5%", "1%", "2%", "5%", "10%")
n_genes <- props*total_genes

sel_compartment_N <- lapply(n_genes, getGOIDs, 
    gene_ont_map = compartment_map, ontology = basic_go)

invisible <- lapply(names(sel_compartment_N), function(n){
    x <- sel_compartment_N[[n]]
    message("Proportion ", n, " contained ", length(unique(x$GO_ID)), " GO ids and ", 
    sum(all_genes %in% unique(x$MGIMarkerAccessionID)), " genes. ", round(mean(all_genes %in% unique(x$MGIMarkerAccessionID))*100, 2), 
    "% of original genes were present in any GO term.")
})
save(sel_compartment_N, compartment_map, file = "results/preprocess/Compartment_gene_map_proportions.Rdata")
