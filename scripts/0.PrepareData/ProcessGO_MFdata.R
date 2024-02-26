#'#################################################################################
#'#################################################################################
#' Process GO MF mice ontology
#'#################################################################################
#'#################################################################################

docker run -v $PWD:$PWD -w $PWD -it gnn_rsession:1.2 bash

## Load libraries
library(tidyverse)
library(ontologyIndex)
library(ontologyPlot)

## Load data
evidences <- read_delim("data/GO/mgi.gaf.gz", 
    comment = "!", 
    delim = "\t",
    col_names = FALSE)
colnames(evidences) <- c("DB", "MGIMarkerID", "MarkerSymbol", "Qualifier", 
    "GO_ID", "DBReference", "EvidenceCode", "From", "Aspect", "MarkerName",
    "Synonyms", "GeneType", "Taxon", "Date", "AssignedBy", 
    "AnnotationExtension", "IsoformID")

basic_go_raw <- get_ontology("data/GO/go-basic.obo")

load("results/preprocess/final_mice_genotype_phenotype.Rdata")

## Get list of markers
all_genos <- unlist(strsplit(genotype_phenotype$MGIMarkerGenotype, ","))
all_genes <- unlist(strsplit(all_genos, "/"))
all_genes <- unique(all_genes[grep("MGI", all_genes)])

## Remove irrelevant columns of evidences
evidences <- evidences[, c("MGIMarkerID", "MarkerSymbol", 
                            "Qualifier", "GO_ID", "EvidenceCode", 
                            "From", "Aspect", "MarkerName",
                            "GeneType")]
message(sum(all_genes %in% unique(evidences$MGIMarkerID)), " genes and ",
round(mean(all_genes %in% unique(evidences$MGIMarkerID))*100, 2), "% of genes are present in any GO term.")

## Process Molecular Funcations - Only select genes enables
evidences_mf <- subset(evidences, Aspect == "F" & Qualifier == "enables")
message(sum(all_genes %in% unique(evidences_mf$MGIMarkerID)), " genes and ",
round(mean(all_genes %in% unique(evidences_mf$MGIMarkerID))*100, 2), "% of genes are present in MF GO terms.")

## Select GO terms with genes present in KO dataset
ko_gos <- subset(evidences_mf, MGIMarkerID %in% all_genes)$GO_ID %>% unique()
evidences_mf_ko <- subset(evidences_mf, GO_ID %in% ko_gos)

## Count genes per GO
genes_go_map_mf <- evidences_mf_ko %>%
    select(MGIMarkerID, GO_ID) %>%
    distinct() 
genes_per_go_mf <- genes_go_map_mf %>%
    group_by(GO_ID) %>%
    summarize(n = n())

## Filter GO tree
basic_go <- lapply(basic_go_raw, function(x) x[!basic_go_raw$obsolete])
class(basic_go) <- "ontology_index"

## Select GO terms based on proportion of genes covered
total_genes <- length(all_genes)
props <- c(0.005, 0.01, 0.02, 0.05, 0.10)
names(props) <- c("0.5%", "1%", "2%", "5%", "10%")
n_genes <- props*total_genes


## Apply the same approach than for GO BP (4.ProcessGOdata.R)
## Copy function to get GOs present in a number of genes from 5.ProcessCompartments.R
sel_mf_N <- lapply(n_genes, getGOIDs, 
    gene_ont_map = genes_go_map_mf, ontology = basic_go)


invisible <- lapply(names(sel_mf_N), function(n){
    x <- sel_mf_N[[n]]
    message("Proportion ", n, " contained ", length(unique(x$GO_ID)), " GO ids and ", 
    sum(all_genes %in% unique(x$MGIMarkerID)), " genes. ", round(mean(all_genes %in% unique(x$MGIMarkerID))*100, 2), 
    "% of original genes were present in any GO term.")
})
save(sel_mf_N, genes_go_map_mf, file = "results/preprocess/GO_MF_gene_map_proportions.Rdata")




## Restrict go terms to high confidence evidences (PMID: 24573882) (Not done. It removes relationships with parents)
good_evidences <- c("CURATED", "IDA", "TAS", "NAS")
sel_evidences <- subset(evidences, EvidenceCode %in% good_evidences)
sel_evidences_bp <- subset(sel_evidences, Aspect == "P")


sum(all_genes %in% unique(sel_evidences$MGIMarkerID))
sum(all_genes %in% unique(sel_evidences_bp$MGIMarkerID))

