#'#################################################################################
#'#################################################################################
#' Process GO BP mice ontology
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

## Start processing Biological Processes - Only select genes involved_in
evidences_bp <- subset(evidences, Aspect == "P" & Qualifier == "involved_in")
message(sum(all_genes %in% unique(evidences_bp$MGIMarkerID)), " genes and ",
round(mean(all_genes %in% unique(evidences_bp$MGIMarkerID))*100, 2), "% of genes are present in BP GO terms.")

## Select GO terms with genes present in KO dataset
ko_gos <- subset(evidences_bp, MGIMarkerID %in% all_genes)$GO_ID %>% unique()
evidences_bp_ko <- subset(evidences_bp, GO_ID %in% ko_gos)

## Count genes per GO
genes_go_map_bp <- evidences_bp_ko %>%
    select(MGIMarkerID, GO_ID) %>%
    distinct() 
genes_per_go_bp <- genes_go_map_bp %>%
    group_by(GO_ID) %>%
    summarize(n = n())

## Select GO terms with at least 5 genes
sel_gos_bp <- subset(genes_per_go_bp, n >= 5)$GO_ID
genes_go_map_bp_sel <- subset(genes_go_map_bp, GO_ID %in%  sel_gos_bp)

message(sum(all_genes %in% unique(genes_go_map_bp$MGIMarkerID)), " genes and ",
round(mean(all_genes %in% unique(genes_go_map_bp$MGIMarkerID))*100, 2), "% of genes are present in any GO term.")


message(sum(all_genes %in% unique(genes_go_map_bp_sel$MGIMarkerID)), " genes and ",
round(mean(all_genes %in% unique(genes_go_map_bp_sel$MGIMarkerID))*100, 2), "% of genes are present in any GO term.")

## Filter GO tree
basic_go <- lapply(basic_go_raw, function(x) x[!basic_go_raw$obsolete])
class(basic_go) <- "ontology_index"

## Select go terms present in our first filter
geneNum_mask <- basic_go$id %in% sel_gos_bp
go_tree_bp_geneNum_filt <- lapply(basic_go, function(x) x[geneNum_mask])
class(go_tree_bp_geneNum_filt) <- "ontology_index"

## Filter parents, ancestors and children slots to GO present in our first filter
go_tree_bp_geneNum_filt$parents <- parallel::mclapply(go_tree_bp_geneNum_filt$parents, 
                                function(x) x[x %in% sel_gos_bp], mc.cores = 20)
go_tree_bp_geneNum_filt$children <-  parallel::mclapply(go_tree_bp_geneNum_filt$children, 
                                function(x) x[x %in% sel_gos_bp], mc.cores = 20)
go_tree_bp_geneNum_filt$ancestors <- parallel::mclapply(go_tree_bp_geneNum_filt$ancestors, 
                                function(x) x[x %in% sel_gos_bp], mc.cores = 20)

leaves_geneNum_filt <- go_tree_bp_geneNum_filt$id[lengths(go_tree_bp_geneNum_filt$children) == 0]
leaves_go_gene <- subset(genes_go_map_bp_sel, GO_ID %in% leaves_geneNum_filt)

message(sum(all_genes %in% unique(leaves_go_gene$MGIMarkerID)), " genes and ",
round(mean(all_genes %in% unique(leaves_go_gene$MGIMarkerID))*100, 2), "% of genes are present in any GO term.")

## Select GO terms that have genes only present in that GO term
leave_genes <- unique(leaves_go_gene$MGIMarkerID)

sing_genes <- leaves_go_gene %>%
    group_by(MGIMarkerID) %>%
    summarize(n = n()) %>%
    filter(n == 1) %>%
    pull(MGIMarkerID)

sing_gos <- unique(subset(leaves_go_gene, MGIMarkerID %in% sing_genes)$GO_ID)

sing_leaves_go_gene <- subset(leaves_go_gene, GO_ID %in% sing_gos)

message(sum(all_genes %in% unique(sing_leaves_go_gene$MGIMarkerID)), " genes and ",
round(mean(all_genes %in% unique(sing_leaves_go_gene$MGIMarkerID))*100, 2), "% of genes are present in any GO term.")

# Check overlap between gene sets
go_gene_list <- lapply(split(sing_leaves_go_gene, f = sing_leaves_go_gene$GO_ID), function(x) x$MGIMarkerID)

comp_matrix <- parallel::mclapply(go_gene_list, function(x) 
    sapply(go_gene_list, function(y) mean(x %in% y)), mc.cores = 20) %>%
    Reduce(f = rbind)
rownames(comp_matrix) <- names(go_gene_list)
comp_matrix2 <- comp_matrix
diag(comp_matrix2) <- 0

## Find GOs with > 80% of overlaps
bad_ind <- which(comp_matrix2 > 0.8, arr.ind = TRUE)
rep_go_tab <- tibble(GO1 = rownames(bad_ind), GO2 = colnames(comp_matrix)[bad_ind[, 2]]) 

lapply(unique(rep_go_tab$GO1), function(go){
    
    sel_terms <- c(go, subset(rep_go_tab, GO1 == go)$GO2)   
    plot_terms <- remove_links(basic_go, get_ancestors(basic_go, sel_terms), hard = TRUE)
    png(paste0("figures/", gsub(":", ".", go, fixed = TRUE), "_GO_sim.png"))
    onto_plot(basic_go, terms = plot_terms, fillcolor = (plot_terms %in% sel_terms) + 2) %>% 
        plot()
    dev.off()
})


## Filter allele combinations to those in selected GOs
all_genos_l <- gsub("\\/[1-2]", "", genotype_phenotype$MGIMarkerGenotype) %>%
    strsplit(",")
sing_leaves_genes <- unique(sing_leaves_go_gene$MGIMarkerID)
genotype_phenotype_go <- genotype_phenotype[unlist(parallel::mclapply(all_genos_l, function(x) all(x %in% sing_leaves_genes), mc.cores = 20 )),]
save(genotype_phenotype_go, file = "results/preprocess/final_mice_genotype_phenotype_go_filtered.Rdata")
save(sing_leaves_go_gene, file = "results/preprocess/GO_gene_map.Rdata")


## Select GO terms based on proportion of genes covered
total_genes <- length(all_genes)
props <- c(0.005, 0.01, 0.02, 0.05, 0.10)
names(props) <- c("0.5%", "1%", "2%", "5%", "10%")
sapply(props, function(x) sum(genes_per_go_bp$n > x*total_genes))
n_genes <- props*total_genes


## Define function to get GOs present in a number of genes
getGOIDs <- function(N){

    sel_gos <- subset(genes_per_go_bp, n >= N)$GO_ID
    map_sel <- subset(genes_go_map_bp, GO_ID %in%  sel_gos)

    ## Filter tree to go terms with the required number of genes
    geneNum_mask <- basic_go$id %in% sel_gos
    go_geneNum_filt <- lapply(basic_go, function(x) x[geneNum_mask])
    class(go_geneNum_filt) <- "ontology_index"

    ## Filter parents, ancestors and children slots to GO present in our first filter
    go_geneNum_filt$parents <- parallel::mclapply(go_geneNum_filt$parents, 
                                function(x) x[x %in% sel_gos], mc.cores = 20)
    go_geneNum_filt$children <-  parallel::mclapply(go_geneNum_filt$children, 
                                function(x) x[x %in% sel_gos], mc.cores = 20)
    go_geneNum_filt$ancestors <- parallel::mclapply(go_geneNum_filt$ancestors, 
                                function(x) x[x %in% sel_gos], mc.cores = 20)

    ## Select GOs without children
    leaves_geneNum_filt <- go_geneNum_filt$id[lengths(go_geneNum_filt$children) == 0]
    leaves_go_gene <- subset(map_sel, GO_ID %in% leaves_geneNum_filt)
    leaves_go_gene
}
sel_go_N <- lapply(n_genes, getGOIDs)

invisible <- lapply(names(sel_go_N), function(n){
    x <- sel_go_N[[n]]
    message("Proportion ", n, " contained ", length(unique(x$GO_ID)), " GO ids and ", 
    sum(all_genes %in% unique(x$MGIMarkerID)), " genes. ", round(mean(all_genes %in% unique(x$MGIMarkerID))*100, 2), 
    "% of original genes were present in any GO term.")
})
save(sel_go_N, genes_go_map_bp, file = "results/preprocess/GO_gene_map_proportions.Rdata")




## Restrict go terms to high confidence evidences (PMID: 24573882) (Not done. It removes relationships with parents)
good_evidences <- c("CURATED", "IDA", "TAS", "NAS")
sel_evidences <- subset(evidences, EvidenceCode %in% good_evidences)
sel_evidences_bp <- subset(sel_evidences, Aspect == "P")


sum(all_genes %in% unique(sel_evidences$MGIMarkerID))
sum(all_genes %in% unique(sel_evidences_bp$MGIMarkerID))

