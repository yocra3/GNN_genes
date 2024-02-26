#'#################################################################################
#'#################################################################################
#' Process MGI genotype data
#'#################################################################################
#'#################################################################################

docker run -v $PWD:$PWD -w $PWD -it gnn_rsession:1.2 bash

## Load libraries
library(tidyverse)

## Load data
allele <- read_delim("data/MGI/MGI_PhenotypicAllele.rpt", 
    comment = "#", 
    col_names = FALSE)
colnames(allele) <- c("MGIAlleleID", "AlleleSymbol", "AlleleName",
    	"AlleleType", "AlleleAttribute", "PubMedID", "MGIMarkerID",	
        "MarkerSymbol",	"MarkerRefSeqID", "MarkerEnsemblID", 
        "HighlevelMammalianPhenotypeID", "Synonyms", "MarkerName")

phenotypes <- read_delim("data/MGI/MGI_PhenoGenoMP.rpt", 
    col_names = FALSE)
colnames(phenotypes) <- c("AllelicComposition", "AlleleSymbol",
	"GeneticBackground", "MammalianPhenotypeID", "PubMedID",
    "MGIMarkerID")



## Select only KO alleles
ko_alleles <- subset(allele, 
    AlleleAttribute %in% c("Null/knockout|Reporter", "Reporter|Null/knockout", "Null/knockout"))

sel_alleles <- unique(ko_alleles$AlleleSymbol)

## Define wt alleles
wt_alleles <- paste0(unique(ko_alleles$MarkerSymbol), "<+>")

ko_alleles %>% 
    group_by(MarkerSymbol) %>%
    summarize(Function = any(!is.na(HighlevelMammalianPhenotypeID))) %>%
    group_by(Function) %>%
    summarize(n = n())

## Select phenotypes due to KO alleles
### Select enabled alleles (KO + WT)
good_alleles <- c(sel_alleles, wt_alleles)
### Convert Allelic Composition to list
phenotypes$AlleleList <- strsplit(phenotypes$AllelicComposition, ",|/")
good_phenotypes <- subset(phenotypes, 
    sapply(phenotypes$AlleleList, function(x) all(x %in% good_alleles)))

marker_len <- lengths(strsplit(good_phenotypes$MGIMarkerID, "|", fixed = TRUE))

## Remove combinations where length marker gene *2 != length alleles
good_phenotypes2 <- subset(good_phenotypes, 
    marker_len*2 == lengths(good_phenotypes$AlleleList))

## Define each marker genotype
mapAlleles <- function(vec){
    res <- sapply(vec, function(x) {

        if (x %in% sel_alleles){
          out <- subset(ko_alleles, AlleleSymbol == x)$MGIMarkerID
        } 
        else if (x %in% wt_alleles){
            out <- ""
        } else{
            out <- NA
        }
        return(out)
    })
    paste(res, collapse = ",")
}

good_phenotypes2 <- mutate(good_phenotypes2,
    MGIMarkerList = unlist(parallel::mclapply(AlleleList, mapAlleles, mc.cores = 10))
)

markerList <- lapply(strsplit(good_phenotypes2$MGIMarkerList, ","), function(x){
    sub <- x[x != ""]
    table(sub)
})
## Check that all the alleles were present in the KO or WT alleles
sum(sapply(markerList, function(x) any(is.na(x))))

createGenotypesVector <- function(x){
        x[order(names(x))] ## Order gene names alphabetically
        paste(paste(names(x), x, sep = "/"), collapse = ",") ## Combine all markers in a single vector
    }

good_phenotypes_final <- mutate(good_phenotypes2,
    MGIMarkerGenotype = sapply(markerList, createGenotypesVector)
)

phenotypes_summary <- group_by(good_phenotypes_final, AllelicComposition, GeneticBackground) %>%
    summarize(MGIMarkerGenotype = unique(MGIMarkerGenotype),
        PhenotypeID = paste(MammalianPhenotypeID, collapse = ","))

## Add KO genes without phenotype. Assume they were tested in homozygosis
no_effect_genes <- ko_alleles %>% 
    group_by(MGIMarkerID) %>%
    summarize(Function = any(!is.na(HighlevelMammalianPhenotypeID))) %>%
    filter(Function == FALSE)

no_effect_alleles <- filter(ko_alleles,
    MGIMarkerID %in%  no_effect_genes$MGIMarkerID)



phenotypes_summary_noeffect <- rbind(phenotypes_summary, 
    tibble(AllelicComposition = paste(no_effect_alleles$AlleleSymbol,no_effect_alleles$AlleleSymbol, sep = ""),
     GeneticBackground = "Unknown",
     MGIMarkerGenotype = paste(no_effect_alleles$MGIMarkerID, rep("/2", nrow(no_effect_alleles)),  sep = ""),
     PhenotypeID = ""
     )
)
dir.create("results/preprocess/")
save(phenotypes_summary, good_phenotypes_final, phenotypes,ko_alleles,
    file = "results/preprocess/MGI_phenotypes_genotypes_intermediates.Rdata")
save(phenotypes_summary_noeffect, file = "results/preprocess/MGI_phenotypes_genotypes_preprocess.Rdata")


