####################################
# Download data from MGI database
####################################

## Create data directory
mkdir data/MGI

## Download data

### KO data
wget -P data/MGI https://www.informatics.jax.org/downloads/reports/MGI_PhenoGenoMP.rpt  https://www.informatics.jax.org/downloads/reports/MGI_PhenotypicAllele.rpt

### Mammalian Phenotype Ontology data (use different formats)
wget -P data/MGI https://www.informatics.jax.org/downloads/reports/MPheno_OBO.ontology  https://www.informatics.jax.org/downloads/reports/mp.json https://www.informatics.jax.org/downloads/reports/mp.owl

### Association of genes to diseases
wget -P data/MGI https://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt