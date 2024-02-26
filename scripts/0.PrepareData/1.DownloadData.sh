####################################
# Download data 
####################################

# MGI database

## Create data directory
mkdir data/MGI

## Download data

### KO data
wget -P data/MGI https://www.informatics.jax.org/downloads/reports/MGI_PhenoGenoMP.rpt  https://www.informatics.jax.org/downloads/reports/MGI_PhenotypicAllele.rpt

### Mammalian Phenotype Ontology data (use different formats)
wget -P data/MGI https://www.informatics.jax.org/downloads/reports/MPheno_OBO.ontology  https://www.informatics.jax.org/downloads/reports/mp.json https://www.informatics.jax.org/downloads/reports/mp.owl

### Association of genes to diseases
wget -P data/MGI https://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt

### InterPRO domains in MGI proteins
wget -P data/MGI https://www.informatics.jax.org/downloads/reports/MGI_InterProDomains.rpt

### Mapping markers to other formats
wget -P data/MGI https://www.informatics.jax.org/downloads/reports/MRK_Sequence.rpt

# Interpro
mkdir data/Interpro
wget -P data/Interpro https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list
wget -P data/Interpro  https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/interpro.xml.gz
wget -P data/Interpro https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/match_complete.xml.gz

# Uniprot
wget -P data/Uniprot https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
wget -P data/Uniprot https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000000589/UP000000589_10090.dat.gz
# GO 
mkdir data/GO
wget -P data/GO https://current.geneontology.org/ontology/subsets/goslim_mouse.obo
wget -P data/GO https://current.geneontology.org/annotations/mgi.gaf.gz
wget -P data/GO https://purl.obolibrary.org/obo/go/go-basic.obo

# JENSEN database
## Create data directory
mkdir data/JENSEN

## Download compartment data (only knowledge channel as in Bioteque)
wget -P data/JENSEN https://download.jensenlab.org/mouse_compartment_knowledge_full.tsv

## Download tissue data (only knowledge channel as in Bioteque)
wget -P data/JENSEN https://download.jensenlab.org/mouse_tissue_knowledge_full.tsv

## Download tissue ontology from BRENDA
wget -O data/JENSEN/tissue.obo https://www.brenda-enzymes.org/ontology/tissue/tree/update/update_files/BrendaTissueOBO