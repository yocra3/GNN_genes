# Docker inheritance
FROM bioconductor/bioconductor_docker:RELEASE_3_18

# Install basic R packages
RUN R -e 'BiocManager::install(c("cowplot", "tidyverse", "GOfuncR", "parallel"))'

# Install package for loading ontologies (v.1.1)
RUN R -e 'BiocManager::install(c("ontologyIndex", "ontologySimilarity", "ontologyPlot"))'

# Install caret for machine learning functions (v.1.2)
RUN R -e 'BiocManager::install("caret")'
