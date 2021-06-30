install.packages('usethis')
install.packages('covr')
install.packages('httr')
install.packages('rversions')
install.packages('devtools')

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("multtest")

devtools::install_github('Seurat')
BiocManager::install("scran")