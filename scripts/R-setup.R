## Step 1: Install Software ----

# Install R locally at version 3.4  https://cran.ma.imperial.ac.uk/
# Install Rstudio Desktop https://www.rstudio.com/products/rstudio/download/#download
# Install pandoc https://pandoc.org/installing.html

## Step 2: Extract Bioconductor packages ----

# Extract bioconductor_packages.zip to ~/Downloads

## Step 3: Open this script in Rstudio and run it ----

required_packages <- c("readr", "here", "glue", "dplyr", "tidyr", "stringr", "ggplot2",
                       "broom", "forcats", "testthat")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

bioconductor_packages <- c("Biostrings")
new_bioconducor <- bioconductor_packages[!(bioconductor_packages %in% installed.packages()[,"Package"])]
# if not proxy issues just install biostrings as normal from bioconductor

if(length(new_packages)){
  install.packages("~/Downloads/01_BiocGenerics_0.24.0.tar.gz")
  install.packages("~/Downloads/02_S4Vectors_0.16.0.tar.gz")
  install.packages("~/Downloads/03_IRanges_2.12.0.tar.gz")
  install.packages("~/Downloads/04_XVector_0.18.0.tar.gz")
  install.packages("~/Downloads/05_Biostrings_2.46.0.tar.gz")
}

# Touch wood that should be it!