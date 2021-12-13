# R Install Script to Install Required Packages/Libraries 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install() # Install core packages
BiocManager::install("DESeq2")
BiocManager::install("DEGreport")
BiocManager::install("apeglm") # Required to compute the log 2 shrinkage 

