# CS510MidtermProject
Content for the Midterm Coding Project for CS 510 (Fall 2021)

This program contains the research I am currently conducting with Dr. Patricia Lopes in analyzing differentially expressed (DE) genes in quail subjects that were subjected to different treatments. For our research/independent study, we are interested in measuring the immune response levels of female quails who were either paired with healthy or sick males. The DE analysis is conducted using BioConductor R packages and the builk of the analysis is contained within this R file/project. 

To run the project, simply run each R chunk. The final output is a formatted table of the top 20 DE genes and some plots. 

NOTE: IT IS IMPORTANT TO DOWNLOAD AND INSTALL ALL OF THE PACKAGES LISTED IN THE FIRST R CHUNK IN ORDER FOR THE PROGRAM TO RUN PROPERLY. The packages DESeq2 and DEGreport must be installed by the BiocManager, the package manager for Bioconductor packages. As such, THE FOLLOWING LINE OF CODE MAY PROVE HELPFUL IN INSTALLING THE BIOCONDUCTOR AND BIOCMANAGER PACKAGES: 

 > if (!requireNamespace("BiocManager", quietly = TRUE)) 
 >    install.packages("BiocManager") 
 > BiocManager::install("DESeq2")
 > BiocManager::install("DEGreport")
 > BiocManager::install("apeglm") # Required for computing the log2-shrinkage in one of the R chunks 
