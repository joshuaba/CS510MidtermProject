# CS510MidtermProject
Content for the Midterm Coding Project for CS 510 (Fall 2021)

This program is related to the research I am currently conducting with Dr. Patricia Lopes in analyzing differentially expressed (DE) genes in quail subjects that were subjected to different treatments. For our research/independent study, we are interested in measuring the immune response levels of female quails who were either paired with healthy or sick males. The DE analysis is conducted using BioConductor R packages and the builk of the analysis is contained within this R file/project.

The goal of this project is to conduct differential expression (DE) analysis on HEK293F mice cells that were subjected to different treatments in order to promote MOV10 gene expression, supress MOV10 gene expression, or maintain regular MOV10 gene expression levels (the control group). Through the DE analysis, the goal is to identify those genes whose expression levels change as a result of the MOV10 gene perturbation (or expression level change). We can then conclude that there is a correlation between MOV10 expression levels and expression levels of these subsequent genes. 

To run the project, simply run each R chunk. The final output is a formatted table of the top 20 DE genes and some plots, in addition to plots to visualize the results. 

The project directory also includes a formal scientific report located in the "results" folder. 

NOTE: IT IS IMPORTANT TO DOWNLOAD AND INSTALL ALL OF THE PACKAGES LISTED IN THE FIRST R CHUNK IN ORDER FOR THE PROGRAM TO RUN PROPERLY. The packages DESeq2 and DEGreport must be installed by the BiocManager, the package manager for Bioconductor packages. As such, THE FOLLOWING LINES OF CODE MAY PROVE HELPFUL IN INSTALLING THE BIOCONDUCTOR AND BIOCMANAGER PACKAGES: 

 > if (!requireNamespace("BiocManager", quietly = TRUE)) 
 >    install.packages("BiocManager")^@ 
 > BiocManager::install("DESeq2")^@
 > BiocManager::install("DEGreport")^@
 > BiocManager::install("apeglm") # Required for computing the log2-shrinkage in one of the R chunks

 Alternatively, uncomment the line 
 > # source("install_script.R", local = knitr::knit_global())

 in the first R code chunk, which will automatically install the BiocManager package manager and the required Bioinformatics packages required for the project. 
