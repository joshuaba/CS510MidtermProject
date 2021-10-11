
## Gene-level differential expression analysis using DESeq2

## Setup
### Bioconductor and CRAN libraries used
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(DEGreport)

## Load in data
data <- read.table("data/Mov10_full_counts.txt", header=T, row.names=1) 

meta <- read.table("meta/Mov10_full_meta.txt", header=T, row.names=1)

### Check classes of the data we just brought in
class(meta)
class(data)

# Viewing the data to double-check before performing analysis 
View(meta)
View(data)

# Let us plot the RNA sequence 

ggplot(data) +
  geom_histogram(aes(x = Mov10_oe_1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

ggplot(data) +
  geom_histogram(aes(x = Mov10_oe_1), stat = "bin", bins = 200) + 
  xlim(-5, 500)  +
  xlab("Raw expression counts") +
  ylab("Number of genes")


# In order to model the RNA-sequences/count data, let us use the negative binomial (NB) distribution, as this model will best fit the actual count data present 
# in our RNA sequences, since mean(count_data) < variance(count_data). Although a Poisson distribution would normally model this type of RNA sequence count
# very well, it does not work well in this case since mean(count_data) != variance(count_data)

# Let's visualize the comments above below: 

# create a data frame with the mean of the RNA sequence counts and the variance of the RNA sequence counts 

mean_counts <- apply(data[, 3:5], 1, mean)
variance_counts <- apply(data[, 3:5], 1, var)
df <- data.frame(mean_counts, variance_counts)

# Plot the data frame created above 

ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  geom_line(aes(x=mean_counts, y=mean_counts, color="red")) +
  scale_y_log10() +
  scale_x_log10()

# The ggplot created above verifies the assertion made earlier that variance tends to be greater than mean RNA sequence counts 



# Raw counts for PD1
PD1 <- c(21, 58, 17, 97, 83, 10)
names(PD1) <- paste0("Sample", 1:6)
PD1 <- data.frame(PD1)
PD1 <- t(PD1)

# Size factors for each sample
size_factors <- c(1.32, 0.70, 1.04, 1.27, 1.11, 0.85)

#TODO: Compute the normalized values given the normalization (size) factors and raw counts. Divide each raw count by the respective normalization (size) factor for the sample
normalized_values <- PD1/size_factors
print(normalized_values)


### Check that sample names match in both files
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

## Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)


# Now that we have created our DESeq2Dataset object, let us normalize the count values by calling the estimateSizeFactors function (see below)
dds <- estimateSizeFactors(dds)

# Let's make sure that the normalized factors are in our data frame 
sizeFactors(dds)

# Let us compute the normalized RNA sequence counts 
normalized_counts <- counts(dds, normalized=TRUE)

# Let's write our normalized values out to a file to save it for further use
write.table(normalized_counts, file="data/normalized_counts.txt", sep="\t", quote=F, col.names=NA)

# Let us now moderate the variance for our normalized counts in order to improve the distances/clustering so we can properly perform/conduct principal component analysis (PCA) 
# and hierarchical clustering

#rlog() function will transform the normalized counts 

rld <- rlog(dds, blind=TRUE)

# Now that we have moderated the variance, let us run the PCA on our normalized RNA sequence counts

plotPCA(rld, intgroup="sampletype")

# Let us now perform hierarchical clustering. Since there is no HC function in R, let us use a heatmap 

### Extract the rlog matrix from the object
rld_mat <- assay(rld)    ## assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2


### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function

head(rld_cor)   ## check the output of cor(), make note of the rownames and colnames

### Plot heatmap
pheatmap(rld_cor)

# Both the heat map constructed in the line above and the PCA analysis performed earlier show the samples properly clustering together according to sample group. 
# This is good indication that our sample data is properly cleaned and we can proceed to DE analysis. 
# Relatively high correlations in our heat map (>0.999) also indicate no outlying samples. 

## Create DESeq object
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)


## Run analysis
dds <- DESeq(dds)

resultsNames(dds)

## Plot dispersion estimates
plotDispEsts(dds)

# We will now perform our hypothesis testing using a Walfowitz test. We need to first define the contrasts, or the variables we will be contrasting in our hypothesis test. 

# Define contrasts, extract results table, and shrink the log2 fold changes

contrast_oe <- c("sampletype", "MOV10_overexpression", "control")

res_tableOE_unshrunken <- results(dds, contrast=contrast_oe, alpha = 0.05)

####### Consult with Dr. Lopes about "type" of lfcShrink to use. Default not working #######
res_tableOE <- lfcShrink(dds, contrast=contrast_oe, res=res_tableOE_unshrunken, type="ashr")


# Let's generate an MA Plot, which will display the means of the normalized genes vs. the log2 foldchanges for MOV10_OE vs. the control group
plotMA(res_tableOE_unshrunken, ylim=c(-2,2))

# Let's view an MA plot of the shrunken results 
plotMA(res_tableOE, ylim=c(-2,2))

class(res_tableOE)

# Let's identify the names of the cols stored in our res_tableOE variable (the DESeq results table/dataframe)
mcols(res_tableOE, use.names=T)


# We'll take a look at our data frame 
res_tableOE
# Alternative Code: res_tableOE %>% data.frame() %>% View()


# ---------------------------- #


# Let's perform a Wald test and DE analysis for the MOV10_KD (knockdown) vs. the control group 

## Define contrasts, extract results table and shrink log2 fold changes
contrast_kd <-  c("sampletype", "MOV10_knockdown", "control")

res_tableKD_unshrunken <- results(dds, contrast=contrast_kd, alpha = 0.05)

res_tableKD <- lfcShrink(dds, contrast=contrast_kd, res=res_tableKD_unshrunken, type="ashr")

# Let's generate an MA Plot, which will display the means of the normalized genes vs. the log2 foldchanges for MOV10_OE vs. the control group
plotMA(res_tableKD_unshrunken, ylim=c(-2,2))

# Let's view an MA plot of the shrunken results 
plotMA(res_tableKD, ylim=c(-2,2))

class(res_tableKD)

mcols(res_tableKD, use.names=TRUE)

res_tableKD


# Let's extract DE genes for overexpression 

## Summarize results
summary(res_tableOE)
summary(res_tableKD)

### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.58


res_tableOE_tb <- res_tableOE %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()


sigOE <- res_tableOE_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)


sigOE

# ------------------------------------------------ #

# Let's extract DE genes for overexpression 

res_tableKD_tb <- res_tableKD %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sigKD <- res_tableKD_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)


sigKD

# TODO: Visualizing results 



