#' # Expression and Methylation Correlation
#' 
#' In this example, we will look at the correlation between mRNAseq-based gene expression and DNA methylation data.  We will do this using two molecular data tables from the isb-cgc:tcga_201510_alpha dataset and a cohort table from the isb-cgc:tcga_cohorts dataset.
#' 
#' NOTE: I think I will rework and/or eliminate this particular example, but am just going through it now to make sure I understand it and it works as expected.
#' 
## ----message=FALSE-------------------------------------------------------
library(dplyr)
library(bigrquery)
library(scales)
library(ggplot2)
library(ISBCGCExamples)

# The directory in which the files containing SQL reside.
#sqlDir = file.path("/PATH/TO/GIT/CLONE/OF/examples-R/inst/",
sqlDir = file.path(system.file(package = "ISBCGCExamples"),
                   "sql")

#' 
## ----eval=FALSE----------------------------------------------------------
## ######################[ TIP ]########################################
## ## Set the Google Cloud Platform project id under which these queries will run.
## ##
## ## If you are using the workshop docker image, this is already
## ## set for you in your .Rprofile and you can skip this step.
## 
## # project = "YOUR-PROJECT-ID"
## #####################################################################

#' 
#' ## Pearson Correlation in BigQuery
#' 
## ----comment=NA----------------------------------------------------------
# Set the desired tables to query.
expressionTable = "isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM"
methylationTable = "isb-cgc:tcga_201510_alpha.DNA_Methylation_betas"
# Add any additional clauses to be applied in WHERE to limit the methylation data further.
# (These specific filters are used here just to make the query run faster.  If a query returns
#  very large results, they may need to be handled differently.  This query should take < 20s)
andWhere = "AND SampleTypeLetterCode = 'TP' AND Study = 'CESC' AND CHR = '9'"
# Do not correlate unless there are at least this many observations available:
minNumObs = 30

# Now we are ready to run the query.  (Should return 6110 rows.)
result = DisplayAndDispatchQuery(
     file.path(sqlDir, "expression-methylation-correlation.sql"),
               project=project,
               replacements=list("_EXPRESSION_TABLE_"=expressionTable,
                                 "_METHYLATION_TABLE_"=methylationTable,
                                 "_AND_WHERE_"=andWhere,
                                 "_MINIMUM_NUMBER_OF_OBSERVATIONS_"=minNumObs))
cat("Number of rows returned by this query: ", nrow(result), "\n")

#' 
#' The result is a table with one row for each (gene,CpG-probe) pair for which at least 30 data values exist that meet the requirements in the "andWhere" clause.  The (gene,CpG-probe) pair is defined by a gene symbol and a CpG-probe ID.  In many cases, there may be multiple CpG probes associated with a single gene.
#' 
## ------------------------------------------------------------------------
# Most negative correlation should be PHYHD1 cg14299940 n=903, correlation = -0.8018487
head(result)

#' 
## ----density, fig.align="center", fig.width=10, message=FALSE, warning=FALSE, comment=NA----

# Histogram overlaid with kernel density curve
ggplot(result, aes(x=correlation)) +
    geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   binwidth=.05,
                   colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot

#' 
#' ## Pearson Correlation in R
#' 
#' Now let's reproduce one of the results directly in R.
#' 
#' ### Retrieve Expression Data
#' 
#' First we retrieve the expression data for a particular gene.
## ----comment=NA----------------------------------------------------------
# Set the desired gene to query.
gene = "PHYHD1"
andWhere = "AND SampleTypeLetterCode = 'TP' AND Study = 'CESC'"
expressionData = DisplayAndDispatchQuery(file.path(sqlDir, "expression-data.sql"),
                                         project=project,
                                         replacements=list("_EXPRESSION_TABLE_"=expressionTable,
                                                           "_GENE_"=gene,
                                                           "_AND_WHERE_"=andWhere))
cat("Number of rows returned by this query: ", nrow(expressionData), "\n")

#' 
#' 
## ------------------------------------------------------------------------
head(expressionData)

#' 
#' ### Retrieve Methylation Data
#' 
#' Then we retrieve the methylation data for a particular probe.
#' 
## ----comment=NA----------------------------------------------------------
# Set the desired probe to query.
probe = "cg14299940"
andWhere = "AND SampleTypeLetterCode = 'TP' AND Study = 'CESC'"
methylationData = DisplayAndDispatchQuery(file.path(sqlDir, "methylation-data.sql"),
                                          project=project,
                                          replacements=list("_METHYLATION_TABLE_"=methylationTable,
                                                            "_AND_WHERE_"=andWhere,
                                                            "_PROBE_"=probe))
cat("Number of rows returned by this query: ", nrow(methylationData), "\n")

#' 
## ------------------------------------------------------------------------
head(methylationData)

#' 
#' ### Perform the correlation
#' 
#' First we take the inner join of this data:
## ------------------------------------------------------------------------
data = inner_join(expressionData, methylationData)
head(data)

#' 
#' And run a pearson correlation on it:
## ------------------------------------------------------------------------
p = round(cor(x=log2(data$normalized_count+1), y=data$Beta_Value, method="pearson"), 3)
qplot(data=data, y=log2(normalized_count), x=Beta_Value, geom=c("point","smooth"),
      xlab="methylation level (beta value)", ylab="mRNA level") +
      geom_text(x = 0.7, y = 11, label = paste("Pearson Corr ", p))
p

#' 
#' And we can see that we have reproduced one of our results from BigQuery.
#' the correlation now is -0.802 on 301 samples...
#' 
#' ## Provenance
## ----provenance, comment=NA----------------------------------------------
sessionInfo()

