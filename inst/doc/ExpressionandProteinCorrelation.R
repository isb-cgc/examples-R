#' # Expression and Protein Correlation
#' 
#' TODO introduction, talk about the tables we will use, etc...
#' 
#' Reproduces some of the analyses in http://watson.nci.nih.gov/~sdavis/tutorials/TCGA_data_integration/
#' 
## ----message=FALSE-------------------------------------------------------
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
#' ## Spearman Correlation in BigQuery
#' 
## ----comment=NA----------------------------------------------------------
# Set the desired tables to query.
expressionTable = "isb-cgc:tcga_data_open.mRNA_UNC_HiSeq_RSEM"
proteinTable = "isb-cgc:tcga_data_open.Protein"
cohortTable = "isb-cgc:test.cohort_14jun2015"

# Do not correlate unless there are at least this many observations available
minimumNumberOfObservations = 30

# Now we are ready to run the query.
result = DisplayAndDispatchQuery(file.path(sqlDir, "protein-mrna-spearman-correlation.sql"),
                                 project=project,
                                 replacements=list("_EXPRESSION_TABLE_"=expressionTable,
                                                   "_PROTEIN_TABLE_"=proteinTable,
                                                   "_COHORT_TABLE_"=cohortTable,
                                                   "_MINIMUM_NUMBER_OF_OBSERVATIONS_"=minimumNumberOfObservations))

#' Number of rows returned by this query: `r nrow(result)`.
#' 
#' The result is one correlation value per row of data, each of which corresponds to  . . . MORE HERE
## ------------------------------------------------------------------------
head(result)

#' 
## ----spearman_density, fig.align="center", fig.width=10, message=FALSE, warning=FALSE, comment=NA----
library(ggplot2)

# Histogram overlaid with kernel density curve
ggplot(result, aes(x=spearman_corr)) + 
    geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   binwidth=.05,
                   colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot

#' 
#' ## Spearman Correlation in R
#' 
#' Now let's reproduce one of the results directly in R.
#' 
#' ### Retrieve Expression Data
#' 
#' First we retrieve the expression data for a particular gene.
## ----comment=NA----------------------------------------------------------
# Set the desired gene to query.
gene = "ANXA1"

expressionData = DisplayAndDispatchQuery(file.path(sqlDir, "expression-data-by-cohort.sql"),
                                         project=project,
                                         replacements=list("_EXPRESSION_TABLE_"=expressionTable,
                                                           "_COHORT_TABLE_"=cohortTable,
                                                           "_GENE_"=gene))

#' Number of rows returned by this query: `r nrow(expressionData)`.
#' 
## ------------------------------------------------------------------------
head(expressionData)

#' 
#' ### Retrieve Protein Data
#' 
#' Then we retrieve the protein data for a particular gene.
#' 
## ----comment=NA----------------------------------------------------------
protein = "Annexin-1"

proteinData = DisplayAndDispatchQuery(file.path(sqlDir, "protein-data-by-cohort.sql"),
                                      project=project,
                                      replacements=list("_PROTEIN_TABLE_"=proteinTable,
                                                        "_COHORT_TABLE_"=cohortTable,
                                                        "_GENE_"=gene,
                                                        "_PROTEIN_"=protein))

#' Number of rows returned by this query: `r nrow(proteinData)`.
#' 
## ------------------------------------------------------------------------
head(proteinData)

#' 
#' ### Perform the correlation
#' 
## ------------------------------------------------------------------------
library(dplyr)

data = inner_join(expressionData, proteinData)
dim(data)
head(data)

#' 
#' First we take the inner join of this data and run a spearman correlation on it.
## ------------------------------------------------------------------------
cor(x=data$normalized_count, y=data$protein_expression, method="spearman")

#' Notice that the value does not match the result from BigQuery.
#' 
#' The reason for this is that the RANK in the BigQuery method is run over all the observations, not just those with a matching observation in both tables.  So let's redo this in R.
#' 
#' First we perform the RANK operation before doing the inner join.
## ------------------------------------------------------------------------
expressionData = mutate(expressionData, expr_rank=rank(normalized_count))
proteinData = mutate(proteinData, prot_rank=rank(protein_expression))

#' 
#' Then we do the join and run a spearman correlation on the ranked values.
## ------------------------------------------------------------------------
data = inner_join(expressionData, proteinData)
dim(data)
head(data)
cor(x=data$expr_rank, y=data$prot_rank, method="pearson")

#' Now the results match those for BigQuery.
#' 
#' But perhaps the BigQuery query should be running the RANK operation only after data has been removed per the inner join.
## ----comment=NA----------------------------------------------------------
resultV2 = DisplayAndDispatchQuery(file.path(sqlDir, "protein-mrna-spearman-correlation-V2.sql"),
                                 project=project,
                                 replacements=list("_EXPRESSION_TABLE_"=expressionTable,
                                                   "_PROTEIN_TABLE_"=proteinTable,
                                                   "_COHORT_TABLE_"=cohortTable,
                                                   "_MINIMUM_NUMBER_OF_OBSERVATIONS_"=minimumNumberOfObservations))

#' Number of rows returned by this query: `r nrow(result)`.
#' 
#' TODO This result does not yet match our earlier spearman correlation via R.
## ------------------------------------------------------------------------
head(resultV2, n=12)

#' 
## ----spearman_density2, fig.align="center", fig.width=10, message=FALSE, warning=FALSE, comment=NA----
library(ggplot2)

# Histogram overlaid with kernel density curve
ggplot(resultV2, aes(x=spearman_corr)) + 
    geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   binwidth=.05,
                   colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot

#' 
#' 
#' ## Provenance
## ----provenance, comment=NA----------------------------------------------
sessionInfo()

