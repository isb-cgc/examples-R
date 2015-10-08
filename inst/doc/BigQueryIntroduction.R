#' # Analyzing Variants with BigQuery
#' 
#' TODO: this is a copy of the introduction used for BioC2015.  It should be re-worked to be more appropriate for ISB-CGC.
#' 
#' 
#' Google Genomics can import variant calls from VCF files or Complete Genomics masterVar files so that you can query them with a simple API as we saw earlier in vignette [Working with Variants](http://bioconductor.org/packages/devel/bioc/vignettes/GoogleGenomics/inst/doc/AnnotatingVariants.html).  You can also export Variants to BigQuery for interactive analysis of these large datasets.  For more detail, see https://cloud.google.com/genomics/v1/managing-variants
#' 
#' In this example we will work with the [Illumina Platinum Genomes](http://googlegenomics.readthedocs.org/en/latest/use_cases/discover_public_data/platinum_genomes.html) dataset.
#' 
#' ### Helpful BigQuery links
#' 
#' For this example, we'll also be working with [Google BigQuery](https://cloud.google.com/bigquery/). It's often helpful to have a [link to the docs](https://cloud.google.com/bigquery/what-is-bigquery) handy, and especially the [query reference](https://cloud.google.com/bigquery/query-reference).
#' 
#' ## Run a query from R
#' 
#' The [bigrquery](https://github.com/hadley/bigrquery) package written by Hadley Wickham implements an R interface to [Google BigQuery](https://cloud.google.com/bigquery/).
#' 
## ----message=FALSE-------------------------------------------------------
library(dplyr)
library(bigrquery)

#' 
## ----eval=FALSE----------------------------------------------------------
## ######################[ TIP ]########################################
## ## Set the Google Cloud Platform project id under which these queries will run.
## ##
## ## If you are using the Google Bioconductor workshop docker image, this is already
## ## set for you in your .Rprofile and you can skip this step.
## 
## # project <- "YOUR-PROJECT-ID"
## #####################################################################

#' 
## ------------------------------------------------------------------------
# Change the table here if you wish to run these queries against a different Variants table.
theTable <- "genomics-public-data:platinum_genomes.variants"
# theTable <- "genomics-public-data:1000_genomes.variants"
# theTable <- "genomics-public-data:1000_genomes_phase_3.variants"

#' 
#' Let's start by just counting the number of records in the table:
## ------------------------------------------------------------------------
querySql <- paste("SELECT COUNT(1) FROM [", theTable, "]", sep="")
querySql

#' 
#' And send the query to the cloud for execution:
## ------------------------------------------------------------------------
result <- query_exec(querySql, project=project)

#' 
#' [bigrquery](https://github.com/hadley/bigrquery) uses the package [httr](https://github.com/hadley/httr) to perform OAuth.
#' 
## ----eval=FALSE----------------------------------------------------------
## ######################[ TIP ]########################################
## ## If you have any trouble with OAuth and need to redo/reset OAuth,
## ## run the following code.
## 
## # if(FALSE != getOption("httr_oauth_cache")) {
## #  file.remove(getOption("httr_oauth_cache"))
## #}
## #message("Restart R to redo/reset OAuth.")
## #####################################################################

#' 
#' And we see that the table has `r result[1,1]` rows - wow!
## ------------------------------------------------------------------------
result

#' 
#' ## Run a query using the BigQuery Web User Interface
#' 
#' So what is actually in this table?  Click on [this link](https://bigquery.cloud.google.com/table/genomics-public-data:platinum_genomes.variants) to view the schema in the BigQuery web user interface.
#' 
#' We can also run the exact same query using the BigQuery web user interface.  In the BigQuery web user interface:
#' 
#'  (1) click on the *"Compose Query"* button
#'  (2) paste in the SQL for the query we just ran via R
#'  (3) click on *"Run Query"*.
#' 
#' ## Run a query stored in a file from R
#' 
#' Instead of typing SQL directly into our R code, we can use a convenience function to read SQL from a file.
## ------------------------------------------------------------------------
library(ISBCGCExamples)
DisplayAndDispatchQuery

#' 
#' This allows queries to be more easily shared among analyses and also reused for different datasets.  For example, in the following file we have a query that will retrieve data from any table exported from a Google Genomics Variant Set.
## ------------------------------------------------------------------------
file.show(file.path(system.file(package = "ISBCGCExamples"),
                    "sql",
                    "variant-level-data-for-brca1.sql"))

#' 
#' Now let's run the query to retrieve variant data for BRCA1:
## ----comment=NA----------------------------------------------------------
result <- DisplayAndDispatchQuery(file.path(system.file(package = "ISBCGCExamples"),
                                            "sql",
                                            "variant-level-data-for-brca1.sql"),
                                  project=project,
                                  replacements=list("_THE_TABLE_"=theTable))

#' Number of rows returned by this query: `r nrow(result)`.
#' 
#' Results from [bigrquery](https://github.com/hadley/bigrquery) are dataframes:
## ------------------------------------------------------------------------
mode(result)
class(result)
summary(result)
head(result)

#' 
#' ## Visualize Query Results
#' 
#' The prior query was basically a data retrieval similar to what we performed earlier in this workshop when we used the GoogleGenomics Bioconductor package to retrieve data from the Google Genomics Variants API.  
#' 
#' But BigQuery really shines when it is used to perform an actual *analysis* - do the heavy-lifting on the big data resident in the cloud, and bring back the result of the analysis to R for further downstream analysis and visualization.
#' 
#' Let's do that now with a query that computes the Transition Transversion ratio for the variants within genomic region windows.
## ----comment=NA----------------------------------------------------------
result <- DisplayAndDispatchQuery(file.path(system.file(package = "ISBCGCExamples"),
                                            "sql",
                                            "ti-tv-ratio.sql"),
                                  project=project,
                                  replacements=list("_THE_TABLE_"=theTable,
                                                    "_WINDOW_SIZE_"=100000))

#' Number of rows returned by this query: `r nrow(result)`.
#' 
## ------------------------------------------------------------------------
summary(result)
head(result)

#' 
#' Since [bigrquery](https://github.com/hadley/bigrquery) results are dataframes, we can make use of all sorts of other great R packages to do our downstream work.  Here we use a few more packages from the Hadleyverse: dplyr for data filtering and ggplot2 for visualization.
#' 
## ------------------------------------------------------------------------
# Change this filter if you want to visualize the result of this analysis for a different chromosome.
chromosomeOneResults <- filter(result, reference_name == "chr1" | reference_name == "1")

#' 
## ----message=FALSE-------------------------------------------------------
library(scales)
library(ggplot2)

#' 
## ----titv, fig.align="center", fig.width=10, message=FALSE, warning=FALSE, comment=NA----
ggplot(chromosomeOneResults, aes(x=window_start, y=titv)) +
  geom_point() +
  stat_smooth() +
  scale_x_continuous(labels=comma) +
  xlab("Genomic Position") +
  ylab("Ti/Tv") +
  ggtitle("Ti/Tv by 100,000 base pair windows on Chromosome 1")

#' 
#' ## Provenance
## ----provenance, comment=NA----------------------------------------------
sessionInfo()

