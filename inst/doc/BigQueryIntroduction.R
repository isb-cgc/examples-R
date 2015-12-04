#' # Exploring the TCGA data in BigQuery
#' 
#' The ISB-CGC (isb-cgc.org) project has aggregated and curated all of the TCGA open-access clinical, biospecimen, and Level-3 molecular data and uploaded it into BigQuery tables that are open to the public.  Here we will show you how you can begin to work with these tables from the familiar R environment.
#' 
#' ### Helpful BigQuery links
#' 
#' For this example, we'll also be working with [Google BigQuery](https://cloud.google.com/bigquery/). It's often helpful to have a [link to the docs](https://cloud.google.com/bigquery/what-is-bigquery) handy, and especially the [query reference](https://cloud.google.com/bigquery/query-reference).
#' 
#' ## Run a query from R
#' 
#' We will start by loading four R packages:
#' - the [bigrquery](https://github.com/hadley/bigrquery) package written by Hadley Wickham implements an R interface to [Google BigQuery](https://cloud.google.com/bigquery/),
#' - the [dplyr](https://github.com/hadley/dplyr) package provides a set of tools for efficiently manipulating datasets in R, and
#' - the [ggplot2](https://github.com/hadley/ggplot2) package for elegant graphics, and
#' - the [scales](https://github.com/hadley/scales) package for visualization-oriented scale functions.
#' 
## ----message=FALSE-------------------------------------------------------
require(dplyr) || install.packages("dplyr")
require(bigrquery) || install.packages("bigrquery")
require(scales) || install.packages("scales")
require(ggplot2) || install.packages("ggplot2")
library(ISBCGCExamples)

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
#' Let's start by working with one of the simplest tables, the Clinical_data table.  The format of a table name in BigQuery is <project_name>:<dataset_name>.<table_name>
#' 
#' 
## ------------------------------------------------------------------------
theTable <- "isb-cgc:tcga_201510_alpha.Clinical_data"

#' 
#' Note that when you send the first query, you will need to go through the authentication flow with BigQuery.  You will be provided with a url to cut and  paste into your browser, and then you will get an authorization code to cut and paste back here.
#' 
#' [bigrquery](https://github.com/hadley/bigrquery) uses the package [httr](https://github.com/hadley/httr) to perform OAuth.
#' 
#' Let's start by just counting the number of records in the table.
#' First we'll just create the query string and echo it:
## ------------------------------------------------------------------------
querySql <- paste("SELECT COUNT(1) FROM [", theTable, "]", sep="")
querySql

#' 
#' And then we'll send the query to the cloud for execution:
## ------------------------------------------------------------------------
result <- query_exec(querySql, project=project)
result

#' 
#' And we see that the table has 1 row - this is the number of unique patients or participants across all of the various TCGA studies.
#' 
## ----eval=FALSE----------------------------------------------------------
## ######################[ TIP ]########################################
## ## If you have any trouble with OAuth and need to redo/reset OAuth,
## ## run the following code.
## 
## #if (FALSE != getOption("httr_oauth_cache")) {
## #  file.remove(getOption("httr_oauth_cache"))
## #}
## 
## ## or maybe it should look like this?
## 
## #if (!is.null(getOption("httr_oauth_cache"))) {
## #  file.remove(getOption("httr_oauth_cache"))
## #}
## 
## #message("Restart R to redo/reset OAuth.")
## #####################################################################

#' 
#' 
#' ## Run a query using the BigQuery Web User Interface
#' 
#' So what is actually in this table?  Click on [this link](https://bigquery.cloud.google.com/table/isb-cgc:tcga_201507_alpha.Clinical_data) to view the schema in the BigQuery web user interface.
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
#' This allows queries to be more easily shared among analyses and also reused for different tables.  For example, in the following file we have a query that will count the number of patients, grouped by disease type in any of our TCGA data tables.
## ------------------------------------------------------------------------
file.show(file.path(system.file(package = "ISBCGCExamples"),
                    "sql",
                    "count-patients-by-study.sql"))

#' 
#' Now let's run the query to see what these counts are in the Clinical_data table:
## ----comment=NA----------------------------------------------------------
result <- DisplayAndDispatchQuery(file.path(system.file(package = "ISBCGCExamples"),
                                            "sql",
                                            "count-patients-by-study.sql"),
                                  project=project,
                                  replacements=list("_THE_TABLE_"=theTable))
cat("Number of rows returned by this query: ", nrow(result), "\n")

#' 
#' Results from [bigrquery](https://github.com/hadley/bigrquery) are returned as R dataframes, meaning that we can make use of all of the regular dataframe functions as well as all sorts of other great R packages to do our downstream work.
## ------------------------------------------------------------------------
mode(result)
class(result)
summary(result)
head(result)

#' 
#' ## Visualize Query Results
#' 
#' Since there are over 30 distinct tumor types within the TCGA project, we may want to filter our results before visualizing.  For example let's look only at tumor types with at least 500 patients in the study:
#' 
## ------------------------------------------------------------------------
subsetResults <- filter(result, n>=500)
subsetResults <- arrange(subsetResults,desc(n))

#' 
#' and then create a barchart of the patient counts:
#' 
## ----titv, fig.align="center", fig.width=10, message=FALSE, warning=FALSE, comment=NA----
ggplot(subsetResults, aes(x=Study, y=n, fill=Study)) +
  geom_bar(stat="identity") +
  ylab("Number of Patients") +
  ggtitle("Study Size")

#' 
#' ## Provenance
## ----provenance, comment=NA----------------------------------------------
sessionInfo()

