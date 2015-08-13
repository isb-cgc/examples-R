#' # Expression and Methylation Correlation
#' 
#' TODO introduction, talk about the tables we will use, etc...
#' 
#' Reproduces some of the analyses in http://watson.nci.nih.gov/~sdavis/tutorials/TCGA_data_integration/
#' 
## ----message=FALSE-------------------------------------------------------
require(ggplot2)
require(ISBCGCExamples)

# The directory in which the files containing SQL reside.
#sqlDir <- file.path("/Users/deflaux/deflaux/examples-R/inst/", 
sqlDir <- file.path(system.file(package = "ISBCGCExamples"),
                    "sql")

#' 
## ----eval=FALSE----------------------------------------------------------
## ######################[ TIP ]########################################
## ## Set the Google Cloud Platform project id under which these queries will run.
## ##
## ## If you are using the workshop docker image, this is already
## ## set for you in your .Rprofile and you can skip this step.
## 
## # project <- "YOUR-PROJECT-ID"
## #####################################################################

#' 
## ----comment=NA----------------------------------------------------------
result <- DisplayAndDispatchQuery(file.path(sqlDir, "expression-methylation-correlation.sql"),
                                  project=project,
                                  replacements=list("_TUMOR_"="CESC",
                                                    "_CHR_"="chr9"))

#' Number of rows returned by this query: `r nrow(result)`.
#' 
#' The result is one correlation value per row of data, each of which corresponds to a methylation probe and
#' its associated expression probe. Note that each expression probe may map to several methylation probes.
## ------------------------------------------------------------------------
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
#' ## Provenance
## ----provenance, comment=NA----------------------------------------------
sessionInfo()

