#' These vignettes require an oauth2 flow and therefore cannot be built as
#' part of the normal package build mechanism.
#'
#' Instead we build them manually and check in the rendered results.
knitAllVignettes <- function() {
  if(FALSE == grepl('doc$', getwd())) {
    stop("be sure to setwd('PATH/TO/inst/doc') before running this command.")
  }
  lapply(c("BigQueryIntroduction.Rmd",
           "ExpressionandMethylationCorrelation.Rmd",
           "ExpressionandProteinCorrelation.Rmd"
           ), function(rmd) {
    purl(rmd, documentation=2)
    knit(rmd)
    knit2html(rmd)
  }
  )
}
