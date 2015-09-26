# These vignettes require an oauth2 flow and therefore cannot be built as
# part of the normal package build mechanism.
#
# Instead we use this internal helper function tobuild them manually and
# check in the rendered results to inst/doc in the package.
knitAllVignettes <- function() {
  if(FALSE == grepl('doc$', getwd())) {
    stop("be sure to setwd('PATH/TO/inst/doc') before running this command.")
  }
  lapply(c("BigQueryIntroduction.Rmd",
           "ExpressionandMethylationCorrelation.Rmd",
           "ExpressionandProteinCorrelation.Rmd"
           ), function(rmd) {
    knitr::purl(rmd, documentation=2)
    knitr::knit(rmd)
    knitr::knit2html(rmd)
  }
  )
}
