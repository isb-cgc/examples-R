# These vignettes require an oauth2 flow and therefore cannot be built as
# part of the normal package build mechanism.
#
# Instead we use this internal helper function tobuild them manually and
# check in the rendered results to inst/doc in the package.
knitAllVignettes <- function() {
  library(knitr)
  opts_chunk$set(error=FALSE)
  if(FALSE == grepl('doc$', getwd())) {
    stop("be sure to setwd('PATH/TO/inst/doc') before running this command.")
  }
  lapply(c("BigQueryIntroduction.Rmd",
           "BCGSC_microRNA_expression.Rmd",
           "Copy_Number_segments.Rmd",
           "creating_cohort_gene_expression_matrices.Rmd",
           "Creating_TCGA_cohorts_part_1.Rmd",
           "Creating_TCGA_cohorts_part_2.Rmd",
           "DESeq2_tutorial.Rmd",
           "DNA_Methylation.Rmd",
           "ExpressionandMethylationCorrelation.Rmd",
           "ExpressionandProteinCorrelation.Rmd",
           "GenomicAndExpression_T_test.Rmd",
           "Protein_expression.Rmd",
           "Somatic_Mutations.Rmd",
           "TCGA_Annotations.Rmd",
           "UNC_HiSeq_mRNAseq_gene_expression_RSEM.Rmd",
           "Working_With_Barcode_Lists.Rmd"
           ), function(rmd) {
    knitr::purl(rmd, documentation=2)
    knitr::knit(rmd)
    knitr::knit2html(rmd)
  }
  )
}
# last comment.


knitOneVignette <- function(rmd) {
    # knit one .rmd file.
    library(knitr)
    opts_chunk$set(error=FALSE)
    if(FALSE == grepl('doc$', getwd())) {
        stop("be sure to setwd('PATH/TO/inst/doc') before running this command.")
    }
    knitr::purl(rmd, documentation=2)
    knitr::knit(rmd)
    knitr::knit2html(rmd)
}
