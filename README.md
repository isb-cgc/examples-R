# examples-R

Analysis examples based on the ISB-CGC hosted TCGA data, using R and R Markdown.

To install:
```
library("devtools")
install_github("isb-cgc/examples-R", build_vignettes=TRUE)
```

To view and run the vignettes.
```
  library(ISBCGCExamples)
  help(package="ISBCGCExamples")
```

### Alpha tables are no longer available!

Please move to the "tcga_201607_beta" dataset, or even better, the newest
GDC datasets "TCGA_hg19_data_v0", "TCGA_hg38_data_v0", and "TCGA_bioclin_v0".

Some of these examples are using the alpha dataset that is now unavailble. If you see
a dataset that begins with "tcga_201510_alpha", then try "tcga_201607_beta", and
it's likely to work. We will be updating these over time.

If you want to move to the newest datasets (recommended), be aware that some of
the most common column names have changed to match the GDC's schemas. For example,
"Study" is now "project_short_name". "ParticipantBarcode" is not "case_barcode".
"SampleBarcode" is now "sample_barcode". Overall the column names have become
all lower case. Please get in touch if you're having trouble.

### OAuth

If you are having trouble with the **OAuth**, see the OAuth section below!

### Authentication on a remote server

To authenticate on a remote server, you need to use out-of-band authentication (OOB). The httr package has an option for this. Set "options(httr_oob_default=TRUE)" after loading bigrquery, but before calling query_exec(), and you should be good to go. [citation: https://support.rstudio.com/hc/en-us/articles/217952868-Generating-OAuth-tokens-from-a-server]

### Vignettes

There are vignettes for each TCGA data type, and more elaborate examples
involving analyzing genomic data, correlating gene expression and methylation,
and correlating protein and mRNA levels.

The vignettes as **R-markdown** can be found in the [examples-R/inst/doc](inst/doc) directory,
which can serve as examples of using builtin BigQuery functions like Pearson
correlation, or even how to implement more complex functions like Spearmans
correlation. Queries can be simple character vectors, or standalone files.
Results are returned as data.frames using the bigrquery package to
interact with the servers.

The **SQL** files used in the vignettes can be found at [examples-R/inst/sql](inst/sql).
These are parsed and dispatched with arguments using the DisplayAndDispatchQuery function,
found in the file of the same name in [examples-R/R](R).

## Intro to the CGC

[Big Query Introduction](inst/doc/BigQueryIntroduction.md)

[TCGA Annotations](inst/doc/TCGA_Annotations.md)

[Creating TCGA cohorts part 1](inst/doc/Creating_TCGA_cohorts_part_1.md)

[Creating TCGA cohorts part 2](inst/doc/Creating_TCGA_cohorts_part_2.md)

[Using the API endpoints to work with barcode lists](inst/doc/Working_With_Barcode_Lists.md)

[Constructing small matrices](inst/doc/creating_cohort_gene_expression_matrices.md)

## Available data types

[microRNA expression](inst/doc/BCGSC_microRNA_expression.md)

[Copy Number segments](inst/doc/Copy_Number_segments.md)

[DNA Methylation](inst/doc/DNA_Methylation.md)

[Protein expression](inst/doc/Protein_expression.md)

[Somatic Mutations](inst/doc/Somatic_Mutations.md)

[mRNAseq gene expression](inst/doc/UNC_HiSeq_mRNAseq_gene_expression_RSEM.md)

## Advanced examples

[DESeq2 workflow on raw data](inst/doc/DESeq2_tutorial.md)

[Expression and Copy Number Correlation](inst/doc/ExpressionAndCopyNumberCorrelation.md)

[Expression and Methylation Correlation](inst/doc/ExpressionandMethylationCorrelation.md)

[Expression and Protein Correlation](inst/doc/ExpressionandProteinCorrelation.md)

[Genomic And Expression T-test](inst/doc/GenomicAndExpression_T_test.md)

## Using Docker

[Processing Raw Data with Bioconductor](inst/doc/Processing_Raw_Data_With_Bioconductor.md)

[Bioconductor](http://www.bioconductor.org/) provides an excellent set of docker containers which include R, RStudio Server, and the sets of Bioconductor packages appropriate for certain use cases.

This R package is also available in a Docker container derived from `bioconductor/release_core`:
```
  b.gcr.io/isb-cgc-public-docker-images/r-examples
```
It can be run like so:
```
  docker run -p 8787:8787 -v YOUR_LOCAL_DIRECTORY:/home/rstudio/data \
    b.gcr.io/isb-cgc-public-docker-images/r-examples:latest
```
and then navigate to http://localhost:8787 on your local machine.

For more details, see [examples-R/inst/docker](inst/docker) and http://www.bioconductor.org/help/docker/.

Then log into Rstudio with username and password 'rstudio', for more details:
https://github.com/rocker-org/rocker/wiki/Using-the-RStudio-image

## OAuth

If you have trouble with the **OAuth**, see [examples-R/inst/doc/BigQueryIntroduction.html](inst/doc/BigQueryIntroduction.md)
for some instructions on resetting it.

## Important note about bigrquery and httr

There was an incompatibility between bigrquery and the httr library. If you are having trouble, try installing the development version of bigrquery or use the prior version of httr (1.0.0).

To install the dev version of bigrquery:
```
   https://github.com/rstats-db/bigrquery
   install.packages('devtools')
   devtools::install_github("rstats-db/bigrquery")
```
