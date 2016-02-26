# examples-R

Analysis examples based on the ISB-CGC hosted TCGA data, using R and R Markdown.

###NOTE: There is an incompatibility between bigrquery and the httr library. Until bigrquery is updated, please use the development branch of bigrquery or use the prior version of httr (1.0.0).

To install the dev version of bigrquery:
```
   https://github.com/rstats-db/bigrquery
   install.packages('devtools')
   devtools::install_github("rstats-db/bigrquery")
```

To install:
```
require(devtools) || install.packages("devtools")
install_github("isb-cgc/examples-R", build_vignettes=TRUE)
```

To view and run the vignettes.
```
  help(package="ISBCGCExamples")
```

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

## Available data types

[microRNA expression](inst/doc/BCGSC_microRNA_expression.md)

[Copy Number segments](inst/doc/Copy_Number_segments.md)

[DNA Methylation](inst/doc/DNA_Methylation.md)

[Protein expression](inst/doc/Protein_expression.md)

[Somatic Mutations](inst/doc/Somatic_Mutations.md)

[mRNAseq gene expression](inst/doc/UNC_HiSeq_mRNAseq_gene_expression_RSEM.md)

## Advanced examples

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
