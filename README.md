# examples-R
Analysis examples based on the ISB-CGC hosted TCGA data, using R and R Markdown.

To install:
```
require(devtools) || install.packages("devtools")
install_github("isb-cgc", "examples-R", build_vignettes=TRUE)
```

To view and run the vignettes.
```
  help(package="ISBCGCExamples")
```

There are four main vignettes using TCGA data, with more in the works. 
They involve analyzing genomic data, correlating gene expression and methylation, 
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

If you have trouble with the **OAuth**, see [examples-R/inst/doc/BigQueryIntroduction.html](inst/doc/BigQueryIntroduction.md) 
for some instructions on resetting it.

## Docker

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
