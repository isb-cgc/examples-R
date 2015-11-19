# examples-R
Analysis examples based on the ISB-CGC hosted TCGA data, using R and R Markdown.

To install:
```
require(devtools) || install.packages("devtools")

# For now
install_local("PATH/TO/YOUR/CLONE/OF/examples-R", build_vignettes=TRUE)

# When public
install_github("isb-cgc", "examples-R", build_vignettes=TRUE)
```

To view and run the vignettes.
```
  help(package="ISBCGCExamples")
```

There are three main vignettes using TCGA data, with more in the works. 
They involve analyzing genomic data, correlating gene expression and methylation, 
and correlating protein and mRNA levels. 

The vignettes as R-markdown can be found in the examples-R/inst/doc directory,
which can serve as examples of using builtin BigQuery functions like Pearson
correlation, or even how to implement more complex functions like Spearmans 
correlation. Queries can be simple character vectors, or standalone files. 
Results are returned as data.frames using the bigrquery package to 
interact with the servers.

The SQL files used in the vignettes can be found at examples-R/inst/sql. 
These are parsed and dispatched with arguments using the DisplayAndDispatchQuery function, 
found in the file of the same name in examples-R/R.

If you have trouble with the OAuth, see examples-R/inst/doc/BigQueryIntroduction.html 
for some instructions on resetting it.

