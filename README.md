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
