# Modelled after https://github.com/Bioconductor/bioc_docker/blob/master/out/devel_sequencing/installFromBiocViews.R
library(BiocInstaller)

pkgs_to_install <- c("devtools")
github_pkgs_to_install <- c("isb-cgc/examples-R", "googlegenomics/bioconductor-workshop-r")

cores <- max(2, parallel::detectCores()-2)
if (parallel::detectCores() == 1)
    cores <- 1
options(list("Ncpus"=cores))

tryCatch({
  biocLite(pkgs_to_install)
  biocLite(github_pkgs_to_install, build_vignettes=TRUE, dependencies=TRUE)
},
warning=function(w){
    if(length(grep("is not available|had non-zero exit status|installation of one or more packages failed", w$message)))
        stop(sprintf("got a fatal warning: %s", w$message))
})

warnings()

if (!is.null(warnings()))
{
    w <- capture.output(warnings())
    if (length(grep(
     "is not available|had non-zero exit status|installation of one or more packages failed", w)))
        quit("no", 1L)
}
