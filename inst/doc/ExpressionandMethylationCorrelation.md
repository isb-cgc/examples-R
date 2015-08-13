# Expression and Methylation Correlation

TODO introduction, talk about the tables we will use, etc...

Reproduces some of the analyses in http://watson.nci.nih.gov/~sdavis/tutorials/TCGA_data_integration/


```r
require(ggplot2)
require(ISBCGCExamples)

# The directory in which the files containing SQL reside.
#sqlDir <- file.path("/Users/deflaux/deflaux/examples-R/inst/", 
sqlDir <- file.path(system.file(package = "ISBCGCExamples"),
                    "sql")
```


```r
######################[ TIP ]########################################
## Set the Google Cloud Platform project id under which these queries will run.
## 
## If you are using the workshop docker image, this is already
## set for you in your .Rprofile and you can skip this step.

# project <- "YOUR-PROJECT-ID"
#####################################################################
```


```r
result <- DisplayAndDispatchQuery(file.path(sqlDir, "expression-methylation-correlation.sql"),
                                  project=project,
                                  replacements=list("_TUMOR_"="CESC",
                                                    "_CHR_"="chr9"))
```

```
# Compute the correlation between expression and methylation data on a particular chromosome for a particular tumor type.
SELECT HGNC_gene_symbol, Probe_ID, CORR(normalized_count, Beta_Value) AS correlation,
FROM (
  # we select the sample-barcode, gene-symbol, gene-expression, probe-id, and beta-value
  # from a "JOIN" of the gene expression data and the methylation data
  SELECT expr.SampleBarcode, HGNC_gene_symbol, normalized_count, Probe_ID, Beta_value
  FROM [isb-cgc:tcga_data_open.mRNA_UNC_HiSeq_RSEM] AS expr
  JOIN EACH (
    FLATTEN ( (
      # we select the sample-barcode, sample-type, study-name, probe-id, beta-value, and gene-symbol
      # from a "JOIN" of the methylation data and the methylation annotation tables
      # which are joined on the CpG probe id that exists in both tables
      # ( for speed we are only working with chr9 for now )
      SELECT SampleBarcode, SampleTypeLetterCode, Study, Probe_ID, Beta_Value, UCSC.RefGene_Name
      FROM [isb-cgc:tcga_data_open.Methylation_chr9] AS methData
      JOIN EACH [isb-cgc:platform_reference.methylation_annotation] AS methAnnot
      ON methData.Probe_ID = methAnnot.Name
      # we require that the gene-symbol not be null, and we are choosing only CESC TP samples for now
      WHERE ( UCSC.RefGene_Name IS NOT null
          AND SampleTypeLetterCode = 'TP'
          AND Study = 'CESC' )
    ), UCSC.RefGene_Name ) ) AS methyl
  ON
    methyl.UCSC.RefGene_Name = expr.HGNC_gene_symbol AND methyl.SampleBarcode = expr.SampleBarcode )
GROUP BY HGNC_gene_symbol, Probe_ID
ORDER BY correlation DESC
```
Number of rows returned by this query: 6127.

The result is one correlation value per row of data, each of which corresponds to a methylation probe and
its associated expression probe. Note that each expression probe may map to several methylation probes.

```r
head(result)
```

```
##   HGNC_gene_symbol   Probe_ID correlation
## 1            GPSM1 cg04305913   0.6015447
## 2            GPSM1 cg14934821   0.5881440
## 3           TMEM8B cg14087413   0.5672279
## 4         ADAMTS13 cg14206140   0.5471293
## 5          PIP5K1B cg13867370   0.5430617
## 6            GPSM1 cg01393841   0.5422848
```


```r
# Histogram overlaid with kernel density curve
ggplot(result, aes(x=correlation)) + 
    geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   binwidth=.05,
                   colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
```

<img src="figure/density-1.png" title="plot of chunk density" alt="plot of chunk density" style="display: block; margin: auto;" />

## Provenance

```r
sessionInfo()
```

```
R version 3.2.0 (2015-04-16)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.10.4 (Yosemite)

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] mgcv_1.8-6         nlme_3.1-120       ISBCGCExamples_0.1
[4] devtools_1.8.0     knitr_1.10.5       ggplot2_1.0.1     
[7] bigrquery_0.1.0   

loaded via a namespace (and not attached):
 [1] Rcpp_0.11.6      formatR_1.2      git2r_0.10.1     plyr_1.8.2      
 [5] tools_3.2.0      digest_0.6.8     jsonlite_0.9.16  evaluate_0.7    
 [9] memoise_0.2.1    gtable_0.1.2     lattice_0.20-31  Matrix_1.2-0    
[13] DBI_0.3.1        rstudioapi_0.3.1 curl_0.9.2       parallel_3.2.0  
[17] proto_0.3-10     dplyr_0.4.1      httr_1.0.0       stringr_1.0.0   
[21] xml2_0.1.1       rversions_1.0.1  grid_3.2.0       R6_2.1.0        
[25] reshape2_1.4.1   magrittr_1.5     scales_0.2.4     MASS_7.3-40     
[29] assertthat_0.1   mime_0.3         colorspace_1.2-6 labeling_0.3    
[33] stringi_0.5-5    lazyeval_0.1.10  munsell_0.4.2    markdown_0.7.7  
```
