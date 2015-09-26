# Expression and Protein Correlation

TODO introduction, talk about the tables we will use, etc...

Reproduces some of the analyses in http://watson.nci.nih.gov/~sdavis/tutorials/TCGA_data_integration/


```r
library(ISBCGCExamples)

# The directory in which the files containing SQL reside.
#sqlDir = file.path("/PATH/TO/GIT/CLONE/OF/examples-R/inst/", 
sqlDir = file.path(system.file(package = "ISBCGCExamples"),
                    "sql")
```


```r
######################[ TIP ]########################################
## Set the Google Cloud Platform project id under which these queries will run.
## 
## If you are using the workshop docker image, this is already
## set for you in your .Rprofile and you can skip this step.

# project = "YOUR-PROJECT-ID"
#####################################################################
```

## Spearman Correlation in BigQuery


```r
# Set the desired tables to query.
expressionTable = "isb-cgc:tcga_data_open.mRNA_UNC_HiSeq_RSEM"
proteinTable = "isb-cgc:tcga_data_open.Protein"
cohortTable = "isb-cgc:test.cohort_14jun2015"

# Now we are ready to run the query.
result = DisplayAndDispatchQuery(file.path(sqlDir, "protein-mrna-spearman-correlation.sql"),
                                 project=project,
                                 replacements=list("_EXPRESSION_TABLE_"=expressionTable,
                                                   "_PROTEIN_TABLE_"=proteinTable,
                                                   "_COHORT_TABLE_"=cohortTable))
```

```
# Correlate the protein quantification data with the mRNA expression data.
SELECT
  feat1.gene AS gene,
  feat1.protein AS protein,
  CORR(feat1.exp_rank, feat2.exp_rank) AS spearman_corr
FROM (
  SELECT
    *,
    RANK() OVER (PARTITION BY protein ORDER BY protein_expression ASC) AS exp_rank
  FROM (
    SELECT
      SampleBarcode,
      Gene_Name AS gene,
      Protein_Name AS protein,
      protein_expression
    FROM
      [isb-cgc:tcga_data_open.Protein]
    WHERE
      SampleBarcode IN (
      SELECT
        sample_barcode
      FROM
        [isb-cgc:test.cohort_14jun2015] ) ) ) feat1
JOIN EACH (
  SELECT
    *,
    RANK() OVER (PARTITION BY gene ORDER BY log2_count ASC) AS exp_rank
  FROM (
    SELECT
      SampleBarcode,
      HGNC_gene_symbol AS gene,
      LOG2(normalized_count+1) AS log2_count
    FROM
      [isb-cgc:tcga_data_open.mRNA_UNC_HiSeq_RSEM]
    WHERE
      SampleBarcode IN (
      SELECT
        sample_barcode
      FROM
        [isb-cgc:test.cohort_14jun2015] ) ) )feat2
ON
  feat1.SampleBarcode = feat2.SampleBarcode
  AND feat1.gene = feat2.gene
GROUP BY
  gene,
  protein
ORDER BY
  spearman_corr
```
Number of rows returned by this query: 247.

The result is one correlation value per row of data, each of which corresponds to  . . . MORE HERE

```r
head(result)
```

```
##     gene     protein spearman_corr
## 1 NFE2L2        Nrf2    -0.9770990
## 2   E2F1        E2F1    -0.8063829
## 3    RET   Ret_pY905    -0.7865162
## 4   TFF1        TFF1    -0.6824969
## 5   TP63         p63    -0.6226925
## 6  ERBB3 HER3_pY1298    -0.5343018
```


```r
library(bigrquery)

# Histogram overlaid with kernel density curve
ggplot(result, aes(x=spearman_corr)) + 
    geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   binwidth=.05,
                   colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
```

<img src="figure/spearman_density-1.png" title="plot of chunk spearman_density" alt="plot of chunk spearman_density" style="display: block; margin: auto;" />

## Spearman Correlation in R

Now let's reproduce one of the results directly in R.

### Retrieve Expression Data

First we retrieve the expression data for a particular gene.

```r
# Set the desired gene to query.
gene = "ERBB3"

expressionData = DisplayAndDispatchQuery(file.path(sqlDir, "expression-data-by-cohort.sql"),
                                         project=project,
                                         replacements=list("_EXPRESSION_TABLE_"=expressionTable,
                                                           "_COHORT_TABLE_"=cohortTable,
                                                           "_GENE_"=gene))
```

```
# Retrieve expression data for a particular gene for samples within a cohort.
SELECT
  SampleBarcode,
  HGNC_gene_symbol,
  normalized_count
FROM [isb-cgc:tcga_data_open.mRNA_UNC_HiSeq_RSEM]
WHERE
  HGNC_gene_symbol = 'ERBB3'
  AND SampleBarcode IN (
  SELECT
    sample_barcode
  FROM
    [isb-cgc:test.cohort_14jun2015] )
ORDER BY
  SampleBarcode
```
Number of rows returned by this query: 472.


```r
head(expressionData)
```

```
##      SampleBarcode HGNC_gene_symbol normalized_count
## 1 TCGA-02-2485-01A            ERBB3          23.3333
## 2 TCGA-05-4244-01A            ERBB3        4069.8141
## 3 TCGA-05-4250-01A            ERBB3        3102.3792
## 4 TCGA-06-0138-01A            ERBB3         365.5840
## 5 TCGA-06-0157-01A            ERBB3          61.3897
## 6 TCGA-06-0158-01A            ERBB3         307.8278
```

### Retrieve Protein Data

Then we retrieve the protein data for a particular gene.


```r
protein = "HER3_pY1298"

proteinData = DisplayAndDispatchQuery(file.path(sqlDir, "protein-data-by-cohort.sql"),
                                      project=project,
                                      replacements=list("_PROTEIN_TABLE_"=proteinTable,
                                                        "_COHORT_TABLE_"=cohortTable,
                                                        "_GENE_"=gene,
                                                        "_PROTEIN_"=protein))
```

```
# Retrieve protein data for a particular gene for samples within a cohort.
SELECT
  SampleBarcode,
  Gene_Name,
  Protein_Name,
  protein_expression
FROM
  [isb-cgc:tcga_data_open.Protein]
WHERE
  Gene_Name = 'ERBB3'
  AND Protein_Name = 'HER3_pY1298'
  AND SampleBarcode IN (
  SELECT
    sample_barcode
  FROM
    [isb-cgc:test.cohort_14jun2015] )
ORDER BY
  SampleBarcode
```

Number of rows returned by this query: 88.


```r
head(proteinData)
```

```
##      SampleBarcode Gene_Name Protein_Name protein_expression
## 1 TCGA-02-0116-01A     ERBB3  HER3_pY1298          -1.380045
## 2 TCGA-02-2485-01A     ERBB3  HER3_pY1298          -1.481069
## 3 TCGA-06-1086-01A     ERBB3  HER3_pY1298          -1.903248
## 4 TCGA-06-2564-01A     ERBB3  HER3_pY1298          -1.245319
## 5 TCGA-06-2565-01A     ERBB3  HER3_pY1298          -1.863513
## 6 TCGA-06-5414-01A     ERBB3  HER3_pY1298          -1.718002
```

### Perform the correlation


```r
library(dplyr)

data = inner_join(expressionData, proteinData)
```

```
## Joining by: "SampleBarcode"
```

```r
head(data)
```

```
##      SampleBarcode HGNC_gene_symbol normalized_count Gene_Name
## 1 TCGA-02-2485-01A            ERBB3          23.3333     ERBB3
## 2 TCGA-06-2564-01A            ERBB3         366.0619     ERBB3
## 3 TCGA-06-2565-01A            ERBB3         461.3095     ERBB3
## 4 TCGA-06-5414-01A            ERBB3         256.5712     ERBB3
## 5 TCGA-12-3650-01A            ERBB3         244.0068     ERBB3
## 6 TCGA-12-3652-01A            ERBB3          57.8406     ERBB3
##   Protein_Name protein_expression
## 1  HER3_pY1298          -1.481069
## 2  HER3_pY1298          -1.245319
## 3  HER3_pY1298          -1.863513
## 4  HER3_pY1298          -1.718002
## 5  HER3_pY1298          -2.131116
## 6  HER3_pY1298          -1.703800
```


```r
# TODO: This isn't quite right, needs a review.
cor(x=data$normalized_count, y=data$protein_expression, method="spearman")
```

```
## [1] -0.4738548
```

## Provenance

```r
sessionInfo()
```

```
R version 3.2.0 (2015-04-16)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.10.5 (Yosemite)

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] mgcv_1.8-6         nlme_3.1-120       ggplot2_1.0.1     
[4] scales_0.2.5       ISBCGCExamples_0.1 bigrquery_0.1.0   
[7] dplyr_0.4.2       

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.0      rstudioapi_0.3.1 knitr_1.10.5     magrittr_1.5    
 [5] MASS_7.3-40      munsell_0.4.2    colorspace_1.2-6 lattice_0.20-31 
 [9] R6_2.1.1         stringr_1.0.0    httr_1.0.0       plyr_1.8.3      
[13] tools_3.2.0      parallel_3.2.0   grid_3.2.0       gtable_0.1.2    
[17] DBI_0.3.1        htmltools_0.2.6  lazyeval_0.1.10  assertthat_0.1  
[21] digest_0.6.8     Matrix_1.2-0     formatR_1.2      reshape2_1.4.1  
[25] curl_0.9.3       mime_0.4         evaluate_0.7.2   rmarkdown_0.7   
[29] labeling_0.3     stringi_0.5-5    jsonlite_0.9.17  markdown_0.7.7  
[33] proto_0.3-10    
```
