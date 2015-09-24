# Expression and Protein Correlation

TODO introduction, talk about the tables we will use, etc...

Reproduces some of the analyses in http://watson.nci.nih.gov/~sdavis/tutorials/TCGA_data_integration/


```r
require(ISBCGCExamples)

# The directory in which the files containing SQL reside.
#sqlDir <- file.path("/PATH/TO/GIT/CLONE/OF/examples-R/inst/", 
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

## Spearman Correlation in BigQuery


```r
# Set the desired tables to query.
expressionTable = "isb-cgc:tcga_data_open.mRNA_UNC_HiSeq_RSEM"
proteinTable = "isb-cgc:tcga_data_open.Protein"
cohortTable = "isb-cgc:test.cohort_14jun2015"

# Now we are ready to run the query.
result <- DisplayAndDispatchQuery(file.path(sqlDir, "protein-mrna-spearman-correlation.sql"),
                                  project=project,
                                  replacements=list("_EXPRESSION_TABLE_"=expressionTable,
                                                    "_PROTEIN_TABLE_"=proteinTable,
                                                    "_COHORT_TABLE_"=cohortTable))
```

```
# Correlate the protein quantification data with the mRNA expression data.
SELECT
  feat1.gene AS gene,
  CORR(feat1.exp_rank, feat2.exp_rank) AS spearman_corr
FROM (
  SELECT
    *,
    RANK() OVER (PARTITION BY gene ORDER BY protein_expression ASC) AS exp_rank
  FROM (
    SELECT
      SampleBarcode,
      Gene_Name AS gene,
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
    RANK() OVER (PARTITION BY gene ORDER BY log2_normalized_count ASC) AS exp_rank
  FROM (
    SELECT
      SampleBarcode,
      HGNC_gene_symbol AS gene,
      IF(0 = normalized_count, 0, LOG2(normalized_count)) AS log2_normalized_count
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
  gene
ORDER BY
  spearman_corr
```
Number of rows returned by this query: 187.

The result is one correlation value per row of data, each of which corresponds to  . . . MORE HERE

```r
head(result)
```

```
##     gene spearman_corr
## 1 NFE2L2    -0.9770990
## 2   E2F1    -0.8063829
## 3    RET    -0.7865162
## 4   TP63    -0.6226925
## 5    SYP    -0.5138228
## 6   TFF1    -0.4943634
```


```r
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
gene = "NFE2L2"

expressionData <- DisplayAndDispatchQuery(file.path(sqlDir, "expression-data-by-cohort.sql"),
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
  HGNC_gene_symbol = 'NFE2L2'
  AND SampleBarcode IN (
  SELECT
    sample_barcode
  FROM
    [isb-cgc:test.cohort_14jun2015] )
```
Number of rows returned by this query: 472.


```r
head(expressionData)
```

```
##      SampleBarcode HGNC_gene_symbol normalized_count
## 1 TCGA-DK-A6AW-01A           NFE2L2        2359.1885
## 2 TCGA-AB-2810-03A           NFE2L2         926.3524
## 3 TCGA-EE-A2GR-06A           NFE2L2        2078.4668
## 4 TCGA-BR-A4CS-01A           NFE2L2        2021.3705
## 5 TCGA-EJ-7784-11A           NFE2L2        2923.7570
## 6 TCGA-CG-4469-01A           NFE2L2        3755.1911
```

### Retrieve Protein Data

Then we retrieve the protein data for a particular gene.


```r
proteinData <- DisplayAndDispatchQuery(file.path(sqlDir, "protein-data-by-cohort.sql"),
                                  project=project,
                                  replacements=list("_PROTEIN_TABLE_"=proteinTable,
                                                    "_COHORT_TABLE_"=cohortTable,
                                                    "_GENE_"=gene))
```

```
# Retrieve protein data for a particular gene for samples within a cohort.
SELECT
  SampleBarcode,
  Gene_Name,
  protein_expression
FROM
  [isb-cgc:tcga_data_open.Protein]
WHERE
  Gene_Name = 'NFE2L2'
  AND SampleBarcode IN (
  SELECT
    sample_barcode
  FROM
    [isb-cgc:test.cohort_14jun2015] )
```

Number of rows returned by this query: 3.


```r
head(proteinData)
```

```
##      SampleBarcode Gene_Name protein_expression
## 1 TCGA-CQ-6227-01A    NFE2L2          0.9564625
## 2 TCGA-55-6981-01A    NFE2L2          1.4299760
## 3 TCGA-75-6205-01A    NFE2L2          1.3022501
```

### Perform the correlation


```r
# TODO: fix package imports so that we can require(dplyr) and drop the package prefix on these calls.
data = dplyr::inner_join(expressionData, proteinData)
```

```
## Joining by: "SampleBarcode"
```

```r
head(data)
```

```
##      SampleBarcode HGNC_gene_symbol normalized_count Gene_Name
## 1 TCGA-75-6205-01A           NFE2L2         2580.644    NFE2L2
## 2 TCGA-55-6981-01A           NFE2L2         2130.798    NFE2L2
## 3 TCGA-CQ-6227-01A           NFE2L2         3692.627    NFE2L2
##   protein_expression
## 1          1.3022501
## 2          1.4299760
## 3          0.9564625
```


```r
# TODO: This isn't quite right, needs a review.
cor(x=data$normalized_count, y=data$protein_expression, method="spearman")
```

```
## [1] -1
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
[1] ISBCGCExamples_0.1 knitr_1.10.5       ggplot2_1.0.1     
[4] bigrquery_0.1.0   

loaded via a namespace (and not attached):
 [1] Rcpp_0.11.6      rstudioapi_0.3.1 magrittr_1.5     MASS_7.3-40     
 [5] munsell_0.4.2    colorspace_1.2-6 R6_2.1.0         stringr_1.0.0   
 [9] httr_1.0.0       plyr_1.8.2       dplyr_0.4.1      tools_3.2.0     
[13] parallel_3.2.0   grid_3.2.0       gtable_0.1.2     DBI_0.3.1       
[17] assertthat_0.1   digest_0.6.8     reshape2_1.4.1   formatR_1.2     
[21] curl_0.9.3       mime_0.4         evaluate_0.7.2   labeling_0.3    
[25] stringi_0.5-5    scales_0.2.4     jsonlite_0.9.17  markdown_0.7.7  
[29] proto_0.3-10    
```
