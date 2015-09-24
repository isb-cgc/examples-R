# Expression and Methylation Correlation

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

## Pearson Correlation in BigQuery


```r
# Set the desired tables to query.
expressionTable = "isb-cgc:tcga_data_open.mRNA_UNC_HiSeq_RSEM"
methylationTable = "isb-cgc:tcga_data_open.Methylation_chr9"
# Add any additional clauses to be applied in WHERE to limit the methylation data further.
andWhere = "AND SampleTypeLetterCode = 'TP' AND Study = 'CESC'"

# Now we are ready to run the query.
result <- DisplayAndDispatchQuery(file.path(sqlDir, "expression-methylation-correlation.sql"),
                                  project=project,
                                  replacements=list("_EXPRESSION_TABLE_"=expressionTable,
                                                    "_METHYLATION_TABLE_"=methylationTable,
                                                    "#_AND_WHERE_"=andWhere))
```

```
# Compute the correlation between expression and methylation data.
SELECT HGNC_gene_symbol, Probe_ID, CORR(normalized_count, Beta_Value) AS correlation,
FROM (
  # We select the sample-barcode, gene-symbol, gene-expression, probe-id, and beta-value
  # from a "JOIN" of the gene expression data and the methylation data.
  SELECT expr.SampleBarcode, HGNC_gene_symbol, normalized_count, Probe_ID, Beta_value
  FROM [isb-cgc:tcga_data_open.mRNA_UNC_HiSeq_RSEM] AS expr
  JOIN EACH (
    FLATTEN ( (
      # We select the sample-barcode, sample-type, study-name, probe-id, beta-value, and gene-symbol
      # from a "JOIN" of the methylation data and the methylation annotation tables
      # which are joined on the CpG probe id that exists in both tables.
      # ( for speed we are only working with chr9 for now )
      SELECT SampleBarcode, SampleTypeLetterCode, Study, Probe_ID, Beta_Value, UCSC.RefGene_Name
      FROM [isb-cgc:tcga_data_open.Methylation_chr9] AS methData
      JOIN EACH [isb-cgc:platform_reference.methylation_annotation] AS methAnnot
      ON methData.Probe_ID = methAnnot.Name
      # We require that the gene-symbol not be null.
      WHERE
        UCSC.RefGene_Name IS NOT null
        # Optionally add clause here to limit the query to a particular
        # sample types and/or studies.
        AND SampleTypeLetterCode = 'TP' AND Study = 'CESC'
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

## Pearson Correlation in R

Now let's reproduce one of the results directly in R.

### Retrieve Expression Data

First we retrieve the expression data for a particular gene.

```r
# Set the desired gene to query.
gene = "GPSM1"

expressionData <- DisplayAndDispatchQuery(file.path(sqlDir, "expression-data.sql"),
                                  project=project,
                                  replacements=list("_EXPRESSION_TABLE_"=expressionTable,
                                                    "_GENE_"=gene))
```

```
# Retrieve expression data for a particular gene.
SELECT
  SampleBarcode,
  HGNC_gene_symbol,
  normalized_count
FROM [isb-cgc:tcga_data_open.mRNA_UNC_HiSeq_RSEM]
WHERE
  HGNC_gene_symbol = 'GPSM1'
```
Number of rows returned by this query: 10252.


```r
head(expressionData)
```

```
##      SampleBarcode HGNC_gene_symbol normalized_count
## 1 TCGA-GM-A2DD-01A            GPSM1         224.8186
## 2 TCGA-SQ-A6I4-11A            GPSM1         472.6918
## 3 TCGA-SQ-A6I4-01A            GPSM1         661.3882
## 4 TCGA-DD-A1EL-01A            GPSM1          54.2662
## 5 TCGA-DD-A1EL-11A            GPSM1          78.4641
## 6 TCGA-IW-A3M5-01A            GPSM1        2611.8053
```

### Retrieve Methylation Data

Then we retrieve the methylation data for a particular probe.


```r
# Set the desired probe to query.
probe = "cg04305913"

# Be sure to apply the same additional clauses to the WHERE to limit the methylation data further.

methylationData <- DisplayAndDispatchQuery(file.path(sqlDir, "methylation-data.sql"),
                                  project=project,
                                  replacements=list("_METHYLATION_TABLE_"=methylationTable,
                                                    "#_AND_WHERE_"=andWhere,
                                                    "_PROBE_"=probe))
```

```
# Retrieve methylation data for a particular probe id.
SELECT
  SampleBarcode,
  SampleTypeLetterCode,
  Study,
  Probe_ID,
  Beta_Value
FROM [isb-cgc:tcga_data_open.Methylation_chr9]
WHERE
  Probe_ID = 'cg04305913'
  # Optionally add clause here to limit the query to a particular
  # sample types and/or studies.
  AND SampleTypeLetterCode = 'TP' AND Study = 'CESC'
```

Number of rows returned by this query: 303.


```r
head(methylationData)
```

```
##      SampleBarcode SampleTypeLetterCode Study   Probe_ID Beta_Value
## 1 TCGA-C5-A1MQ-01A                   TP  CESC cg04305913       0.83
## 2 TCGA-MU-A5YI-01A                   TP  CESC cg04305913       0.50
## 3 TCGA-VS-A9UJ-01A                   TP  CESC cg04305913       0.79
## 4 TCGA-JW-A852-01A                   TP  CESC cg04305913       0.52
## 5 TCGA-FU-A3HY-01A                   TP  CESC cg04305913       0.31
## 6 TCGA-C5-A7UI-01A                   TP  CESC cg04305913       0.45
```

### Perform the correlation


```r
# TODO: fix package imports so that we can require(dplyr) and drop the package prefix on these calls.
data = dplyr::inner_join(expressionData, methylationData)
```

```
## Joining by: "SampleBarcode"
```

```r
head(data)
```

```
##      SampleBarcode HGNC_gene_symbol normalized_count SampleTypeLetterCode
## 1 TCGA-VS-A9UH-01A            GPSM1         419.5329                   TP
## 2 TCGA-IR-A3LA-01A            GPSM1         437.9900                   TP
## 3 TCGA-FU-A23K-01A            GPSM1         224.8253                   TP
## 4 TCGA-WL-A834-01A            GPSM1          55.8795                   TP
## 5 TCGA-VS-A9UU-01A            GPSM1         271.3057                   TP
## 6 TCGA-UC-A7PF-01A            GPSM1         251.5625                   TP
##   Study   Probe_ID Beta_Value
## 1  CESC cg04305913       0.76
## 2  CESC cg04305913       0.85
## 3  CESC cg04305913       0.50
## 4  CESC cg04305913       0.30
## 5  CESC cg04305913       0.49
## 6  CESC cg04305913       0.61
```


```r
cor(x=data$normalized_count, y=data$Beta_Value, method="pearson")
```

```
## [1] 0.6015447
```

And we can see that we have reproduced one of our results from BigQuery.

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
[21] curl_0.9.3       evaluate_0.7.2   labeling_0.3     stringi_0.5-5   
[25] scales_0.2.4     jsonlite_0.9.17  markdown_0.7.7   proto_0.3-10    
```
