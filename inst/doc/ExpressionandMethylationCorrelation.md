# Expression and Methylation Correlation

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

## Pearson Correlation in BigQuery


```r
# Set the desired tables to query.
expressionTable = "isb-cgc:tcga_data_open.mRNA_UNC_HiSeq_RSEM"
methylationTable = "isb-cgc:tcga_data_open.Methylation_chr9"
# Add any additional clauses to be applied in WHERE to limit the methylation data further.
andWhere = "AND SampleTypeLetterCode = 'TP' AND Study = 'CESC'"
# Do not correlate unless there are at least this many observations available
minimumNumberOfObservations = 30

# Now we are ready to run the query.
result = DisplayAndDispatchQuery(file.path(sqlDir, "expression-methylation-correlation.sql"),
                                 project=project,
                                 replacements=list("_EXPRESSION_TABLE_"=expressionTable,
                                                   "_METHYLATION_TABLE_"=methylationTable,
                                                   "#_AND_WHERE_"=andWhere,
                                                   "_MINIMUM_NUMBER_OF_OBSERVATIONS_"=minimumNumberOfObservations))
```

```
  # Compute the correlation between expression and methylation data.
SELECT
  HGNC_gene_symbol,
  Probe_ID,
  COUNT(1) AS num_observations,
  CORR(normalized_count, Beta_Value) AS correlation,
FROM (
    # We select the sample-barcode, gene-symbol, gene-expression, probe-id, and beta-value
    # from a "JOIN" of the gene expression data and the methylation data.
  SELECT
    expr.SampleBarcode,
    HGNC_gene_symbol,
    normalized_count,
    Probe_ID,
    Beta_value
  FROM
    [isb-cgc:tcga_data_open.mRNA_UNC_HiSeq_RSEM] AS expr
  JOIN EACH ( FLATTEN ( (
          # We select the sample-barcode, sample-type, study-name, probe-id, beta-value, and gene-symbol
          # from a "JOIN" of the methylation data and the methylation annotation tables
          # which are joined on the CpG probe id that exists in both tables.
          # ( for speed we are only working with chr9 for now )
        SELECT
          SampleBarcode,
          SampleTypeLetterCode,
          Study,
          Probe_ID,
          Beta_Value,
          UCSC.RefGene_Name
        FROM
          [isb-cgc:tcga_data_open.Methylation_chr9] AS methData
        JOIN EACH [isb-cgc:platform_reference.methylation_annotation] AS methAnnot
        ON
          methData.Probe_ID = methAnnot.Name
          # We require that the gene-symbol not be null.
        WHERE
          UCSC.RefGene_Name IS NOT NULL
          # Optionally add clause here to limit the query to a particular
          # sample types and/or studies.
          AND SampleTypeLetterCode = 'TP' AND Study = 'CESC'
          ), UCSC.RefGene_Name ) ) AS methyl
  ON
    methyl.UCSC.RefGene_Name = expr.HGNC_gene_symbol
    AND methyl.SampleBarcode = expr.SampleBarcode )
GROUP BY
  HGNC_gene_symbol,
  Probe_ID
HAVING
  num_observations >= 30
ORDER BY
  correlation DESC
```
Number of rows returned by this query: 6127.

The result is one correlation value per row of data, each of which corresponds to a methylation probe and
its associated expression probe. Note that each expression probe may map to several methylation probes.

```r
head(result)
```

```
##   HGNC_gene_symbol   Probe_ID num_observations correlation
## 1            GPSM1 cg04305913              602   0.6015447
## 2            GPSM1 cg14934821              602   0.5881440
## 3           TMEM8B cg14087413              602   0.5672279
## 4         ADAMTS13 cg14206140             1204   0.5471293
## 5          PIP5K1B cg13867370              301   0.5430617
## 6            GPSM1 cg01393841              602   0.5422848
```


```r
library(bigrquery)

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

expressionData = DisplayAndDispatchQuery(file.path(sqlDir, "expression-data.sql"),
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
ORDER BY
  SampleBarcode
```
Number of rows returned by this query: 10252.


```r
head(expressionData)
```

```
##      SampleBarcode HGNC_gene_symbol normalized_count
## 1 TCGA-02-0047-01A            GPSM1        1217.0011
## 2 TCGA-02-0055-01A            GPSM1         719.2712
## 3 TCGA-02-2483-01A            GPSM1        4626.0686
## 4 TCGA-02-2485-01A            GPSM1        1257.1429
## 5 TCGA-02-2486-01A            GPSM1         392.6829
## 6 TCGA-04-1348-01A            GPSM1          98.4136
```

### Retrieve Methylation Data

Then we retrieve the methylation data for a particular probe.


```r
# Set the desired probe to query.
probe = "cg04305913"

# Be sure to apply the same additional clauses to the WHERE to limit the methylation data further.

methylationData = DisplayAndDispatchQuery(file.path(sqlDir, "methylation-data.sql"),
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
ORDER BY
  SampleBarcode,
  SampleTypeLetterCode,
  Study
```

Number of rows returned by this query: 303.


```r
head(methylationData)
```

```
##      SampleBarcode SampleTypeLetterCode Study   Probe_ID Beta_Value
## 1 TCGA-2W-A8YY-01A                   TP  CESC cg04305913       0.76
## 2 TCGA-4J-AA1J-01A                   TP  CESC cg04305913       0.26
## 3 TCGA-BI-A0VR-01A                   TP  CESC cg04305913       0.50
## 4 TCGA-BI-A0VS-01A                   TP  CESC cg04305913       0.40
## 5 TCGA-BI-A20A-01A                   TP  CESC cg04305913       0.75
## 6 TCGA-C5-A0TN-01A                   TP  CESC cg04305913       0.81
```

### Perform the correlation


```r
library(dplyr)

data = inner_join(expressionData, methylationData)
```

```
## Joining by: "SampleBarcode"
```

```r
head(data)
```

```
##      SampleBarcode HGNC_gene_symbol normalized_count SampleTypeLetterCode
## 1 TCGA-2W-A8YY-01A            GPSM1         277.5393                   TP
## 2 TCGA-4J-AA1J-01A            GPSM1         119.6070                   TP
## 3 TCGA-BI-A0VR-01A            GPSM1          86.0595                   TP
## 4 TCGA-BI-A0VS-01A            GPSM1          91.6372                   TP
## 5 TCGA-BI-A20A-01A            GPSM1         256.8048                   TP
## 6 TCGA-C5-A0TN-01A            GPSM1         439.5386                   TP
##   Study   Probe_ID Beta_Value
## 1  CESC cg04305913       0.76
## 2  CESC cg04305913       0.26
## 3  CESC cg04305913       0.50
## 4  CESC cg04305913       0.40
## 5  CESC cg04305913       0.75
## 6  CESC cg04305913       0.81
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
[1] mgcv_1.8-6         nlme_3.1-120       ggplot2_1.0.1     
[4] scales_0.2.5       bigrquery_0.1.0    dplyr_0.4.2       
[7] ISBCGCExamples_0.1

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.0      rstudioapi_0.3.1 knitr_1.10.5     magrittr_1.5    
 [5] MASS_7.3-40      munsell_0.4.2    lattice_0.20-31  colorspace_1.2-6
 [9] R6_2.1.1         stringr_1.0.0    httr_1.0.0       plyr_1.8.3      
[13] tools_3.2.0      parallel_3.2.0   grid_3.2.0       gtable_0.1.2    
[17] DBI_0.3.1        lazyeval_0.1.10  assertthat_0.1   digest_0.6.8    
[21] Matrix_1.2-0     reshape2_1.4.1   formatR_1.2      curl_0.9.3      
[25] mime_0.4         evaluate_0.7.2   labeling_0.3     stringi_0.5-5   
[29] jsonlite_0.9.17  markdown_0.7.7   proto_0.3-10    
```
