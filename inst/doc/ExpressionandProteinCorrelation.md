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
expressionTable = "isb-cgc:tcga_201507_alpha.mRNA_UNC_HiSeq_RSEM"
proteinTable = "isb-cgc:tcga_201507_alpha.Protein_RPPA_data"
cohortTable = "isb-cgc:tcga_cohorts.BRCA"

# Do not correlate unless there are at least this many observations available
minimumNumberOfObservations = 30

# Now we are ready to run the query.
result = DisplayAndDispatchQuery(file.path(sqlDir, "protein-mrna-spearman-correlation.sql"),
                                 project=project,
                                 replacements=list("_EXPRESSION_TABLE_"=expressionTable,
                                                   "_PROTEIN_TABLE_"=proteinTable,
                                                   "_COHORT_TABLE_"=cohortTable,
                                                   "_MINIMUM_NUMBER_OF_OBSERVATIONS_"=minimumNumberOfObservations))
```

```
# Correlate the protein quantification data with the mRNA expression data.
SELECT
  gene,
  protein,
  COUNT(1) AS num_observations,
  CORR(expr_rank, prot_rank) AS spearman_corr
FROM (
  SELECT
    barcode,
    gene,
    protein,
    RANK() OVER (PARTITION BY gene, protein ORDER BY log2_count ASC) AS expr_rank,
    RANK() OVER (PARTITION BY gene, protein ORDER BY protein_expression ASC) AS prot_rank,
  FROM (
    SELECT
      feat1.SampleBarcode AS barcode,
      Gene_Name AS gene,
      Protein_Name AS protein,
      protein_expression,
      log2_count,
    FROM (
      SELECT
        SampleBarcode,
        Gene_Name,
        Protein_Name,
        protein_expression
      FROM
        [isb-cgc:tcga_201507_alpha.Protein_RPPA_data]
      WHERE
        SampleBarcode IN (
        SELECT
          SampleBarcode
        FROM
          [isb-cgc:tcga_cohorts.BRCA] ) ) feat1
    JOIN EACH (
      SELECT
        SampleBarcode,
        HGNC_gene_symbol,
        LOG2(normalized_count+1) AS log2_count
      FROM
        [isb-cgc:tcga_201507_alpha.mRNA_UNC_HiSeq_RSEM]
      WHERE
        SampleBarcode IN (
        SELECT
          SampleBarcode
        FROM
          [isb-cgc:tcga_cohorts.BRCA] ) ) feat2
    ON
      feat1.SampleBarcode = feat2.SampleBarcode
      AND feat1.Gene_Name = feat2.HGNC_gene_symbol))
GROUP BY
  gene,
  protein
HAVING
  num_observations >= 30
ORDER BY
  spearman_corr DESC
```
Number of rows returned by this query: 133.

The result is one correlation value per row of data, each of which corresponds to  . . . MORE HERE

```r
head(result)
```

```
##     gene  protein num_observations spearman_corr
## 1   ESR1 ER-alpha              403     0.9075414
## 2    PGR       PR              403     0.8718714
## 3   BCL2    Bcl-2              403     0.8438356
## 4  GATA3    GATA3              403     0.8320966
## 5 IGFBP2   IGFBP2              403     0.7975304
## 6   ASNS     ASNS              403     0.7523398
```


```r
library(ggplot2)

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
gene = "ESR1"

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
FROM [isb-cgc:tcga_201507_alpha.mRNA_UNC_HiSeq_RSEM]
WHERE
  HGNC_gene_symbol = 'ESR1'
  AND SampleBarcode IN (
  SELECT
    SampleBarcode
  FROM
    [isb-cgc:tcga_cohorts.BRCA] )
ORDER BY
  SampleBarcode
```
Number of rows returned by this query: 1202.


```r
head(expressionData)
```

```
##      SampleBarcode HGNC_gene_symbol normalized_count
## 1 TCGA-3C-AAAU-01A             ESR1        3457.9620
## 2 TCGA-3C-AALI-01A             ESR1          68.5155
## 3 TCGA-3C-AALJ-01A             ESR1        7482.3209
## 4 TCGA-3C-AALK-01A             ESR1        2485.3124
## 5 TCGA-4H-AAAK-01A             ESR1        5518.2979
## 6 TCGA-5L-AAT0-01A             ESR1        6592.1689
```

### Retrieve Protein Data

Then we retrieve the protein data for a particular gene.


```r
protein = "ER-alpha"

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
  [isb-cgc:tcga_201507_alpha.Protein_RPPA_data]
WHERE
  Gene_Name = 'ESR1'
  AND Protein_Name = 'ER-alpha'
  AND SampleBarcode IN (
  SELECT
    SampleBarcode
  FROM
    [isb-cgc:tcga_cohorts.BRCA] )
ORDER BY
  SampleBarcode
```
Number of rows returned by this query: 404.


```r
head(proteinData)
```

```
##      SampleBarcode Gene_Name Protein_Name protein_expression
## 1 TCGA-A1-A0SH-01A      ESR1     ER-alpha         -1.6136630
## 2 TCGA-A1-A0SJ-01A      ESR1     ER-alpha          0.2359291
## 3 TCGA-A1-A0SK-01A      ESR1     ER-alpha         -3.7201419
## 4 TCGA-A1-A0SO-01A      ESR1     ER-alpha         -4.0307670
## 5 TCGA-A2-A04N-01A      ESR1     ER-alpha         -0.7388394
## 6 TCGA-A2-A04P-01A      ESR1     ER-alpha         -1.3375791
```

### Perform the correlation

First we take the inner join of this data:

```r
library(dplyr)

data = inner_join(expressionData, proteinData)
```

```
## Joining by: "SampleBarcode"
```

```r
dim(data)
```

```
## [1] 403   6
```

```r
head(arrange(data, normalized_count))
```

```
##      SampleBarcode HGNC_gene_symbol normalized_count Gene_Name
## 1 TCGA-E2-A158-01A             ESR1           4.4555      ESR1
## 2 TCGA-AN-A0G0-01A             ESR1           5.8621      ESR1
## 3 TCGA-AR-A1AH-01A             ESR1           6.0368      ESR1
## 4 TCGA-AR-A0TP-01A             ESR1           6.1550      ESR1
## 5 TCGA-A2-A0D0-01A             ESR1           6.1779      ESR1
## 6 TCGA-A2-A04U-01A             ESR1           6.7541      ESR1
##   Protein_Name protein_expression
## 1     ER-alpha          -3.570038
## 2     ER-alpha          -3.873648
## 3     ER-alpha          -3.604598
## 4     ER-alpha          -3.338214
## 5     ER-alpha          -3.886696
## 6     ER-alpha          -2.897148
```

And run a spearman correlation on it:

```r
cor(x=data$normalized_count, y=data$protein_expression, method="spearman")
```

```
## [1] 0.9075414
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
 [5] MASS_7.3-40      munsell_0.4.2    colorspace_1.2-6 lattice_0.20-31 
 [9] R6_2.1.1         stringr_1.0.0    httr_1.0.0       plyr_1.8.3      
[13] tools_3.2.0      parallel_3.2.0   grid_3.2.0       gtable_0.1.2    
[17] DBI_0.3.1        lazyeval_0.1.10  assertthat_0.1   digest_0.6.8    
[21] Matrix_1.2-0     formatR_1.2      reshape2_1.4.1   curl_0.9.3      
[25] mime_0.4         evaluate_0.7.2   labeling_0.3     stringi_0.5-5   
[29] jsonlite_0.9.17  markdown_0.7.7   proto_0.3-10    
```
