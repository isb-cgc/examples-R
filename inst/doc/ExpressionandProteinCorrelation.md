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
  feat1.gene AS gene,
  feat1.protein AS protein,
  COUNT(1) AS num_observations,
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
HAVING
  num_observations >= 30
ORDER BY
  spearman_corr DESC
```
Number of rows returned by this query: 213.

The result is one correlation value per row of data, each of which corresponds to  . . . MORE HERE

```r
head(result)
```

```
##     gene         protein num_observations spearman_corr
## 1   GAB2            GAB2               60     0.7476542
## 2    SYK             Syk              180     0.7455379
## 3 NOTCH3          Notch3               49     0.7135100
## 4  ANXA1       Annexin-1               94     0.6839287
## 5  PEA15          PEA-15               49     0.6837365
## 6  PRKCA PKC-alpha_pS657              180     0.6740750
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
gene = "ANXA1"

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
  HGNC_gene_symbol = 'ANXA1'
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
## 1 TCGA-02-2485-01A            ANXA1         5974.762
## 2 TCGA-05-4244-01A            ANXA1        17585.828
## 3 TCGA-05-4250-01A            ANXA1        10903.424
## 4 TCGA-06-0138-01A            ANXA1         6590.878
## 5 TCGA-06-0157-01A            ANXA1         3268.640
## 6 TCGA-06-0158-01A            ANXA1         6452.934
```

### Retrieve Protein Data

Then we retrieve the protein data for a particular gene.


```r
protein = "Annexin-1"

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
  Gene_Name = 'ANXA1'
  AND Protein_Name = 'Annexin-1'
  AND SampleBarcode IN (
  SELECT
    sample_barcode
  FROM
    [isb-cgc:test.cohort_14jun2015] )
ORDER BY
  SampleBarcode
```
Number of rows returned by this query: 96.


```r
head(proteinData)
```

```
##      SampleBarcode Gene_Name Protein_Name protein_expression
## 1 TCGA-2Z-A9J7-01A     ANXA1    Annexin-1         0.12343802
## 2 TCGA-4A-A93X-01A     ANXA1    Annexin-1        -0.57261007
## 3 TCGA-A4-8311-01A     ANXA1    Annexin-1         0.04999981
## 4 TCGA-B9-5155-01A     ANXA1    Annexin-1         0.92203965
## 5 TCGA-B9-7268-01A     ANXA1    Annexin-1         0.16892302
## 6 TCGA-B9-A69E-01A     ANXA1    Annexin-1        -0.42145909
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
dim(data)
```

```
## [1] 94  6
```

```r
head(data)
```

```
##      SampleBarcode HGNC_gene_symbol normalized_count Gene_Name
## 1 TCGA-2Z-A9J7-01A            ANXA1        7509.1798     ANXA1
## 2 TCGA-4A-A93X-01A            ANXA1         209.4447     ANXA1
## 3 TCGA-A4-8311-01A            ANXA1        5779.9858     ANXA1
## 4 TCGA-B9-5155-01A            ANXA1       12123.4097     ANXA1
## 5 TCGA-B9-7268-01A            ANXA1        6643.2432     ANXA1
## 6 TCGA-B9-A69E-01A            ANXA1        1137.3156     ANXA1
##   Protein_Name protein_expression
## 1    Annexin-1         0.12343802
## 2    Annexin-1        -0.57261007
## 3    Annexin-1         0.04999981
## 4    Annexin-1         0.92203965
## 5    Annexin-1         0.16892302
## 6    Annexin-1        -0.42145909
```

First we take the inner join of this data and run a spearman correlation on it.

```r
cor(x=data$normalized_count, y=data$protein_expression, method="spearman")
```

```
## [1] 0.6976339
```
Notice that the value does not match the result from BigQuery.

The reason for this is that the RANK in the BigQuery method is run over all the observations, not just those with a matching observation in both tables.  So let's redo this in R.

First we perform the RANK operation before doing the inner join.

```r
expressionData = mutate(expressionData, expr_rank=rank(normalized_count))
proteinData = mutate(proteinData, prot_rank=rank(protein_expression))
```

Then we do the join and run a spearman correlation on the ranked values.

```r
data = inner_join(expressionData, proteinData)
```

```
## Joining by: "SampleBarcode"
```

```r
dim(data)
```

```
## [1] 94  8
```

```r
head(data)
```

```
##      SampleBarcode HGNC_gene_symbol normalized_count expr_rank Gene_Name
## 1 TCGA-2Z-A9J7-01A            ANXA1        7509.1798       360     ANXA1
## 2 TCGA-4A-A93X-01A            ANXA1         209.4447        13     ANXA1
## 3 TCGA-A4-8311-01A            ANXA1        5779.9858       313     ANXA1
## 4 TCGA-B9-5155-01A            ANXA1       12123.4097       422     ANXA1
## 5 TCGA-B9-7268-01A            ANXA1        6643.2432       343     ANXA1
## 6 TCGA-B9-A69E-01A            ANXA1        1137.3156       138     ANXA1
##   Protein_Name protein_expression prot_rank
## 1    Annexin-1         0.12343802        53
## 2    Annexin-1        -0.57261007        16
## 3    Annexin-1         0.04999981        47
## 4    Annexin-1         0.92203965        90
## 5    Annexin-1         0.16892302        55
## 6    Annexin-1        -0.42145909        26
```

```r
cor(x=data$expr_rank, y=data$prot_rank, method="pearson")
```

```
## [1] 0.6839287
```
Now the results match those for BigQuery.

But perhaps the BigQuery query should be running the RANK operation only after data has been removed per the inner join.

```r
resultV2 = DisplayAndDispatchQuery(file.path(sqlDir, "protein-mrna-spearman-correlation-V2.sql"),
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
    RANK() OVER (PARTITION BY gene ORDER BY log2_count ASC) AS expr_rank,
    RANK() OVER (PARTITION BY protein ORDER BY protein_expression ASC) AS prot_rank,
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
        [isb-cgc:tcga_data_open.Protein]
      WHERE
        SampleBarcode IN (
        SELECT
          sample_barcode
        FROM
          [isb-cgc:test.cohort_14jun2015] ) ) feat1
    JOIN EACH (
      SELECT
        SampleBarcode,
        HGNC_gene_symbol,
        LOG2(normalized_count+1) AS log2_count
      FROM
        [isb-cgc:tcga_data_open.mRNA_UNC_HiSeq_RSEM]
      WHERE
        SampleBarcode IN (
        SELECT
          sample_barcode
        FROM
          [isb-cgc:test.cohort_14jun2015] ) ) feat2
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
Number of rows returned by this query: 213.

TODO This result does not yet match our earlier spearman correlation via R.

```r
head(resultV2, n=12)
```

```
##      gene         protein num_observations spearman_corr
## 1     SYK             Syk              180     0.7438830
## 2    GAB2            GAB2               60     0.7407510
## 3  NOTCH3          Notch3               49     0.6939796
## 4   ANXA1       Annexin-1               94     0.6886659
## 5   PEA15          PEA-15               49     0.6813926
## 6   PRKCA PKC-alpha_pS657              180     0.6784358
## 7     SRC             Src              180     0.6746813
## 8    BCL2           Bcl-2              180     0.6397852
## 9   PREX1           PREX1              117     0.6279054
## 10     AR              AR              180     0.6277496
## 11   ESR1        ER-alpha              180     0.6253135
## 12  MYH11           MYH11              125     0.5679508
```


```r
library(ggplot2)

# Histogram overlaid with kernel density curve
ggplot(resultV2, aes(x=spearman_corr)) + 
    geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   binwidth=.05,
                   colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
```

<img src="figure/spearman_density2-1.png" title="plot of chunk spearman_density2" alt="plot of chunk spearman_density2" style="display: block; margin: auto;" />


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
[1] knitr_1.10.5       ggplot2_1.0.1      bigrquery_0.1.0   
[4] dplyr_0.4.2        ISBCGCExamples_0.1

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.0      rstudioapi_0.3.1 magrittr_1.5     MASS_7.3-40     
 [5] munsell_0.4.2    colorspace_1.2-6 R6_2.1.1         stringr_1.0.0   
 [9] httr_1.0.0       plyr_1.8.2       tools_3.2.0      parallel_3.2.0  
[13] grid_3.2.0       gtable_0.1.2     DBI_0.3.1        lazyeval_0.1.10 
[17] assertthat_0.1   digest_0.6.8     reshape2_1.4.1   formatR_1.2     
[21] curl_0.9.3       mime_0.4         evaluate_0.7.2   labeling_0.3    
[25] stringi_0.5-5    scales_0.2.4     jsonlite_0.9.17  markdown_0.7.7  
[29] proto_0.3-10    
```
