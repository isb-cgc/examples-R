# Expression and Methylation Correlation

In this example, we will look at the correlation between mRNAseq-based gene expression and DNA methylation data.  We will do this using two molecular data tables from the isb-cgc:tcga_201507_alpha dataset and a cohort table from the isb-cgc:tcga_cohorts dataset.

NOTE: I think I will rework and/or eliminate this particular example, but am just going through it now to make sure I understand it and it works as expected.


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
expressionTable = "isb-cgc:tcga_201507_alpha.mRNA_UNC_HiSeq_RSEM"
methylationTable = "isb-cgc:tcga_201507_alpha.DNA_Methylation_betas"
# Add any additional clauses to be applied in WHERE to limit the methylation data further.
# (These specific filters are used here just to make the query run faster.  If a query returns
#  very large results, they may need to be handled differently.  This query should take < 20s)
andWhere = "AND SampleTypeLetterCode = 'TP' AND Study = 'CESC' AND CHR = '9'"
# Do not correlate unless there are at least this many observations available:
minNumObs = 30

# Now we are ready to run the query.  (Should return 6110 rows.)
result = DisplayAndDispatchQuery(
     file.path(sqlDir, "expression-methylation-correlation.sql"),
               project=project,
               replacements=list("_EXPRESSION_TABLE_"=expressionTable,
                                 "_METHYLATION_TABLE_"=methylationTable,
                                 "_AND_WHERE_"=andWhere,
                                 "_MINIMUM_NUMBER_OF_OBSERVATIONS_"=minNumObs))
```

```
# Compute the correlation between expression and methylation data.

SELECT
  HGNC_gene_symbol,
  Probe_ID,
  COUNT(1) AS num_observations,
  CORR(log2_count, Beta_Value) AS correlation,
FROM (
  # We select the sample-barcode, gene-symbol, gene-expression, probe-id, and beta-value
  # from a "JOIN" of the gene expression data and the methylation data.  Note that we log-
  # transform the expression since the value in the table is a normalized_count value. 
  SELECT
    expr.SampleBarcode,
    HGNC_gene_symbol,
    LOG2(normalized_count+1) AS log2_count,
    Probe_ID,
    Beta_value
  FROM
    [isb-cgc:tcga_201507_alpha.mRNA_UNC_HiSeq_RSEM] AS expr
  JOIN EACH ( FLATTEN ( (
        # We select the sample-barcode, sample-type, study-name, probe-id, beta-value, and gene-symbol
        # from the results of a "JOIN" of the methylation data and the methylation annotation tables
        # which are joined on the CpG probe id that exists in both tables.  Note that we need to 
        # FLATTEN this because the UCSC.RefGene information is a (potentially) repeated field.
        SELECT
          SampleBarcode,
          SampleTypeLetterCode,
          Study,
          Probe_ID,
          Beta_Value,
	  CHR,
          UCSC.RefGene_Name
        FROM
          [isb-cgc:tcga_201507_alpha.DNA_Methylation_betas] AS methData
        JOIN EACH [isb-cgc:platform_reference.methylation_annotation] AS methAnnot
        ON
          methData.Probe_ID = methAnnot.Name
          # We require that the gene-symbol not be null.
        WHERE
          UCSC.RefGene_Name IS NOT NULL
          # Optionally add clause here to limit the query to a particular
          # sample types and/or study and/or chromosome.
          AND SampleTypeLetterCode = 'TP' AND Study = 'CESC' AND CHR = '9'
        ), UCSC.RefGene_Name ) ) AS methyl
  ON
    methyl.UCSC.RefGene_Name = expr.HGNC_gene_symbol
    AND methyl.SampleBarcode = expr.SampleBarcode )
GROUP BY
  HGNC_gene_symbol,
  Probe_ID
HAVING
  num_observations >= 30 AND correlation > -2.
ORDER BY
  correlation ASC
```
Number of rows returned by this query: 6110.

The result is a table with one row for each (gene,CpG-probe) pair for which at least 30 data values exist that meet the requirements in the "andWhere" clause.  The (gene,CpG-probe) pair is defined by a gene symbol and a CpG-probe ID.  In many cases, there may be multiple CpG probes associated with a single gene.


```r
# Most negative correlation should be PHYHD1 cg14299940 n=903, correlation = -0.8018487
head(result)
```

```
##   HGNC_gene_symbol   Probe_ID num_observations correlation
## 1           PHYHD1 cg14299940              903  -0.8018487
## 2            INSL6 cg13504907              301  -0.7800666
## 3            BICD2 cg02929681              602  -0.7252331
## 4             CRAT cg22192879              903  -0.6979521
## 5            BICD2 cg13683626              602  -0.6953693
## 6            BICD2 cg14181777              602  -0.6806368
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
gene = "PHYHD1"

#FIXME: there probably needs to be some additional "AND WHERE" filtering here... 
#because right now this is returning 10252 rows...
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
FROM [isb-cgc:tcga_201507_alpha.mRNA_UNC_HiSeq_RSEM]
WHERE
  HGNC_gene_symbol = 'PHYHD1'
ORDER BY
  SampleBarcode
```
Number of rows returned by this query: 10252.


```r
head(expressionData)
```

```
##      SampleBarcode HGNC_gene_symbol normalized_count
## 1 TCGA-02-0047-01A           PHYHD1          84.8213
## 2 TCGA-02-0055-01A           PHYHD1         220.8830
## 3 TCGA-02-2483-01A           PHYHD1          67.9683
## 4 TCGA-02-2485-01A           PHYHD1          39.0476
## 5 TCGA-02-2486-01A           PHYHD1         252.4390
## 6 TCGA-04-1348-01A           PHYHD1        1181.5004
```

### Retrieve Methylation Data

Then we retrieve the methylation data for a particular probe.


```r
# Set the desired probe to query.
probe = "cg14299940"

# Be sure to apply the same additional clauses to the WHERE to limit the methylation data further.

#FIXME the andWhere clause breaks here because the previously it was being applied to the result of a JOIN
methylationData = DisplayAndDispatchQuery(file.path(sqlDir, "methylation-data.sql"),
                                          project=project,
                                          replacements=list("_METHYLATION_TABLE_"=methylationTable,
                                                            "_AND_WHERE_"="",
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
FROM [isb-cgc:tcga_201507_alpha.DNA_Methylation_betas]
WHERE
  Probe_ID = 'cg14299940'
  # Optionally add clause here to limit the query to a particular
  # sample types and/or studies.
  
ORDER BY
  SampleBarcode,
  SampleTypeLetterCode,
  Study
```

Number of rows returned by this query: 9626.


```r
head(methylationData)
```

```
##      SampleBarcode SampleTypeLetterCode Study   Probe_ID Beta_Value
## 1 TCGA-05-4384-01A                   TP  LUAD cg14299940       0.17
## 2 TCGA-05-4390-01A                   TP  LUAD cg14299940       0.09
## 3 TCGA-05-4396-01A                   TP  LUAD cg14299940       0.49
## 4 TCGA-05-4405-01A                   TP  LUAD cg14299940       0.13
## 5 TCGA-05-4410-01A                   TP  LUAD cg14299940       0.38
## 6 TCGA-05-4415-01A                   TP  LUAD cg14299940       0.09
```

### Perform the correlation

First we take the inner join of this data:

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
## 1 TCGA-05-4384-01A           PHYHD1         263.2865                   TP
## 2 TCGA-05-4390-01A           PHYHD1         515.4080                   TP
## 3 TCGA-05-4396-01A           PHYHD1          15.9212                   TP
## 4 TCGA-05-4405-01A           PHYHD1         412.1566                   TP
## 5 TCGA-05-4410-01A           PHYHD1         365.4434                   TP
## 6 TCGA-05-4415-01A           PHYHD1         144.0000                   TP
##   Study   Probe_ID Beta_Value
## 1  LUAD cg14299940       0.17
## 2  LUAD cg14299940       0.09
## 3  LUAD cg14299940       0.49
## 4  LUAD cg14299940       0.13
## 5  LUAD cg14299940       0.38
## 6  LUAD cg14299940       0.09
```

And run a pearson correlation on it:

```r
cor(x=data$normalized_count, y=data$Beta_Value, method="pearson")
```

```
## [1] -0.4054692
```

And we can see that we have reproduced one of our results from BigQuery.
#FIXME but we have not ;)  the correlation now is -0.405 on 8756 samples rather than just 603...

## Provenance

```r
sessionInfo()
```

```
R version 3.2.2 (2015-08-14)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.10.5 (Yosemite)

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] knitr_1.11         scales_0.3.0       ggplot2_1.0.1     
[4] bigrquery_0.1.0    dplyr_0.4.3        ISBCGCExamples_0.1

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.1      magrittr_1.5     MASS_7.3-43      munsell_0.4.2   
 [5] colorspace_1.2-6 R6_2.1.1         stringr_1.0.0    httr_1.0.0      
 [9] plyr_1.8.3       tools_3.2.2      parallel_3.2.2   grid_3.2.2      
[13] gtable_0.1.2     DBI_0.3.1        lazyeval_0.1.10  assertthat_0.1  
[17] digest_0.6.8     reshape2_1.4.1   formatR_1.2.1    curl_0.9.3      
[21] mime_0.4         evaluate_0.8     labeling_0.3     stringi_0.5-5   
[25] jsonlite_0.9.17  markdown_0.7.7   proto_0.3-10    
```
