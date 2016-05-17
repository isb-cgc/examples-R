# UNC HiSeq mRNAseq gene expression (RSEM)
The goal of this notebook is to introduce you to the mRNAseq gene expression BigQuery table.
This table contains all available TCGA Level-3 gene expression data produced by UNC's RNAseqV2 pipeline using the Illumina HiSeq platform, as of October 2015. (Actual archive dates range from January 2013 to June 2015.) The most recent archive (eg unc.edu_BRCA.IlluminaHiSeq_RNASeqV2.Level_3.1.11.0) for each of the 33 tumor types was downloaded from the DCC, and data extracted from all files matching the pattern %.rsem.genes.normalized_results. Each of these raw “RSEM genes normalized results” files has two columns: gene_id and normalized_count. The gene_id string contains two parts: the gene symbol, and the Entrez gene ID, separated by | eg: TP53|7157. During ETL, the gene_id string is split and the gene symbol is stored in the original_gene_symbol field, and the Entrez gene ID is stored in the gene_id field. In addition, the Entrez ID is used to look up the current HGNC approved gene symbol, which is stored in the HGNC_gene_sybmol field.

In order to work with BigQuery, you need to import the bigrquery library and you need to know the name(s) of the table(s) you are going to be working with:


```r
require(bigrquery) || install.packages("bigrquery")
```

```
## [1] TRUE
```

```r
rnaTable <- "[isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM]"
```

From now on, occasionally, we will refer to this table using this variable rnaTable, but we could just as well explicitly give the table name each time.

Let's start by taking a look at the table schema:


```r
querySql <- paste("SELECT * FROM ",rnaTable," limit 1", sep="")
result <- query_exec(querySql, project=project)
data.frame(Columns=colnames(result))
```

```
##                 Columns
## 1    ParticipantBarcode
## 2         SampleBarcode
## 3        AliquotBarcode
## 4                 Study
## 5  SampleTypeLetterCode
## 6              Platform
## 7  original_gene_symbol
## 8      HGNC_gene_symbol
## 9               gene_id
## 10     normalized_count
```

Now let's count up the number of unique patients, samples and aliquots mentioned in this table. We will do this by defining a very simple parameterized query. (Note that when using a variable for the table name in the FROM clause, you should not also use the square brackets that you usually would if you were specifying the table name as a string.)


```r
for (x in c("ParticipantBarcode", "SampleBarcode", "AliquotBarcode")) {
  querySql <- paste("SELECT COUNT(DISTINCT(",x,"), 25000) AS n ",
                    "FROM ",rnaTable)
  result <- query_exec(querySql, project=project)
  cat(x, ": ", result[[1]], "\n")
}
```

```
## ParticipantBarcode :  9566 
## SampleBarcode :  10314 
## AliquotBarcode :  10343
```

We can do the same thing to look at how many unique gene symbols and gene ids exist in the table:


```r
for (x in c("original_gene_symbol", "HGNC_gene_symbol", "gene_id")) {
  querySql <- paste("SELECT COUNT(DISTINCT(",x,"), 25000) AS n ",
                    "FROM ",rnaTable)
  result <- query_exec(querySql, project=project)
  cat(x, ": ", result[[1]], "\n")
}
```

```
## original_gene_symbol :  20501 
## HGNC_gene_symbol :  20182 
## gene_id :  20531
```

Based on the counts, we can see that there are a few instances where the original gene symbol (from the underlying TCGA data file), or the HGNC gene symbol or the gene id (also from the original TCGA data file) is missing, but for the majority of genes, all three values should be available and for the most part the original gene symbol and the HGNC gene symbol that was added during ETL should all match up. This next query will generate the complete list of genes for which none of the identifiers are null, and where the original gene symbol and the HGNC gene symbol match. This list has over 18000 genes in it.


```r
querySql <- "
SELECT
  HGNC_gene_symbol,
  original_gene_symbol,
  gene_id
FROM
  [isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM]
WHERE
  ( original_gene_symbol IS NOT NULL
    AND HGNC_gene_symbol IS NOT NULL
    AND original_gene_symbol=HGNC_gene_symbol
    AND gene_id IS NOT NULL )
GROUP BY
  original_gene_symbol,
  HGNC_gene_symbol,
  gene_id
ORDER BY
  HGNC_gene_symbol"

result <- query_exec(querySql, project=project)
```

```r
head(result)
```

```
##   HGNC_gene_symbol original_gene_symbol gene_id
## 1             A1BG                 A1BG       1
## 2             A1CF                 A1CF   29974
## 3              A2M                  A2M       2
## 4            A2ML1                A2ML1  144568
## 5           A4GALT               A4GALT   53947
## 6            A4GNT                A4GNT   51146
```

We might also want to know how often the gene symbols do not agree:


```r
querySql <- "
SELECT
  HGNC_gene_symbol,
  original_gene_symbol,
  gene_id
FROM
  [isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM]
WHERE
  ( original_gene_symbol IS NOT NULL
    AND HGNC_gene_symbol IS NOT NULL
    AND original_gene_symbol!=HGNC_gene_symbol
    AND gene_id IS NOT NULL )
GROUP BY
  original_gene_symbol,
  HGNC_gene_symbol,
  gene_id
ORDER BY
  HGNC_gene_symbol"

result <- query_exec(querySql, project=project)
head(result)
```

```
##   HGNC_gene_symbol original_gene_symbol gene_id
## 1         A1BG-AS1           NCRNA00181  503538
## 2          A2M-AS1            LOC144571  144571
## 3           AACSP1                AACSL  729522
## 4          AADACP1            LOC201651  201651
## 5            AAED1              C9orf21  195827
## 6            AAMDC             C11orf67   28971
```

BigQuery is not just a "look-up" service -- you can also use it to perform calculations. In this next query, we take a look at the mean, standard deviation, and coefficient of variation for the expression of EGFR, within each tumor-type, as well as the number of primary tumor samples that went into each summary statistic.


```r
querySql <- "
SELECT
  Study,
  n,
  exp_mean,
  exp_sigma,
  (exp_sigma/exp_mean) AS exp_cv
FROM (
  SELECT
    Study,
    AVG(LOG2(normalized_count+1)) AS exp_mean,
    STDDEV_POP(LOG2(normalized_count+1)) AS exp_sigma,
    COUNT(AliquotBarcode) AS n
  FROM
    [isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM]
  WHERE
    ( SampleTypeLetterCode='TP'
      AND HGNC_gene_symbol='EGFR' )
  GROUP BY
    Study )
ORDER BY
  exp_sigma DESC"

result <- query_exec(querySql, project=project)
```

We can also easily move the gene-symbol out of the WHERE clause and into the SELECT and GROUP BY clauses and have BigQuery do this same calculation over all genes and all tumor types.


```r
querySql <- "
SELECT
  Study,
  HGNC_gene_symbol,
  n,
  exp_mean,
  exp_sigma,
  (exp_sigma/exp_mean) AS exp_cv
FROM (
  SELECT
    Study,
    HGNC_gene_symbol,
    AVG(LOG2(normalized_count+1)) AS exp_mean,
    STDDEV_POP(LOG2(normalized_count+1)) AS exp_sigma,
    COUNT(AliquotBarcode) AS n
  FROM
    [isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM]
  WHERE
    ( SampleTypeLetterCode='TP' )
  GROUP BY
    Study,
    HGNC_gene_symbol )
ORDER BY
  exp_sigma DESC
LIMIT 2000"

result <- query_exec(querySql, project=project)
head(result)
```

```
##   Study HGNC_gene_symbol   n  exp_mean exp_sigma    exp_cv
## 1  SKCM             KRT1 104  7.284893  6.374956 0.8750927
## 2  SKCM            KRT16 104  9.211635  6.234682 0.6768269
## 3  SKCM            KRT6A 104  9.866586  6.225797 0.6309981
## 4  PAAD            PNLIP 156  8.743529  6.204216 0.7095781
## 5  ESCA             KRT5 182 12.252972  6.059924 0.4945677
## 6  SKCM            KRT6C 104  8.383869  6.050728 0.7217107
```

```r
subResult <- result[result$exp_mean > 6 & result$n >= 200 & result$exp_cv > 0.5,]

head( subResult[order(subResult$exp_cv, decreasing=T),] )
```

```
##     Study HGNC_gene_symbol    n exp_mean exp_sigma    exp_cv
## 28   SARC           RPS4Y1  254 6.123702  5.546050 0.9056694
## 38   COAD             XIST  288 6.202007  5.456103 0.8797318
## 79   COAD           RPS4Y1  288 6.152626  5.143234 0.8359413
## 110  BLCA            GSTM1  407 6.007833  5.009420 0.8338148
## 56   KIRC             XIST  519 6.440970  5.275562 0.8190632
## 77   BRCA             CPB1 1095 6.531856  5.153108 0.7889194
```

Since the result of the previous query was quite large (over 600,000 rows representing ~20,000 genes x ~30 tumor types), we could load that table into bigquery for subsequent work. Instead we can explore the results in R.
