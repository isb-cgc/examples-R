# microRNA expression (BCGSC RPKM)

The goal of this notebook is to introduce you to the microRNA expression BigQuery table.
This table contains all available TCGA Level-3 microRNA expression data produced by BCGSC's microRNA pipeline using the Illumina HiSeq platform, as of October 2015. (Actual archive dates range from April 2013 to June 2015.) The most recent archive (eg bcgsc.ca_THCA.IlluminaHiSeq_miRNASeq.Level_3.1.9.0) for each of the 32 tumor types was downloaded from the DCC, and data extracted from all files matching the pattern %.isoform.quantification.txt. The isoform-quantification values were then processed through a Perl script provided by BCGSC which produces normalized expression levels for mature microRNAs. Each of these mature microRNAs is identified by name (eg hsa-mir-21) and by MIMAT accession number (eg MIMAT0000076).

In order to work with BigQuery, you need to import the bigrquery library and you need to know the name(s) of the table(s) you are going to be working with:


```r
require(bigrquery) || install.packages("bigrquery")
```

```
## [1] TRUE
```

```r
require(ggplot2) || install.packages("ggplot2")
```

```
## [1] TRUE
```

```r
library(ISBCGCExamples)

miRNATable <- "[isb-cgc:tcga_201510_alpha.miRNA_expression]"
```

From now on, we will refer to this table using this variable ($miRNA_BQtable), but we could just as well explicitly give the table name each time.
Let's start by taking a look at the table schema:


```r
querySql <- paste("SELECT * FROM ",miRNATable," limit 1", sep="")
result <- query_exec(querySql, project=project)
data.frame(Columns=colnames(result))
```

```
##                Columns
## 1   ParticipantBarcode
## 2        SampleBarcode
## 3       AliquotBarcode
## 4 SampleTypeLetterCode
## 5                Study
## 6             Platform
## 7             mirna_id
## 8      mirna_accession
## 9     normalized_count
```

Now let's count up the number of unique patients, samples and aliquots mentioned in this table. We will do this by defining a very simple parameterized query. (Note that when using a variable for the table name in the FROM clause, you should not also use the square brackets that you usually would if you were specifying the table name as a string.)


```r
for (x in c("ParticipantBarcode", "SampleBarcode", "AliquotBarcode")) {
  querySql <- paste("SELECT COUNT(DISTINCT(",x,"), 25000) AS n ",
                    "FROM ",miRNATable)
  result <- query_exec(querySql, project=project)
  cat(x, ": ", result[[1]], "\n")
}
```

```
## ParticipantBarcode :  9992 
## SampleBarcode :  10712 
## AliquotBarcode :  10773
```

We can do the same thing to look at how many unique microRNAs exist in the table:


```r
for (x in c("mirna_id", "mirna_accession")) {
  querySql <- paste("SELECT COUNT(DISTINCT(",x,"), 25000) AS n ",
                    "FROM ",miRNATable)
  result <- query_exec(querySql, project=project)
  cat(x, ": ", result[[1]], "\n")
}
```

```
## mirna_id :  965 
## mirna_accession :  1222
```


These counts show that the mirna_id field is not a unique identifier and should be used in combination with the MIMAT accession number.
Another thing to note about this table is that these expression values are obtained from two different platforms -- approximately 15% of the data is from the Illumina GA platform, and 85% from the Illumina HiSeq:


```r
querySql <- "SELECT
              Platform,
              COUNT(*) AS n
            FROM
              [isb-cgc:tcga_201510_alpha.miRNA_expression]
            GROUP BY
              Platform
            ORDER BY
              n DESC"
query_exec(querySql, project=project)
```

```
##        Platform        n
## 1 IlluminaHiSeq 11489244
## 2    IlluminaGA  1994304
```

```r
sessionInfo()
```
