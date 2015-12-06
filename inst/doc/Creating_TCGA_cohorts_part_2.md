# Creating TCGA cohorts (part 2)

This notebook will provide another example of building a cohort, this time based on the types of samples available.


```r
require(bigrquery) || install.packages("bigrquery")
```

```
## [1] TRUE
```

```r
bigrquery::list_tables("isb-cgc", "tcga_201510_alpha")
```

```
##  [1] "Annotations"            "Biospecimen_data"      
##  [3] "Clinical_data"          "Copy_Number_segments"  
##  [5] "DNA_Methylation_betas"  "Protein_RPPA_data"     
##  [7] "Somatic_Mutation_calls" "mRNA_BCGSC_HiSeq_RPKM" 
##  [9] "mRNA_UNC_HiSeq_RSEM"    "miRNA_expression"
```

Many different types of samples were obtained from the TCGA participants, and details about these samples are available in the Biospecimen data table. This next query shows how many samples exist of each type, as well as the full names and abbreviations of each type:


```r
querySql <- "
SELECT
  SampleType,
  SampleTypeLetterCode,
  COUNT(*) AS n
FROM
  [isb-cgc:tcga_201510_alpha.Biospecimen_data]
GROUP BY
  SampleType,
  SampleTypeLetterCode,
ORDER BY
  n DESC"

result <- query_exec(querySql, project=project)
head(result)
```

```
##                                        SampleType SampleTypeLetterCode
## 1                             Primary solid Tumor                   TP
## 2                            Blood Derived Normal                   NB
## 3                             Solid Tissue Normal                   NT
## 4                                      Metastatic                   TM
## 5 Primary Blood Derived Cancer - Peripheral Blood                   TB
## 6                           Recurrent Solid Tumor                   TR
##       n
## 1 10787
## 2  9382
## 3  2685
## 4   396
## 5   356
## 6    60
```

Note that there are many types of tumor samples: primary, metastatic, recurrent, etc, although the vast majority are samples from primary tumors. In the TCGA project, almost all tumor samples were assayed on multiple platforms for mRNA and miRNA expression, DNA methylation, DNA copy-number, and either exome- or whole-genome DNA sequence. For some tumor samples, protein activity was also measured using RPPA arrays. When available, adjacent "normal" tissue samples were also assayed on a subset of these platforms. The "blood normal" samples were primarily used only as a reference source of germline DNA in order to call somatic mutations.

We can do a similar counting exercise of the sample types represented in one of the molecular data tables, using one of the mRNA expression data tables:


```r
querySql <- "
SELECT
  SampleTypeLetterCode,
  COUNT(*) AS n
FROM (
  SELECT
    SampleBarcode,
    SampleTypeLetterCode
  FROM
    [isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM]
  GROUP BY
    SampleBarcode,
    SampleTypeLetterCode )
GROUP BY
  SampleTypeLetterCode
ORDER BY
  n DESC"

result <- query_exec(querySql, project=project)
head(result)
```

```
##   SampleTypeLetterCode    n
## 1                   TP 8976
## 2                   NT  715
## 3                   TM  393
## 4                   TB  173
## 5                   TR   45
## 6                  TAP   11
```

In this example, let's assume that we would like to do a study that requires a primary tumor sample and a matched-normal (adjacent) tissue sample. In order to find out which patients provided which types of samples, we need to query the Biospecimen data table. This next query module uses two sub-queries, one to get all patients with TP samples and another to get all patients with NT samples. The final query joins these two and returns a single list of patients.


```r
# get the list of patients who have provided TP samples
querySql <- "
SELECT
  ParticipantBarcode
FROM
  [isb-cgc:tcga_201510_alpha.Biospecimen_data]
WHERE
  ( SampleTypeLetterCode='TP' )
GROUP BY
  ParticipantBarcode
ORDER BY
  ParticipantBarcode"

tp_result <- query_exec(querySql, project=project)
head(tp_result)
```

```
##   ParticipantBarcode
## 1       TCGA-02-0001
## 2       TCGA-02-0003
## 3       TCGA-02-0004
## 4       TCGA-02-0006
## 5       TCGA-02-0007
## 6       TCGA-02-0009
```


```r
# now get a list of patients who have provided NT samples
querySql <- "
SELECT
  ParticipantBarcode
FROM
  [isb-cgc:tcga_201510_alpha.Biospecimen_data]
WHERE
  ( SampleTypeLetterCode='NT' )
GROUP BY
  ParticipantBarcode
ORDER BY
  ParticipantBarcode"

nt_result <- query_exec(querySql, project=project)
head(nt_result)
```

```
##   ParticipantBarcode
## 1       TCGA-01-0628
## 2       TCGA-01-0630
## 3       TCGA-01-0631
## 4       TCGA-01-0633
## 5       TCGA-01-0636
## 6       TCGA-01-0637
```


```r
# and finally join these two lists to get the intersection of these two lists
patients_both <- (intersect(tp_result$ParticipantBarcode, nt_result$ParticipantBarcode))
head(patients_both)
```

```
## [1] "TCGA-04-1335" "TCGA-04-1336" "TCGA-04-1337" "TCGA-04-1338"
## [5] "TCGA-04-1342" "TCGA-04-1346"
```

It might be interesting to find out what the distribution of tumor types is for this list of patients with matched tumor-normal sample pairs. We can define a new SQL module that refers to the results of a previously defined query as long as we pass that reference in when we call bq.Query():


```r
# now we'll use this list to find what types of tumors these patients
# belong to:

wrapSingleQuotes <- function(x) {
  paste("\'", paste(x, collapse="\',\'"), "\'", sep="")
}

querySql <- paste(
"SELECT
  Study,
  COUNT(*) AS n
FROM
  [isb-cgc:tcga_201510_alpha.Clinical_data]
WHERE
  ParticipantBarcode IN (", wrapSingleQuotes(patients_both) ,
  ")
GROUP BY
  Study
ORDER BY
  n DESC", sep="")

result <- query_exec(querySql, project=project)
head(result)
```

```
##   Study   n
## 1  KIRC 442
## 2  LUSC 253
## 3  LUAD 211
## 4  BRCA 162
## 5    OV 126
## 6  PRAD 118
```
