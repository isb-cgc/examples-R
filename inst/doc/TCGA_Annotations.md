# TCGA Annotations
The goal of this notebook is to introduce you to the TCGA Annotations BigQuery table. You can find more detail about Annotations on the TCGA Wiki, but the key things to know are:
an annotation can refer to any "type" of TCGA "item" (eg patient, sample, portion, slide, analyte or aliquot), and
each annotation has a "classification" and a "category", both of which are drawn from controlled vocabularies.
The current set of annotation classifications includes: Redaction, Notification, CenterNotification, and Observation. The authority for Redactions and Notifications is the BCR (Biospecimen Core Resource), while CenterNotifications can come from any of the data-generating centers (GSC or GCC), and Observations from any authorized TCGA personnel. Within each classification type, there are several categories.
We will look at these further by querying directly on the Annotations table.
Note that annotations about patients, samples, and aliquots are separate from the clinical, biospecimen, and molecular data, and most patients, samples, and aliquots do not in fact have any annotations associated with them. It can be important, however, when creating a cohort or analyzing the molecular data associated with a cohort, to check for the existence of annotations.

In order to work with BigQuery, you need to import the bigrquery library and you need to know the name(s) of the table(s) you are going to be working with:


```r
require(bigrquery) || install.packages("bigrquery")
```

```
## [1] TRUE
```

```r
annotationsTable <- "[isb-cgc:tcga_201510_alpha.Annotations]"
```

From now on, we will occationally refer to this table using this variable annotationsTable, but we could just as well explicitly give the table name each time.
Let's start by taking a look at the table schema:


```r
querySql <- paste("SELECT * FROM ",annotationsTable," limit 1", sep="")
result <- query_exec(querySql, project=project)
data.frame(Columns=colnames(result))
```

```
##                     Columns
## 1              annotationId
## 2      annotationCategoryId
## 3    annotationCategoryName
## 4  annotationClassification
## 5        annotationNoteText
## 6                     Study
## 7              itemTypeName
## 8               itemBarcode
## 9            AliquotBarcode
## 10       ParticipantBarcode
## 11            SampleBarcode
## 12                dateAdded
## 13              dateCreated
## 14               dateEdited
```

## Item Types
Most of the schema fields come directly from the TCGA Annotations. First and foremost, an annotation is associated with an itemType, as described above. This can be a patient, an aliquot, etc. Let's see what the breakdown is of annotations according to item-type:


```r
querySql <- "
SELECT itemTypeName, COUNT(*) AS n
FROM [isb-cgc:tcga_201510_alpha.Annotations]
GROUP BY itemTypeName
ORDER BY n DESC"

result <- query_exec(querySql, project=project)
head(result)
```

```
##      itemTypeName     n
## 1         Aliquot 12928
## 2 Shipped Portion  1749
## 3         Patient  1378
## 4         Analyte   789
## 5           Slide   552
## 6          Sample   114
```

The length of the barcode in the itemBarcode field will depend on the value in the itemTypeName field: if the itemType is 'Patient', then the barcode will be something like TCGA-E2-A15J, whereas if the itemType is "Aliquot", the barcode will be a full-length barcode, eg TCGA-E2-A15J-10A-01D-a12N-01.

## Annotation Classifications and Categories
The next most important pieces of information about an annotation are the "classification" and "category". Each of these comes from a controlled vocabulary and each "classification" has a specific set of allowed "categories".
One important thing to understand is that if an aliquot carries some sort of disqualifying annotation, in general all other data from other samples or aliquots associated with that same patient should still be usable. On the other hand, if a patient carries some sort of disqualifying annotation, then that information should be considered prior to using any of the samples or aliquots derived from that patient.

To illustrate this, let's look at the most frequent annotation classifications and categories when the itemType is Patient:


```r
querySql <- "
SELECT
  annotationClassification,
  annotationCategoryName,
  COUNT(*) AS n
FROM
  [isb-cgc:tcga_201510_alpha.Annotations]
WHERE
  ( itemTypeName='Patient' )
GROUP BY
  annotationClassification,
  annotationCategoryName
HAVING ( n >= 50 )
ORDER BY
  n DESC"

result <- query_exec(querySql, project=project)
head(result)
```

```
##   annotationClassification
## 1             Notification
## 2             Notification
## 3             Notification
## 4             Notification
## 5             Notification
## 6             Notification
##                                                        annotationCategoryName
## 1                                                            Prior malignancy
## 2                                                   Alternate sample pipeline
## 3 History of unacceptable prior treatment related to a prior/other malignancy
## 4                                                      Synchronous malignancy
## 5                                                         Neoadjuvant therapy
## 6                                                        Item is noncanonical
##     n
## 1 407
## 2 200
## 3 139
## 4 110
## 5 102
## 6  81
```

The results of the previous query indicate that the majority of patient-level annotations are "Notifications", most frequently regarding prior malignancies. In most TCGA publications, "history of unacceptable prior treatment" and "item is noncanonical" notifications are treated as disqualifying annotations, and all data associated with those patients is not used in any analysis.
Let's make a slight modification to the last query to see what types of annotation categories and classifications we see when the item type is not patient:


```r
querySql <- "
SELECT
  annotationClassification,
  annotationCategoryName,
  itemTypeName,
  COUNT(*) AS n
FROM
  [isb-cgc:tcga_201510_alpha.Annotations]
WHERE
  ( itemTypeName!='Patient' )
GROUP BY
  annotationClassification,
  annotationCategoryName,
  itemTypeName
HAVING ( n >= 50 )
ORDER BY
  n DESC"

result <- query_exec(querySql, project=project)
head(result)
```

```
##   annotationClassification annotationCategoryName    itemTypeName     n
## 1       CenterNotification       Item flagged DNU         Aliquot 12303
## 2             Notification   Item is noncanonical Shipped Portion  1741
## 3             Notification   Item is noncanonical           Slide   541
## 4             Notification   Item is noncanonical         Analyte   464
## 5              Observation                General         Analyte   179
## 6       CenterNotification       Center QC failed         Aliquot   149
```

The results of the previous query indicate that the vast majority of annotations are at the aliquot level, and more specifically were submitted by one of the data-generating centers, indicating that the data derived from that aliquot is "DNU" (Do Not Use). In general, this should not affect any other aliquots derived from the same sample or any other samples derived from the same patient.
We see in the output of the previous query that a Notification that an "Item is noncanonical" can be applied to different types of items (eg slides and analytes). Let's investigate this a little bit further, for example let's count up these types of annotations by study (ie tumor-type):


```r
querySql <- "
SELECT
  Study,
  COUNT(*) AS n
FROM
  [isb-cgc:tcga_201510_alpha.Annotations]
WHERE
  ( annotationCategoryName='Item is noncanonical' )
GROUP BY
  Study
ORDER BY
  n DESC"

result <- query_exec(querySql, project=project)
head(result)
```

```
##   Study   n
## 1    OV 743
## 2   GBM 519
## 3  KIRC 455
## 4  COAD 314
## 5  LUAD 238
## 6  LUSC 231
```

and now let's pick one of these tumor types, and delve a little bit further:


```r
querySql <- "
SELECT
  itemTypeName,
  COUNT(*) AS n
FROM
  [isb-cgc:tcga_201510_alpha.Annotations]
WHERE
  ( annotationCategoryName='Item is noncanonical'
    AND Study='OV' )
GROUP BY
  itemTypeName
ORDER BY
  n DESC"

result <- query_exec(querySql, project=project)
head(result)
```

```
##      itemTypeName   n
## 1           Slide 409
## 2 Shipped Portion 220
## 3         Analyte 110
## 4         Patient   3
## 5          Sample   1
```

## Barcodes
As described above, an annotation is specific to a single TCGA "item" and the fields itemTypeName and itemBarcode are the most important keys to understanding which TCGA item carries the annotation. Because we use the fields ParticipantBarcode, SampleBarcode, and AliquotBarcode throughout our other TCGA BigQuery tables, we have added them to this table as well, but they should be interpreted with some care: when an annotation is specific to an aliquot (ie itemTypeName="Aliquot"), the ParticipantBarcode, SampleBarcode, and AliquotBarcode fields will all be set, but this should not be interpreted to mean that the annotation applies to all data derived from that patient.

This will be illustrated with the following two queries which extract information pertaining to a few specific patients:


```r
querySql <- "
SELECT
 Study,
 itemTypeName,
 itemBarcode,
 annotationCategoryName,
 annotationClassification,
 ParticipantBarcode,
 SampleBarcode,
 AliquotBarcode,
 LENGTH(itemBarcode) AS n
FROM
  [isb-cgc:tcga_201510_alpha.Annotations]
WHERE
  ( ParticipantBarcode='TCGA-61-1916' )
ORDER BY n ASC"

result <- query_exec(querySql, project=project)
head(result)
```

```
##   Study itemTypeName          itemBarcode annotationCategoryName
## 1    OV      Patient         TCGA-61-1916 Item in special subset
## 2    OV      Analyte TCGA-61-1916-01A-01T   Item is noncanonical
## 3    OV      Analyte TCGA-61-1916-01A-01W   Item is noncanonical
## 4    OV      Analyte TCGA-61-1916-02A-01W   Item is noncanonical
## 5    OV      Analyte TCGA-61-1916-11A-01W   Item is noncanonical
## 6    OV      Analyte TCGA-61-1916-01A-01G   Item is noncanonical
##   annotationClassification ParticipantBarcode    SampleBarcode
## 1             Notification       TCGA-61-1916             <NA>
## 2             Notification       TCGA-61-1916 TCGA-61-1916-01A
## 3             Notification       TCGA-61-1916 TCGA-61-1916-01A
## 4             Notification       TCGA-61-1916 TCGA-61-1916-02A
## 5             Notification       TCGA-61-1916 TCGA-61-1916-11A
## 6             Notification       TCGA-61-1916 TCGA-61-1916-01A
##   AliquotBarcode  n
## 1           <NA> 12
## 2           <NA> 20
## 3           <NA> 20
## 4           <NA> 20
## 5           <NA> 20
## 6           <NA> 20
```



```r
querySql <- "
SELECT
 Study,
 itemTypeName,
 itemBarcode,
 annotationCategoryName,
 annotationClassification,
 ParticipantBarcode,
 SampleBarcode,
 AliquotBarcode,
 LENGTH(itemBarcode) AS n
FROM
  [isb-cgc:tcga_201510_alpha.Annotations]
WHERE
  ( ParticipantBarcode='TCGA-GN-A261' )
ORDER BY n ASC"

result <- query_exec(querySql, project=project)
head(result)
```

```
##   Study itemTypeName  itemBarcode        annotationCategoryName
## 1  SKCM      Patient TCGA-GN-A261 Tumor tissue origin incorrect
## 2  SKCM      Patient TCGA-GN-A261           Neoadjuvant therapy
##   annotationClassification ParticipantBarcode SampleBarcode AliquotBarcode
## 1                Redaction       TCGA-GN-A261          <NA>           <NA>
## 2             Notification       TCGA-GN-A261          <NA>           <NA>
##    n
## 1 12
## 2 12
```

As you can see in the results returned from the previous two queries, the SampleBarcode and the AliquotBarcode fields may or may not be filled in, depending on the itemTypeName.


```r
querySql <- "
SELECT
 Study,
 itemTypeName,
 itemBarcode,
 annotationCategoryName,
 annotationClassification,
 annotationNoteText,
 ParticipantBarcode,
 SampleBarcode,
 AliquotBarcode,
 LENGTH(itemBarcode) AS n
FROM
  [isb-cgc:tcga_201510_alpha.Annotations]
WHERE
  ( ParticipantBarcode='TCGA-RS-A6TP' )
ORDER BY n ASC"

result <- query_exec(querySql, project=project)
head(result)
```

```
##   Study itemTypeName          itemBarcode annotationCategoryName
## 1  HNSC      Analyte TCGA-RS-A6TP-10A-01D                General
##   annotationClassification
## 1              Observation
##                                                                                                                                                                                                            annotationNoteText
## 1 DNA analyte UUID: 8304F61F-C217-4B9F-BA64-6486DA54E6C8 was involved in an extraction protocol deviation wherein an additional column purification step was used as a means of buffer exchange on the column-eluted analyte.
##   ParticipantBarcode    SampleBarcode AliquotBarcode  n
## 1       TCGA-RS-A6TP TCGA-RS-A6TP-10A           <NA> 20
```

In this example, there is just one annotation relevant to this particular patient, and one has to look at the annotationNoteText to find out what the potential issue may be with this particular analyte. Any aliquots derived from this blood-normal analyte might need to be used with care.
