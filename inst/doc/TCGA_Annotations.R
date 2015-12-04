#' # TCGA Annotations
#' The goal of this notebook is to introduce you to the TCGA Annotations BigQuery table. You can find more detail about Annotations on the TCGA Wiki, but the key things to know are:
#' an annotation can refer to any "type" of TCGA "item" (eg patient, sample, portion, slide, analyte or aliquot), and
#' each annotation has a "classification" and a "category", both of which are drawn from controlled vocabularies.
#' The current set of annotation classifications includes: Redaction, Notification, CenterNotification, and Observation. The authority for Redactions and Notifications is the BCR (Biospecimen Core Resource), while CenterNotifications can come from any of the data-generating centers (GSC or GCC), and Observations from any authorized TCGA personnel. Within each classification type, there are several categories.
#' We will look at these further by querying directly on the Annotations table.
#' Note that annotations about patients, samples, and aliquots are separate from the clinical, biospecimen, and molecular data, and most patients, samples, and aliquots do not in fact have any annotations associated with them. It can be important, however, when creating a cohort or analyzing the molecular data associated with a cohort, to check for the existence of annotations.
#' 
#' In order to work with BigQuery, you need to import the bigrquery library and you need to know the name(s) of the table(s) you are going to be working with:
#' 
## ------------------------------------------------------------------------
require(bigrquery) || install.packages("bigrquery")

annotationsTable <- "[isb-cgc:tcga_201510_alpha.Annotations]"

#' 
#' From now on, we will occationally refer to this table using this variable annotationsTable, but we could just as well explicitly give the table name each time.
#' Let's start by taking a look at the table schema:
#' 
## ------------------------------------------------------------------------
querySql <- paste("SELECT * FROM ",annotationsTable," limit 1", sep="")
result <- query_exec(querySql, project=project)
data.frame(Columns=colnames(result))

#' 
#' ## Item Types
#' Most of the schema fields come directly from the TCGA Annotations. First and foremost, an annotation is associated with an itemType, as described above. This can be a patient, an aliquot, etc. Let's see what the breakdown is of annotations according to item-type:
#' 
## ------------------------------------------------------------------------
querySql <- "
SELECT itemTypeName, COUNT(*) AS n
FROM [isb-cgc:tcga_201510_alpha.Annotations]
GROUP BY itemTypeName
ORDER BY n DESC"

result <- query_exec(querySql, project=project)
head(result)

#' 
#' The length of the barcode in the itemBarcode field will depend on the value in the itemTypeName field: if the itemType is 'Patient', then the barcode will be something like TCGA-E2-A15J, whereas if the itemType is "Aliquot", the barcode will be a full-length barcode, eg TCGA-E2-A15J-10A-01D-a12N-01.
#' 
#' ## Annotation Classifications and Categories
#' The next most important pieces of information about an annotation are the "classification" and "category". Each of these comes from a controlled vocabulary and each "classification" has a specific set of allowed "categories".
#' One important thing to understand is that if an aliquot carries some sort of disqualifying annotation, in general all other data from other samples or aliquots associated with that same patient should still be usable. On the other hand, if a patient carries some sort of disqualifying annotation, then that information should be considered prior to using any of the samples or aliquots derived from that patient.
#' 
#' To illustrate this, let's look at the most frequent annotation classifications and categories when the itemType is Patient:
#' 
## ------------------------------------------------------------------------
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

#' 
#' The results of the previous query indicate that the majority of patient-level annotations are "Notifications", most frequently regarding prior malignancies. In most TCGA publications, "history of unacceptable prior treatment" and "item is noncanonical" notifications are treated as disqualifying annotations, and all data associated with those patients is not used in any analysis.
#' Let's make a slight modification to the last query to see what types of annotation categories and classifications we see when the item type is not patient:
#' 
## ------------------------------------------------------------------------
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

#' 
#' The results of the previous query indicate that the vast majority of annotations are at the aliquot level, and more specifically were submitted by one of the data-generating centers, indicating that the data derived from that aliquot is "DNU" (Do Not Use). In general, this should not affect any other aliquots derived from the same sample or any other samples derived from the same patient.
#' We see in the output of the previous query that a Notification that an "Item is noncanonical" can be applied to different types of items (eg slides and analytes). Let's investigate this a little bit further, for example let's count up these types of annotations by study (ie tumor-type):
#' 
## ------------------------------------------------------------------------
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

#' 
#' and now let's pick one of these tumor types, and delve a little bit further:
#' 
## ------------------------------------------------------------------------
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

#' 
#' ## Barcodes
#' As described above, an annotation is specific to a single TCGA "item" and the fields itemTypeName and itemBarcode are the most important keys to understanding which TCGA item carries the annotation. Because we use the fields ParticipantBarcode, SampleBarcode, and AliquotBarcode throughout our other TCGA BigQuery tables, we have added them to this table as well, but they should be interpreted with some care: when an annotation is specific to an aliquot (ie itemTypeName="Aliquot"), the ParticipantBarcode, SampleBarcode, and AliquotBarcode fields will all be set, but this should not be interpreted to mean that the annotation applies to all data derived from that patient.
#' 
#' This will be illustrated with the following two queries which extract information pertaining to a few specific patients:
#' 
## ------------------------------------------------------------------------
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

#' 
#' 
## ------------------------------------------------------------------------
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

#' 
#' As you can see in the results returned from the previous two queries, the SampleBarcode and the AliquotBarcode fields may or may not be filled in, depending on the itemTypeName.
#' 
## ------------------------------------------------------------------------
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

#' 
#' In this example, there is just one annotation relevant to this particular patient, and one has to look at the annotationNoteText to find out what the potential issue may be with this particular analyte. Any aliquots derived from this blood-normal analyte might need to be used with care.
