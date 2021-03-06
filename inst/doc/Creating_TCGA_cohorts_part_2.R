#' # Creating TCGA cohorts (part 2)
#' 
#' This notebook will provide another example of building a cohort, this time based on the types of samples available.
#' 
## ------------------------------------------------------------------------
require(bigrquery) || install.packages("bigrquery")

bigrquery::list_tables("isb-cgc", "tcga_201510_alpha")

#' 
#' Many different types of samples were obtained from the TCGA participants, and details about these samples are available in the Biospecimen data table. This next query shows how many samples exist of each type, as well as the full names and abbreviations of each type:
#' 
## ------------------------------------------------------------------------
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

#' 
#' Note that there are many types of tumor samples: primary, metastatic, recurrent, etc, although the vast majority are samples from primary tumors. In the TCGA project, almost all tumor samples were assayed on multiple platforms for mRNA and miRNA expression, DNA methylation, DNA copy-number, and either exome- or whole-genome DNA sequence. For some tumor samples, protein activity was also measured using RPPA arrays. When available, adjacent "normal" tissue samples were also assayed on a subset of these platforms. The "blood normal" samples were primarily used only as a reference source of germline DNA in order to call somatic mutations.
#' 
#' We can do a similar counting exercise of the sample types represented in one of the molecular data tables, using one of the mRNA expression data tables:
#' 
## ------------------------------------------------------------------------
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

#' 
#' In this example, let's assume that we would like to do a study that requires a primary tumor sample and a matched-normal (adjacent) tissue sample. In order to find out which patients provided which types of samples, we need to query the Biospecimen data table. This next query module uses two sub-queries, one to get all patients with TP samples and another to get all patients with NT samples. The final query joins these two and returns a single list of patients.
#' 
## ------------------------------------------------------------------------
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

#' 
## ------------------------------------------------------------------------
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

#' 
## ------------------------------------------------------------------------
# and finally join these two lists to get the intersection of these two lists
patients_both <- (intersect(tp_result$ParticipantBarcode, nt_result$ParticipantBarcode))
head(patients_both)

#' 
#' It might be interesting to find out what the distribution of tumor types is for this list of patients with matched tumor-normal sample pairs. We can define a new SQL module that refers to the results of a previously defined query as long as we pass that reference in when we call bq.Query():
#' 
## ------------------------------------------------------------------------
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

