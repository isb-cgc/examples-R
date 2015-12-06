#' # DNA Methylation (JHU-USC beta values)
#' The goal of this notebook is to introduce you to the DNA methylation BigQuery table.
#' 
#' This table contains all available TCGA Level-3 DNA methylation data produced by the JHU-USC methylation pipeline using the Illumina Infinium Human Methylation 27k and 450k platforms, as of October 2015. (Actual archive dates range from October 2012 to August 2015.) The most recent archives (eg jhu-usc.edu_HNSC.HumanMethylation450.Level_3.18.8.0) for each of the 33 tumor types were downloaded from the DCC, and data extracted from all files matching the pattern jhu-usc.edu_%.HumanMethylation%.lvl-3.%.txt. Each of these text files has five columns. The first two columns contain the CpG probe id and the methylation beta value. The additional columns contain annotation information (gene symbol(s), and chromosome and genomic coordinate for the CpG probe). Only the CpG probe id and the beta value were extracted during ETL and stored in this BigQuery table, along with the aliquot ID (which can be found both in the text filename, and in the SDRF file in the mage-tab archive).
#' 
#' WARNING: This BigQuery table contains almost 4 billion rows of data and is over 400 GB in size. When experimenting with new queries, be sure to put a "LIMIT" on the results to avoid accidentally launching a query that might either take a very very long time or produce a very large results table!
#' 
#' In order to work with BigQuery, you need to import the bigrquery library and you need to know the name(s) of the table(s) you are going to be working with:
#' 
## ------------------------------------------------------------------------
methTable <- "[isb-cgc:tcga_201510_alpha.DNA_Methylation_betas]"

#' 
#' From now on, we will refer to this table using this variable methTable, but we could just as well explicitly give the table name each time.
#' Let's start by taking a look at the table schema:
#' 
## ------------------------------------------------------------------------
querySql <- paste("SELECT * FROM ",methTable," limit 1", sep="")
result <- query_exec(querySql, project=project)
data.frame(Columns=colnames(result))

#' 
#' Let's count up the number of unique patients, samples and aliquots mentioned in this table. Using the same approach, we can count up the number of unique CpG probes. We will do this by defining a very simple parameterized query. (Note that when using a variable for the table name in the FROM clause, you should not also use the square brackets that you usually would if you were specifying the table name as a string.)
#' 
#' 
## ------------------------------------------------------------------------
for (x in c("ParticipantBarcode", "SampleBarcode", "AliquotBarcode", "Probe_Id")) {
  querySql <- paste("SELECT COUNT(DISTINCT(",x,"), 25000) AS n ",
                    "FROM ",methTable)
  result <- query_exec(querySql, project=project)
  cat(x, ": ", result[[1]], "\n")
}

#' 
#' As mentioned above, two different platforms were used to measure DNA methylation. The annotations from Illumina are also available in a BigQuery table:
#' 
## ------------------------------------------------------------------------
methRef <- '[isb-cgc:platform_reference.methylation_annotation]'
querySql <- paste("SELECT * FROM ",methRef," limit 1", sep="")
result <- query_exec(querySql, project=project)
data.frame(Columns=colnames(result))

#' 
#' Given the coordinates for a gene of interest, we can find the associated methylation probes.
#' 
## ------------------------------------------------------------------------
geneChr   = "\'3\'"
geneStart = 37034841 - 2500
geneStop  = 37092337 + 2500

querySql <- paste("
SELECT
  IlmnID, Methyl27_Loci, CHR, MAPINFO
FROM
  ",methRef,"
WHERE
  ( CHR=",geneChr,"
    AND ( MAPINFO> ", geneStart," AND MAPINFO< ", geneStop, ") )
ORDER BY
  Methyl27_Loci DESC,
  MAPINFO ASC", sep="")

mlh1Probes = query_exec(querySql, project=project)
mlh1Probes

#' 
#' There are a total of 50 methlyation probes in and near the MLH1 gene, although only 6 of them are on both the 27k and the 450k versions of the platform.
#' 
## ------------------------------------------------------------------------
#write.table(mlh1Probes, sep=",", quote=F, row.names=F, file="mlh1Probes.csv")

#' 
#' Then let's use this data to create a new table. To do that, I used my browser to surf over to the bigquery web interface. There, under my project ID, I used the menu to create a new data set. I named the dataset DG after my initials. Then, within that dataset, I used the menu to create a new table using the csv file. After being walked through the process (which includes identifying each column), the table is uploaded and we can query it. I also took the opportunity, while uploading the table, to rename some of the columns (IlmnID -> ProbeID and MAPINFO -> Pos).
#' 
## ----dna_meth_fig1-------------------------------------------------------
buildQuery <- function(sampleType) {
paste("SELECT
  cpg.ProbeID AS Probe_Id,
  cpg.Methyl27 AS Methyl27_Loci,
  cpg.Chr AS Chr,
  cpg.Pos AS Position,
  data.beta_stdev AS beta_stdev,
  data.beta_mean AS beta_mean,
  data.beta_min AS beta_min,
  data.beta_max AS beta_max
FROM (
  SELECT *
  FROM [DG.mlh1Probes]
) AS cpg
JOIN (
  SELECT
    Probe_Id,
    STDDEV(beta_value) beta_stdev,
    AVG(beta_value) beta_mean,
    MIN(beta_value) beta_min,
    MAX(beta_value) beta_max
    FROM [isb-cgc:tcga_201510_alpha.DNA_Methylation_betas]
    WHERE ( SampleTypeLetterCode= \'",sampleType, "\' )
    GROUP BY Probe_Id
) AS data
ON
  cpg.ProbeID = data.Probe_Id
ORDER BY
  Position ASC", sep="")
}

rTP = query_exec(buildQuery("TP"), project=project) # tumor samples
rNT = query_exec(buildQuery("NT"), project=project) # normal tissue samples

rTP$Type <- "Tumor"
rNT$Type <- "Normal"
rdat <- rbind(rTP, rNT)
ggplot(rdat, aes(beta_mean, color=Type))+ geom_freqpoly(aes(group = Type))
ggplot(rdat, aes(beta_stdev, color=Type))+ geom_freqpoly(aes(group = Type))

#' 
#' From this figure, we can see that, with the exception of the CpG probes near the 3' end of MLH1, the primary tumor samples have a slightly higher average methylation, with significantly greater variability.
