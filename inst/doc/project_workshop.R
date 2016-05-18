#' ---
#' title: "Differential gene expression associated with HPV integration"
#' output: html_document
#' ---
#' 
#' In this example, we will study the effect of HPV integration on the expression of recurrent target genes in CESC and HNSC tumors.
#' This example demonstrates using R to issue BigQuery queries involving multiple tables across multiple data sets. We will also show users
#' how to bring in their own data to use in conjunction with the TCGA data already available as BigQuery tables. In this exercise,
#' we will reproduce some figures from Tang et. al. [1] to visualize altered expression of host genes frequently targeted by HPV.
#' 
#' References:
#' 1. Tang et. al. The landscape of viral expression and host gene fusion and adaptation in human cancer. Nature Communications 4, Article number:2513|doi:10.1038/ncomms3513
#' 
#' -------------------------------------------------------------------------
#' 
#' ## Analysis workflow
#' 
#' Let's get started then! Let's first load all the required libraries and initialize all global variables
## ------------------------------------------------------------------------
require(bigrquery,quietly = TRUE) || install.packages('bigrquery',verbose = FALSE)
require(tidyr,quietly = TRUE) || install.packages('tidyr',verbose = FALSE)
require(dplyr,quietly = TRUE) || install.packages('dplyr',verbose = FALSE)
library(ggplot2,quietly = TRUE) || install.packages('ggplot2',verbose = FALSE)
library(broom,quietly = TRUE) || install.packages('broom',verbose = FALSE))

#' 
#' Specify cloud project name(s)
## ------------------------------------------------------------------------
#cloud_project_workshop = "your project"

#' 
#' Specify BigQuery datasets you want to work with
## ------------------------------------------------------------------------
tcga_ds = "tcga_201510_alpha"
workshop_ds = "workspace"

#' 
#' First let's make sure everything is workign and list tables in the TCGA dataset.
#' 
## ------------------------------------------------------------------------
bigrquery::list_tables(cloud_project_main,tcga_ds)

#' 
#' Tables we will be using in this example...
#' 
## ------------------------------------------------------------------------
clinical_table = "[isb-cgc:tcga_201510_alpha.Clinical_data]"
gexp_table     = "[isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM]"
ncomms_gene_table = "[isb-cgc:workshop.ncomms3513_s3]"

#' 
#' In this analysis, we are going to be looking at two studies.
#' 
## ------------------------------------------------------------------------
study=c('CESC','HNSC')

#' 
#' Now, let's gather the relevant data from BigQuery
#' 
#' 1. Get all CESC and HNSC samples and their hpv status from clinical data
## ------------------------------------------------------------------------

# this is an example of building queries programmatically.

sqlQuery = paste("SELECT ParticipantBarcode, Study, hpv_calls, hpv_status ",
                 "FROM ", clinical_table,
                 " WHERE Study in (",paste(shQuote(study),collapse = ','),")",sep="")

sqlQuery

hpv_table = query_exec(sqlQuery,project = cloud_project_workshop)

dim(hpv_table)

head(hpv_table)

# We can do some quality control ...
# Assert that if hpv_calls is NA, hpv_status is Negative
stopifnot((is.na(hpv_table$hpv_calls) && hpv_table$hpv_status=="Negative") || !is.na(hpv_table$hpv_calls))

# Let's explore the cohort
ggplot(data=hpv_table, aes(x=hpv_status, fill=Study)) + geom_bar(stat="count", position=position_dodge())

#' 
#' 2. TCGA data or BBT analysis does not give us the location of HPV integration into host sequences,
#' So we'll get a list of frequently targeted genes published with this paper:
#' Ka-Wei Tang et. al. The Landscape of viral expression and host gene fusion and adaptation in human cancer. doi:10.1038/ncomms3513
#' 
#' (Supplementary Data 2: Integration analysis results)
#' 
#' We will access the data from our cloud bucket by either using the command line or from the browser.
#' 
#' Using the google command line tool:
#' gsutil cp gs://isb-cgc-workshop-data/ncomms3513-s3.tsv .
#' gsutil cp gs://isb-cgc-workshop-data/ncomms3513-s3_Schema.json .
#' 
#' Using the cloud console, go to https://console.cloud.google.com and find the
#' workshop bucket.
#' 
#' Then, to load the data into a BQ table, we use the 'bq' command. Make sure to change the table directory (the 'DG')!
#' bq load --source_format CSV --field_delimiter "\t"  --schema ncomms3513-s3_Schema.json  DG.ncomms3513_s3 ncomms3513-s3.tsv
#' 
#' Now we can directly query 'our' own data, and start to combine it with other tables.
#' 
## ------------------------------------------------------------------------
sqlQuery = "
SELECT
  Overlapping_genes,
  Cancer
FROM
  [isb-cgc-02-0001:DG.ncomms3513_s3]
WHERE
  Cancer IN ('CESC',
    'HNSC')
  AND Overlapping_genes <> 'Intergenic'
GROUP BY
  Cancer,
  Overlapping_genes
  "

affected_genes = query_exec(sqlQuery,project = cloud_project_workshop)

head(affected_genes)

table(affected_genes$Cancer)

#' 3. Now, we want to get gene expression data for affected_genes for the tumor types they are affected in
#' 
## ------------------------------------------------------------------------
sqlQuery = "
SELECT
  ParticipantBarcode,
  SampleBarcode,
  Study,
  HGNC_gene_symbol,
  normalized_count
FROM
  [isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM]
WHERE
  Study IN ('CESC','HNSC')
  AND SampleTypeLetterCode = 'TP'
  AND HGNC_gene_symbol IN (
  SELECT
    Overlapping_genes as HGNC_gene_symbol
  FROM
    [isb-cgc-04-0030:workspace.ncomms3513_s3]
  WHERE
    Cancer IN ('CESC','HNSC')
    AND Overlapping_genes <> 'Intergenic'
  GROUP BY
    HGNC_gene_symbol )
"

gexp_affected_genes = query_exec(sqlQuery,project = cloud_project_workshop)

#view results
head(gexp_affected_genes)

# a couple different ways to look at the results
#qplot(data=gexp_affected_genes, x=Study, y=normalized_count, col=HGNC_gene_symbol, geom="boxplot")
#qplot(data=gexp_affected_genes, x=Study, y=log2(normalized_count), col=HGNC_gene_symbol, geom="boxplot")
qplot(data=gexp_affected_genes, x=log2(normalized_count+1), col=HGNC_gene_symbol, geom="density") + facet_wrap(~ Study)

#' 
#' Not all samples listed in the clinical data have gene expression data.
#' Let's filter the hpv_table to match the samples to those in gexp_affected_genes
#' 
## ------------------------------------------------------------------------
# let's get rid of 'indeterminate' samples
hpv_table = dplyr::filter(hpv_table, hpv_status != "Indeterminate", ParticipantBarcode %in% gexp_affected_genes$ParticipantBarcode)

#' 
#' Then, we are going to use a couple nice libraries to perform t.tests on expression by hpv_status
#' 
## ------------------------------------------------------------------------
gxps <- merge(x=gexp_affected_genes, y=hpv_table, by=c("Study","ParticipantBarcode"))

# Performing a t-test between hpv+ and hpv- by study and gene
res0 <- gxps %>%
group_by(Study, HGNC_gene_symbol) %>%
do(tidy(t.test(log2(normalized_count+1) ~ hpv_status, data=.))) %>%
ungroup() %>%
arrange(desc(statistic))

# These are the top 5 results ...
top5 <- select(top_n(res0, 5, statistic), Study, HGNC_gene_symbol)

# Let's subset the data by the top 5 results...
res1 <- merge(x=top5, y=gxps) %>% mutate( Study_Gene = paste0(Study, "_", HGNC_gene_symbol))

# now we can plot the results...
ggplot(res1, aes(x=Study_Gene, y=log2(normalized_count+1), fill=hpv_status)) + geom_boxplot()

#' 
#' ## Making BigQueries
#' 
#' Now, we previously downloaded data, and performed some work on it. But another way
#' is to perform as much work as possible in the cloud, and use R to visualize summary results.
#' 
#' First we will compute some statistics on gene expression data.
#' 
#' 
#' "SELECT ParticipantBarcode, Study, hpv_calls, hpv_status FROM [isb-cgc:tcga_201510_alpha.Clinical_data] WHERE Study in ('CESC','HNSC')"
#' 
## ------------------------------------------------------------------------
sqlQuery = "
SELECT
  ParticipantBarcode,
  SampleBarcode,
  Study,
  HGNC_gene_symbol,
  normalized_count
FROM
  [isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM]
WHERE
  Study = 'CESC'
  AND SampleTypeLetterCode = 'TP'
  AND ParticipantBarcode IN (
  SELECT
    ParticipantBarcode
  FROM
    [isb-cgc:tcga_201510_alpha.Clinical_data]
  WHERE
    hpv_status = 'Positive' )
  AND HGNC_gene_symbol IN (
  SELECT
    Overlapping_genes AS HGNC_gene_symbol
  FROM
    [isb-cgc-04-0030:workspace.ncomms3513_s3]
  WHERE
    Cancer = 'CESC'
    AND Overlapping_genes <> 'Intergenic'
  GROUP BY
    HGNC_gene_symbol )
"
q1 = query_exec(sqlQuery,project = cloud_project_workshop)
dim(q1)

#' 
#' Now lets make a small change, and get gene expression for subjects that are
#' hpv negative.
#' 
## ------------------------------------------------------------------------
sqlQuery = "
SELECT
  ParticipantBarcode,
  SampleBarcode,
  Study,
  HGNC_gene_symbol,
  normalized_count
FROM
  [isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM]
WHERE
  Study = 'CESC'
  AND SampleTypeLetterCode = 'TP'
  AND ParticipantBarcode IN (
  SELECT
    ParticipantBarcode
  FROM
    [isb-cgc:tcga_201510_alpha.Clinical_data]
  WHERE
    hpv_status = 'Negative' )
  AND HGNC_gene_symbol IN (
  SELECT
    Overlapping_genes AS HGNC_gene_symbol
  FROM
    [isb-cgc-04-0030:workspace.ncomms3513_s3]
  WHERE
    Cancer = 'CESC'
    AND Overlapping_genes <> 'Intergenic'
  GROUP BY
    HGNC_gene_symbol )
"

q2 <- query_exec(sqlQuery,project = cloud_project_workshop)
dim(q2)

#' 
#' Now we will merge the previous two queries.
#' 
## ------------------------------------------------------------------------
q <- "
SELECT
  p.HGNC_gene_symbol AS gene,
  p.study AS study,
  p.x AS x,
  p.sx2 AS sx2,
  p.nx AS nx,
  o.y AS y,
  o.sy2 AS sy2,
  o.ny AS ny,
  (p.x-o.y) / SQRT((p.sx2/p.nx) + (o.sy2/o.ny)) AS T
FROM (
  SELECT
    Study,
    HGNC_gene_symbol,
    AVG(LOG2(normalized_count+1)) AS y,
    POW(STDDEV(LOG2(normalized_count+1)),2) AS sy2,
    COUNT(ParticipantBarcode) AS ny
  FROM
    [isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM]
  WHERE
    Study = 'CESC'
    AND SampleTypeLetterCode = 'TP'
    AND ParticipantBarcode IN (
    SELECT
      ParticipantBarcode
    FROM
      [isb-cgc:tcga_201510_alpha.Clinical_data]
    WHERE
      hpv_status = 'Positive' )
    AND HGNC_gene_symbol IN (
    SELECT
      Overlapping_genes AS HGNC_gene_symbol
    FROM
      [isb-cgc-04-0030:workspace.ncomms3513_s3]
    WHERE
      Cancer = 'CESC'
      AND Overlapping_genes <> 'Intergenic'
    GROUP BY
      HGNC_gene_symbol )
  GROUP BY
    Study,
    HGNC_gene_symbol) AS o
JOIN (
  SELECT
    Study,
    HGNC_gene_symbol,
    AVG(LOG2(normalized_count+1)) AS x,
    POW(STDDEV(LOG2(normalized_count+1)),2) AS sx2,
    COUNT(ParticipantBarcode) AS nx
  FROM
    [isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM]
  WHERE
    Study = 'CESC'
    AND SampleTypeLetterCode = 'TP'
    AND ParticipantBarcode IN (
    SELECT
      ParticipantBarcode
    FROM
      [isb-cgc:tcga_201510_alpha.Clinical_data]
    WHERE
      hpv_status = 'Negative' )
    AND HGNC_gene_symbol IN (
    SELECT
      Overlapping_genes AS HGNC_gene_symbol
    FROM
      [isb-cgc-04-0030:workspace.ncomms3513_s3]
    WHERE
      Cancer = 'CESC'
      AND Overlapping_genes <> 'Intergenic'
    GROUP BY
      HGNC_gene_symbol )
  GROUP BY
    Study,
    HGNC_gene_symbol) AS p
ON
  p.HGNC_gene_symbol = o.HGNC_gene_symbol
  AND p.Study = o.Study
GROUP BY
  gene,
  Study,
  x,
  sx2,
  nx,
  y,
  sy2,
  ny,
  T
 "
 t_test_result <- query_exec(q, project = cloud_project_workshop)
 head(t_test_result)

#' 
#' ## Extras
#' 
#' Transform gexp_affected_genes_df into a gexp-by-samples feature matrix
#' 
## ------------------------------------------------------------------------
gexp_fm = tidyr::spread(gexp_affected_genes,HGNC_gene_symbol,normalized_count)
gexp_fm[1:5,1:5]

#' 
#' 
## ------------------------------------------------------------------------
# ng-chm

#library(NGCHM)
#library(ISBCHM)
#library(magrittr)
# NOTE: ip address of the NGCHM server is hard coded. Make usre it's the correct
#       one by checking the VM instance on the Google Cloud Console.
#chmCreateManagedServer('cloud','ng-chm','104.154.59.99')
#options(cloudproject='isb-cgc')

#chm= exprCHM('GEXP_Hpv Status',study,getStudyCohort(study),affected_genes_df[,],
#'Comparison of mRNA expression levels between HPV positive and HPV negative CESC samples')

#exprCHM() from ISBCHM ends here    
#plot(chm)

