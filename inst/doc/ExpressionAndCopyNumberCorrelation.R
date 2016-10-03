#' # Expression and Copy Number Correlation
#' 
#' In this example, we will look at the correlation between mRNAseq-based gene expression and copy number data.  We will do this using two data tables from the isb-cgc:tcga_201510_alpha dataset and a genome table from the isb-cgc:genome_reference dataset.
#' 
#' 
## ----message=FALSE-------------------------------------------------------
library(dplyr)
library(bigrquery)
library(ggplot2)
library(stringr)
library(ISBCGCExamples)

# The directory in which the files containing SQL reside.
#sqlDir = file.path("/PATH/TO/GIT/CLONE/OF/examples-R/inst/",
sqlDir = file.path(system.file(package = "ISBCGCExamples"), "sql")

#' 
## ----eval=FALSE----------------------------------------------------------
## ######################[ TIP ]########################################
## ## Set the Google Cloud Platform project id under which these queries will run.
## ##
## ## If you are using the workshop docker image, this is already
## ## set for you in your .Rprofile and you can skip this step.
## 
## # project = "YOUR-PROJECT-ID"
## #####################################################################

#' 
#' ## Getting gene information
#' 
#' We're going to use the genome data set to query some information about
#' a given gene.
#' 
## ------------------------------------------------------------------------
q <- "select
         seqname,
         MIN(start) as start,
         MAX(end) as end,
         gene_name
      from
         [isb-cgc:genome_reference.GENCODE_v19]
      where
         gene_name = 'TP53'
      group by
         seqname,
         gene_name"

gene_info <- query_exec(q, project)
gene_info

#' 
#' 
#' ## Pearson Correlation in BigQuery
#' 
#' This sql example builds a table where for each sample barcode, we have
#' the RNA-seq level and the mean of copy number changes for segments that
#' overlap the gene region. Pearson correlation is computed over expr and
#' copy number columns.
#' 
## ----comment=NA----------------------------------------------------------
# Set the desired tables to query.
expressionTable = "isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM"
copyNumberTable = "isb-cgc:tcga_201510_alpha.Copy_Number_segments"

# Now we are ready to run the query.
result = DisplayAndDispatchQuery(
     file.path(sqlDir, "copy_number_and_expr_corr.sql"),
               project=project,
               replacements=list("_STUDY_"="BRCA",
                                 "_EXPRESSION_TABLE_"=expressionTable,
                                 "_COPY_NUMBER_TABLE_"=copyNumberTable,
                                 "_GENE_"=gene_info$gene_name,
                                 "_TP_"="TP",
                                 "_CHR_"=str_sub(gene_info$seqname, start=4),
                                 "_START_"=gene_info$start,
                                 "_END_"=gene_info$end))
result

#' 
#' The result is a table with the study, gene, and pearson correlation for
#' expression with mean of copy number segments overlapping the gene of interest.
#' 
#' But now we can run it over all studies.
#' 
## ------------------------------------------------------------------------

# Now we are ready to run the query.
result = DisplayAndDispatchQuery(
     file.path(sqlDir, "copy_number_and_expr_corr_all_studies.sql"),
               project=project,
               replacements=list("_EXPRESSION_TABLE_"=expressionTable,
                                 "_COPY_NUMBER_TABLE_"=copyNumberTable,
                                 "_GENE_"=gene_info$gene_name,
                                 "_TP_"="TP",
                                 "_CHR_"=str_sub(gene_info$seqname, start=4),
                                 "_START_"=gene_info$start,
                                 "_END_"=gene_info$end))
result

#' 
## ----barplot, fig.align="center", fig.width=10, message=FALSE, warning=FALSE, comment=NA----

# Plot comparing the CN ~ Expr correlation across studies.
ggplot(result, aes(x=factor(Study, ordered=T, levels=Study), y=correlation)) +
geom_point() +
theme(axis.text.x=element_text(angle=90, size=8)) +
xlab("Study") + ylab("Pearson correlation of gene expr and copy number data") +
ggtitle(gene_info$gene_name)

#' 
#' Looks like the MESO has one of the strongest correlations. Let's check the data.
#' 
## ------------------------------------------------------------------------
q <- "
SELECT
  expr.Study as Study,
  expr.HGNC_gene_symbol as Gene,
  expr.SampleBarcode as SampleBarcode,
  log2(expr.normalized_count+1) as expr,
  cn.mean_segment_mean as mean_cn
FROM
  [isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM] AS expr
JOIN (
  SELECT
    Study,
    SampleBarcode,
    Num_Probes,
    AVG(Segment_Mean) AS mean_segment_mean,
  FROM
    [isb-cgc:tcga_201510_alpha.Copy_Number_segments]
  WHERE
    SampleTypeLetterCode = 'TP'
    AND Study = 'LIHC'
    AND Chromosome = '17'
    AND ((start <= 7565097 AND END >= 7590856)
      OR (start >= 7565097 AND start <= 7590856)
      OR (END >= 7565097   AND END <= 7590856)
      OR (start >= 7565097 AND END <= 7590856))
  GROUP BY
    Study,
    SampleBarcode,
    Num_Probes ) AS cn
ON
  expr.SampleBarcode = cn.SampleBarcode
  AND expr.Study = cn.Study
WHERE
  expr.HGNC_gene_symbol = 'TP53'
  AND expr.Study = 'LIHC'
GROUP BY
  Study,
  Gene,
  SampleBarcode,
  expr,
  mean_cn
"

data <- query_exec(q, project)
head(data)

# Plot comparing the CN ~ Expr correlation across studies.
ggplot(data, aes(x=expr, y=mean_cn)) +
geom_point() +
geom_smooth(method="lm") +
xlab("Log2 RNA-seq expression level") + ylab("Mean copy number change") +
ggtitle(gene_info$gene_name)


#' 
#' 
#' ## Provenance
## ----provenance, comment=NA----------------------------------------------
sessionInfo()

