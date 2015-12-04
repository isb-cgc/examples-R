#' # UNC HiSeq mRNAseq gene expression (RSEM)
#' The goal of this notebook is to introduce you to the mRNAseq gene expression BigQuery table.
#' This table contains all available TCGA Level-3 gene expression data produced by UNC's RNAseqV2 pipeline using the Illumina HiSeq platform, as of October 2015. (Actual archive dates range from January 2013 to June 2015.) The most recent archive (eg unc.edu_BRCA.IlluminaHiSeq_RNASeqV2.Level_3.1.11.0) for each of the 33 tumor types was downloaded from the DCC, and data extracted from all files matching the pattern %.rsem.genes.normalized_results. Each of these raw “RSEM genes normalized results” files has two columns: gene_id and normalized_count. The gene_id string contains two parts: the gene symbol, and the Entrez gene ID, separated by | eg: TP53|7157. During ETL, the gene_id string is split and the gene symbol is stored in the original_gene_symbol field, and the Entrez gene ID is stored in the gene_id field. In addition, the Entrez ID is used to look up the current HGNC approved gene symbol, which is stored in the HGNC_gene_sybmol field.
#' 
#' In order to work with BigQuery, you need to import the bigrquery library and you need to know the name(s) of the table(s) you are going to be working with:
#' 
## ------------------------------------------------------------------------
require(bigrquery) || install.packages("bigrquery")

rnaTable <- "[isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM]"

#' 
#' From now on, occasionally, we will refer to this table using this variable rnaTable, but we could just as well explicitly give the table name each time.
#' 
#' Let's start by taking a look at the table schema:
#' 
## ------------------------------------------------------------------------
querySql <- paste("SELECT * FROM ",rnaTable," limit 1", sep="")
result <- query_exec(querySql, project=project)
data.frame(Columns=colnames(result))

#' 
#' Now let's count up the number of unique patients, samples and aliquots mentioned in this table. We will do this by defining a very simple parameterized query. (Note that when using a variable for the table name in the FROM clause, you should not also use the square brackets that you usually would if you were specifying the table name as a string.)
#' 
## ------------------------------------------------------------------------
for (x in c("ParticipantBarcode", "SampleBarcode", "AliquotBarcode")) {
  querySql <- paste("SELECT COUNT(DISTINCT(",x,"), 25000) AS n ",
                    "FROM ",rnaTable)
  result <- query_exec(querySql, project=project)
  cat(x, ": ", result[[1]], "\n")
}

#' 
#' We can do the same thing to look at how many unique gene symbols and gene ids exist in the table:
#' 
## ------------------------------------------------------------------------
for (x in c("original_gene_symbol", "HGNC_gene_symbol", "gene_id")) {
  querySql <- paste("SELECT COUNT(DISTINCT(",x,"), 25000) AS n ",
                    "FROM ",rnaTable)
  result <- query_exec(querySql, project=project)
  cat(x, ": ", result[[1]], "\n")
}

#' 
#' Based on the counts, we can see that there are a few instances where the original gene symbol (from the underlying TCGA data file), or the HGNC gene symbol or the gene id (also from the original TCGA data file) is missing, but for the majority of genes, all three values should be available and for the most part the original gene symbol and the HGNC gene symbol that was added during ETL should all match up. This next query will generate the complete list of genes for which none of the identifiers are null, and where the original gene symbol and the HGNC gene symbol match. This list has over 18000 genes in it.
#' 
## ------------------------------------------------------------------------
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
head(result)

#' 
#' We might also want to know how often the gene symbols do not agree:
#' 
## ------------------------------------------------------------------------
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

#' 
#' BigQuery is not just a "look-up" service -- you can also use it to perform calculations. In this next query, we take a look at the mean, standard deviation, and coefficient of variation for the expression of EGFR, within each tumor-type, as well as the number of primary tumor samples that went into each summary statistic.
#' 
## ------------------------------------------------------------------------
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

#' 
#' We can also easily move the gene-symbol out of the WHERE clause and into the SELECT and GROUP BY clauses and have BigQuery do this same calculation over all genes and all tumor types.
#' 
## ------------------------------------------------------------------------
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

subResult <- result[result$exp_mean > 6 & result$n >= 200 & result$exp_cv > 0.5,]

head( subResult[order(subResult$exp_cv, decreasing=T),] )

#' 
#' Since the result of the previous query was quite large (over 600,000 rows representing ~20,000 genes x ~30 tumor types), we could load that table into bigquery for subsequent work. Instead we can explore the results in R.
