# This query performs a Spearman (ranked) correlation between the protein- and the mRNA-expression data
# partitioned based on gene and protein.  The correlation is done by the BigQuery CORR function, but first
# we need to use the RANK function to turn expression values into ranks.

#### QUESTION: should the COUNT(1) be COUNT(2) instead ???

SELECT
  gene,
  protein,
  COUNT(1) AS num_observations,
  CORR(expr_rank, prot_rank) AS spearman_corr
FROM (
  # in order to do a ranke-based (Spearman) correlation, we need to turn the expression values into ranks
  SELECT
    barcode,
    gene,
    protein,
    RANK() OVER (PARTITION BY gene, protein ORDER BY log2_count ASC) AS expr_rank,
    RANK() OVER (PARTITION BY gene, protein ORDER BY protein_expression ASC) AS prot_rank,
  FROM (
    # here we need the sample identifier, the gene symbol, the protein name, and the two expression values
    # that come out of the JOIN defined below
    SELECT
      feat1.SampleBarcode AS barcode,
      Gene_Name AS gene,
      Protein_Name AS protein,
      protein_expression,
      log2_count,
    FROM (
      # on the left side of the JOIN, we SELECT the sample identifier, the gene symbol, the protein name,
      # and the protein expression value, but only for samples that are in our specified "cohort"
      SELECT
        SampleBarcode,
        Gene_Name,
        Protein_Name,
        protein_expression
      FROM
        [_PROTEIN_TABLE_]
      WHERE
        SampleBarcode IN (
        SELECT
          SampleBarcode
        FROM
          [_COHORT_TABLE_] ) ) feat1
    JOIN EACH (
      # on the right side of the JOIN, we SELECT the sample identifier, the gene symbol, and the gene
      # expression value (to which we apply a LOG2() operation), again only for samples in our "cohort"
      SELECT
        SampleBarcode,
        HGNC_gene_symbol,
        LOG2(normalized_count+1) AS log2_count
      FROM
        [_EXPRESSION_TABLE_]
      WHERE
        SampleBarcode IN (
        SELECT
          SampleBarcode
        FROM
          [_COHORT_TABLE_] ) ) feat2
    ON
      # the JOIN needs to line up rows where the sample identifier and the gene symbol match
      feat1.SampleBarcode = feat2.SampleBarcode
      AND feat1.Gene_Name = feat2.HGNC_gene_symbol))
GROUP BY
  gene,
  protein
HAVING
  # although we will compute correlations for all genes and proteins, we only want to keep values where
  # the correlation was based on at least 30 observations
  num_observations >= _MINIMUM_NUMBER_OF_OBSERVATIONS_
ORDER BY
  # and finally we want to order the ouputs from largest positive correlation, descending to the most negative
  # correlation
  spearman_corr DESC
