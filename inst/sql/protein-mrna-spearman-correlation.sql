# Correlate the protein quantification data with the mRNA expression data.
SELECT
  feat1.gene AS gene,
  feat1.protein AS protein,
  COUNT(1) AS num_observations,
  CORR(feat1.exp_rank, feat2.exp_rank) AS spearman_corr
FROM (
  SELECT
    *,
    RANK() OVER (PARTITION BY protein ORDER BY protein_expression ASC) AS exp_rank
  FROM (
    SELECT
      SampleBarcode,
      Gene_Name AS gene,
      Protein_Name AS protein,
      protein_expression
    FROM
      [_PROTEIN_TABLE_]
    WHERE
      SampleBarcode IN (
      SELECT
        sample_barcode
      FROM
        [_COHORT_TABLE_] ) ) ) feat1
JOIN EACH (
  SELECT
    *,
    RANK() OVER (PARTITION BY gene ORDER BY log2_count ASC) AS exp_rank
  FROM (
    SELECT
      SampleBarcode,
      HGNC_gene_symbol AS gene,
      LOG2(normalized_count+1) AS log2_count
    FROM
      [_EXPRESSION_TABLE_]
    WHERE
      SampleBarcode IN (
      SELECT
        sample_barcode
      FROM
        [_COHORT_TABLE_] ) ) )feat2
ON
  feat1.SampleBarcode = feat2.SampleBarcode
  AND feat1.gene = feat2.gene
GROUP BY
  gene,
  protein
HAVING
  num_observations >= _MINIMUM_NUMBER_OF_OBSERVATIONS_
ORDER BY
  spearman_corr DESC
