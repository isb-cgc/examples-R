# Correlate the protein quantification data with the mRNA expression data.
SELECT
  feat1.gene AS gene,
  CORR(feat1.exp_rank, feat2.exp_rank) AS spearman_corr
FROM (
  SELECT
    *,
    RANK() OVER (PARTITION BY gene ORDER BY protein_expression ASC) AS exp_rank
  FROM (
    SELECT
      SampleBarcode,
      Gene_Name AS gene,
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
    RANK() OVER (PARTITION BY gene ORDER BY log2_normalized_count ASC) AS exp_rank
  FROM (
    SELECT
      SampleBarcode,
      HGNC_gene_symbol AS gene,
      IF(0 = normalized_count, 0, LOG2(normalized_count)) AS log2_normalized_count
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
  gene
ORDER BY
  spearman_corr
