# Retrieve expression data for a particular gene for samples within a cohort.
SELECT
  SampleBarcode,
  HGNC_gene_symbol,
  normalized_count
FROM [_EXPRESSION_TABLE_]
WHERE
  HGNC_gene_symbol = '_GENE_'
  AND SampleBarcode IN (
  SELECT
    sample_barcode
  FROM
    [_COHORT_TABLE_] )
ORDER BY
  SampleBarcode
