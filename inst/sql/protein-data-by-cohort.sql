# Retrieve protein data for a particular gene for samples within a cohort.
SELECT
  SampleBarcode,
  Gene_Name,
  Protein_Name,
  protein_expression
FROM
  [_PROTEIN_TABLE_]
WHERE
  Gene_Name = '_GENE_'
  AND Protein_Name = '_PROTEIN_'
  AND SampleBarcode IN (
  SELECT
    sample_barcode
  FROM
    [_COHORT_TABLE_] )
ORDER BY
  SampleBarcode
