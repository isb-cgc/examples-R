# Retrieve expression data for a particular gene.
SELECT
  SampleBarcode,
  HGNC_gene_symbol,
  normalized_count
FROM [_EXPRESSION_TABLE_]
WHERE
  HGNC_gene_symbol = '_GENE_'
