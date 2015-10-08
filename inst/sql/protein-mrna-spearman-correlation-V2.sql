# Correlate the protein quantification data with the mRNA expression data.
SELECT
  gene,
  protein,
  COUNT(1) AS num_observations,
  CORR(expr_rank, prot_rank) AS spearman_corr
FROM (
  SELECT
    barcode,
    gene,
    protein,
    RANK() OVER (PARTITION BY gene ORDER BY log2_count ASC) AS expr_rank,
    RANK() OVER (PARTITION BY protein ORDER BY protein_expression ASC) AS prot_rank,
  FROM (
    SELECT
      feat1.SampleBarcode AS barcode,
      Gene_Name AS gene,
      Protein_Name AS protein,
      protein_expression,
      log2_count,
    FROM (
      SELECT
        SampleBarcode,
        Gene_Name,
        Protein_Name,
        protein_expression
      FROM
        [isb-cgc:tcga_data_open.Protein]
      WHERE
        SampleBarcode IN (
        SELECT
          sample_barcode
        FROM
          [isb-cgc:test.cohort_14jun2015] ) ) feat1
    JOIN EACH (
      SELECT
        SampleBarcode,
        HGNC_gene_symbol,
        LOG2(normalized_count+1) AS log2_count
      FROM
        [isb-cgc:tcga_data_open.mRNA_UNC_HiSeq_RSEM]
      WHERE
        SampleBarcode IN (
        SELECT
          sample_barcode
        FROM
          [isb-cgc:test.cohort_14jun2015] ) ) feat2
    ON
      feat1.SampleBarcode = feat2.SampleBarcode
      AND feat1.Gene_Name = feat2.HGNC_gene_symbol))
GROUP BY
  gene,
  protein
HAVING
  num_observations >= _MINIMUM_NUMBER_OF_OBSERVATIONS_
ORDER BY
  spearman_corr DESC
