# Compute the correlation between expression and methylation data.

SELECT
  HGNC_gene_symbol,
  Probe_ID,
  COUNT(DISTINCT(expr.SampleBarcode)) AS num_observations,
  CORR(log2_count, Beta_Value) AS correlation,
FROM (
  # We select the sample-barcode, gene-symbol, gene-expression, probe-id, and beta-value
  # from a "JOIN" of the gene expression data and the methylation data.  Note that we log-
  # transform the expression since the value in the table is a normalized_count value. 
  SELECT
    expr.SampleBarcode,
    HGNC_gene_symbol,
    LOG2(normalized_count+1) AS log2_count,
    Probe_ID,
    Beta_value
  FROM
    [_EXPRESSION_TABLE_] AS expr
  JOIN EACH ( FLATTEN ( (
        # We select the sample-barcode, sample-type, study-name, probe-id, beta-value, and gene-symbol
        # from the results of a "JOIN" of the methylation data and the methylation annotation tables
        # which are joined on the CpG probe id that exists in both tables.  Note that we need to 
        # FLATTEN this because the UCSC.RefGene information is a (potentially) repeated field.
        SELECT
          SampleBarcode,
          SampleTypeLetterCode,
          Study,
          Probe_ID,
          Beta_Value,
	        CHR,
          UCSC.RefGene_Name
        FROM
          [_METHYLATION_TABLE_] AS methData
        JOIN EACH [isb-cgc:platform_reference.methylation_annotation] AS methAnnot
        ON
          methData.Probe_ID = methAnnot.Name
          # We require that the gene-symbol not be null.
        WHERE
          UCSC.RefGene_Name IS NOT NULL
          # Optionally add clause here to limit the query to a particular
          # sample types and/or study and/or chromosome.
          _AND_WHERE_
        ), UCSC.RefGene_Name ) ) AS methyl
  ON
    methyl.UCSC.RefGene_Name = expr.HGNC_gene_symbol
    AND methyl.SampleBarcode = expr.SampleBarcode )
GROUP BY
  HGNC_gene_symbol,
  Probe_ID
HAVING
  num_observations >= _MINIMUM_NUMBER_OF_OBSERVATIONS_ AND correlation > -2.
ORDER BY
  correlation ASC
