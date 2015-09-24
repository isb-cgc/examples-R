# Compute the correlation between expression and methylation data.
SELECT HGNC_gene_symbol, Probe_ID, CORR(normalized_count, Beta_Value) AS correlation,
FROM (
  # We select the sample-barcode, gene-symbol, gene-expression, probe-id, and beta-value
  # from a "JOIN" of the gene expression data and the methylation data.
  SELECT expr.SampleBarcode, HGNC_gene_symbol, normalized_count, Probe_ID, Beta_value
  FROM [_EXPRESSION_TABLE_] AS expr
  JOIN EACH (
    FLATTEN ( (
      # We select the sample-barcode, sample-type, study-name, probe-id, beta-value, and gene-symbol
      # from a "JOIN" of the methylation data and the methylation annotation tables
      # which are joined on the CpG probe id that exists in both tables.
      # ( for speed we are only working with chr9 for now )
      SELECT SampleBarcode, SampleTypeLetterCode, Study, Probe_ID, Beta_Value, UCSC.RefGene_Name
      FROM [_METHYLATION_TABLE_] AS methData
      JOIN EACH [isb-cgc:platform_reference.methylation_annotation] AS methAnnot
      ON methData.Probe_ID = methAnnot.Name
      # We require that the gene-symbol not be null.
      WHERE
        UCSC.RefGene_Name IS NOT null
        # Optionally add clause here to limit the query to a particular
        # sample types and/or studies.
        #_AND_WHERE_
    ), UCSC.RefGene_Name ) ) AS methyl
  ON
    methyl.UCSC.RefGene_Name = expr.HGNC_gene_symbol AND methyl.SampleBarcode = expr.SampleBarcode )
GROUP BY HGNC_gene_symbol, Probe_ID
ORDER BY correlation DESC
