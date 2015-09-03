# Compute the correlation between expression and methylation data on a particular chromosome for a particular tumor type.
SELECT HGNC_gene_symbol, Probe_ID, CORR(normalized_count, Beta_Value) AS correlation,
FROM (
  # we select the sample-barcode, gene-symbol, gene-expression, probe-id, and beta-value
  # from a "JOIN" of the gene expression data and the methylation data
  SELECT expr.SampleBarcode, HGNC_gene_symbol, normalized_count, Probe_ID, Beta_value
  FROM [isb-cgc:tcga_data_open.mRNA_UNC_HiSeq_RSEM] AS expr
  JOIN EACH (
    FLATTEN ( (
      # we select the sample-barcode, sample-type, study-name, probe-id, beta-value, and gene-symbol
      # from a "JOIN" of the methylation data and the methylation annotation tables
      # which are joined on the CpG probe id that exists in both tables
      # ( for speed we are only working with chr9 for now )
      SELECT SampleBarcode, SampleTypeLetterCode, Study, Probe_ID, Beta_Value, UCSC.RefGene_Name
      FROM [isb-cgc:tcga_data_open.Methylation__CHR_] AS methData
      JOIN EACH [isb-cgc:platform_reference.methylation_annotation] AS methAnnot
      ON methData.Probe_ID = methAnnot.Name
      # we require that the gene-symbol not be null, and we are choosing only CESC TP samples for now
      WHERE ( UCSC.RefGene_Name IS NOT null
          AND SampleTypeLetterCode = 'TP'
          AND Study = '_TUMOR_' )
    ), UCSC.RefGene_Name ) ) AS methyl
  ON
    methyl.UCSC.RefGene_Name = expr.HGNC_gene_symbol AND methyl.SampleBarcode = expr.SampleBarcode )
GROUP BY HGNC_gene_symbol, Probe_ID
ORDER BY correlation DESC
