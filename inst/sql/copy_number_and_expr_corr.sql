
SELECT
  expr.Study as Study,
  expr.HGNC_gene_symbol as Gene,
  CORR(log2(expr.normalized_count+1), cn.mean_segment_mean) as correlation
FROM
  [_EXPRESSION_TABLE_] AS expr
JOIN (
  SELECT
    Study,
    SampleBarcode,
    Num_Probes,
    AVG(Segment_Mean) AS mean_segment_mean,
  FROM
    [_COPY_NUMBER_TABLE_]
  WHERE
    SampleTypeLetterCode = '_TP_'
    AND Chromosome = '_CHR_'
    AND ((start <= _START_ AND END >= _END_)
      OR (start >= _START_ AND start <= _END_)
      OR (END >= _START_   AND END <= _END_)
      OR (start >= _START_ AND END <= _END_))
  GROUP BY
    Study,
    SampleBarcode,
    Num_Probes ) AS cn
ON
  expr.SampleBarcode = cn.SampleBarcode
  AND expr.Study = cn.Study
WHERE
  expr.HGNC_gene_symbol = '_GENE_'
  AND expr.Study = '_STUDY_'
GROUP BY
  Study,
  Gene
ORDER BY
  correlation
