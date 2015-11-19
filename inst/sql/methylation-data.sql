# Retrieve methylation data for a particular probe id.
SELECT
  SampleBarcode,
  SampleTypeLetterCode,
  Study,
  Probe_ID,
  Beta_Value
FROM [_METHYLATION_TABLE_]
WHERE
  Probe_ID = '_PROBE_'
  _AND_WHERE_
ORDER BY
  SampleBarcode,
  SampleTypeLetterCode,
  Study
