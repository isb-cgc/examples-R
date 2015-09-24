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
  # Optionally add clause here to limit the query to a particular
  # sample types and/or studies.
  #_AND_WHERE_
