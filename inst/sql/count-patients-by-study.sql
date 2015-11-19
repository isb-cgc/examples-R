# all of the TCGA molecular data tables contain the fields ParticipantBarcode and Study

SELECT Study, COUNT(*) AS n
FROM (
    SELECT ParticipantBarcode, Study
    FROM [_THE_TABLE_]
    GROUP BY ParticipantBarcode, Study
) GROUP BY Study
ORDER BY n DESC
