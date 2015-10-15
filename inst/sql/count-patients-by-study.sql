# all of the TCGA molecular data tables contain the fields ParticipantBarcode and Disease_Code
# TODO: this is actually being re-standardized to be "Study" FIXME
SELECT Disease_Code, COUNT(*) AS n
FROM (
    SELECT ParticipantBarcode, Disease_Code
    FROM [_THE_TABLE_]
    GROUP BY ParticipantBarcode, Disease_Code
) GROUP BY Disease_Code
ORDER BY n DESC
