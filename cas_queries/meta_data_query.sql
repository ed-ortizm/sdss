SELECT
  specobjid,
  mjd, plate, fiberid, run2d,
  ra, dec,
  snMedian,
  z, zErr, zWarning,
  z_noqso, zErr_noqso, zWarning_noqso,
  class, subClass,
  targetType,
  programname,
  instrument
FROM SpecObj
WHERE survey='sdss'
ORDER BY snMedian DESC
