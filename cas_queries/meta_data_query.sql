SELECT
  
  specobjid,  
  mjd,  
  plate,  
  fiberid,  
  run2d,  
  ra,  
  dec,  
  z,  
  zErr,  
  zWarning,  
  class,  
  subClass,  
  z_noqso,  
  zErr_noqso,
  zWarning_noqso,
  targetType,
  programname,
  instrument,  
  snMedian
  
FROM SpecObj
WHERE survey='sdss'
ORDER BY snMedian DESC
