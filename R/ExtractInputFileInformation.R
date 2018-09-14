# Load input files, check their structure, extract information and check input logic

ExtractInputFileInformation <- function(ElementFile, MoleculeFile,
                                        MeasurementFile, UltraHighRes,
                                        CorrectTracerImpurity, logEnvironment,
                                        verbose) {
  ElementData <-
    checkElementDataStructure(
      data = readFileByExt(ElementFile, verbose = verbose), logEnvironment =
        logEnvironment, verbose = verbose
    )
  MoleculeData <-
    checkMoleculeDataStructure(
      data = readFileByExt(MoleculeFile, verbose = verbose),
      logEnvironment = logEnvironment, verbose = verbose
    )
  RawData <-
    checkRawData(
      inputFile = MeasurementFile, logEnvironment = logEnvironment,
      verbose = verbose
    )

  # Extract user input from the element information file

  ElementInfo <-
    ElementInfoExtraction(
      ElementData = ElementData,
      UltraHighRes = UltraHighRes,
      logEnvironment = logEnvironment,
      verbose = verbose
    )

  # Extract user input from the molecule information file

  MoleculeInfo <- MoleculeInfoExtraction(
    MoleculeData = MoleculeData,
    ElementInfo = ElementInfo,
    UltraHighRes = UltraHighRes,
    CorrectTracerImpurity = CorrectTracerImpurity,
    MoleculesFound = RawData$MoleculesFound,
    logEnvironment = logEnvironment,
    verbose = verbose
  )

  # Automatic determination of the transitions/measured isotopologues expected
  # according to the molecule parameters entered in the molecule information file

  if (!UltraHighRes) {
    MoleculeInfo <- CalculateTransitionsNR(
      MoleculeInfo = MoleculeInfo,
      ElementInfo = ElementInfo,
      verbose = verbose
    )
  } else {
    MoleculeInfo <- CalculateTransitionsUHR(
      MoleculeInfo = MoleculeInfo,
      ElementInfo = ElementInfo,
      verbose = verbose
    )
  }

  # Extraction of measurement information from the measurement file provided

  resultRawDataExtractionList <-
    RawDataExtraction(
      data = RawData[[1]],
      MoleculeArray = MoleculeInfo,
      logEnvironment = logEnvironment,
      verbose = verbose
    )

  MoleculeInfo <- resultRawDataExtractionList[["MoleculeInfo"]]

  MeasurementInfo <- resultRawDataExtractionList[["dataRaw"]]

  return(list(
    "MoleculeInfo" = MoleculeInfo,
    "ElementInfo" = ElementInfo,
    "MeasurementInfo" = MeasurementInfo
  ))
}
