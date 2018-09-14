# DATA CORRECTION USING THE PROBABILITY MATRIX
#
# In this function, for each sample and molecule from the the input raw data the
# measured data will be corrected using the probability matrix generated before.

RawDataCorrectionOnCompleteDataset <- function(MeasurementInfo, MoleculeInfo, CombinedProbList, MoleculesTotal, SamplesTotal, NamesCorrectedOutput, roundDigits, logEnvironment, verbose) {
  if (verbose) {
    message(
      date(),
      " :: preparing ",
      MoleculesTotal,
      " molecules for correction [OK]"
    )

    message("\n", rep("#", options()$width))
    message(
      date(),
      " :: Starting data correction for <",
      SamplesTotal,
      " samples> ..."
    )
    message(rep("#", options()$width), "\n")
  }

  CorrectionResultList <- list()

  for (SampleNo in seq_len(SamplesTotal)) {
    if (verbose) {
      message(
        date(),
        " :: processing Sample #",
        SampleNo,
        " [",
        colnames(MeasurementInfo)[SampleNo],
        "] "
      )
    }

    tmpCorrectionResult <- list()

    for (MoleculeNo in seq_len(MoleculesTotal)) {
      ProbMatrix <- CombinedProbList[[MoleculeNo]]

      MoleculeData <- MoleculeInfo[[MoleculeNo]]
      MoleculeName <- names(MoleculeInfo[MoleculeNo])

      TransitionLocation <- unlist(MoleculeData[["TransitionLocation"]])

      # For the molecule defined by MoleculeNo, get the vector of uncorrected values. The location of the value corresponding to a given transition is given by
      # TransitionLocation and based on finding measurement name tags in the input data (RawDataExtraction).

      UncorrectedData <- MeasurementInfo[TransitionLocation, SampleNo]

      # The function 'RawDataCorrection' performs the actual
      # data correction.
      CorrectionResult <- RawDataCorrection(
        UncorrectedData = UncorrectedData,
        MoleculeData = MoleculeData,
        MoleculeName = MoleculeName,
        ProbMatrix = ProbMatrix,
        MoleculeNo = MoleculeNo,
        SampleNo = SampleNo,
        SampleName = colnames(MeasurementInfo)[SampleNo],
        roundDigits = roundDigits,
        logEnvironment = logEnvironment, verbose = verbose
      )

      tmpCorrectionResult[[MoleculeNo]] <- CorrectionResult
    } # MoleculeNo

    names(tmpCorrectionResult) <- names(MoleculeInfo)

    CorrectionResultList[[SampleNo]] <- tmpCorrectionResult
  } # SampleNo

  names(CorrectionResultList) <- colnames(MeasurementInfo)

  x <- convert2matrix(input = CorrectionResultList)

  CorrectedDataOutput <- list()

  CorrectedDataOutput$RawData <- as.data.frame(MeasurementInfo)

  for (name in NamesCorrectedOutput$CorrectedDataNamesBase) {
    CorrectedDataOutput[[name]] <- as.data.frame(x[[name]])
  }

  return(list("CorrectedDataOutput" = CorrectedDataOutput, "CorrectionResultList" = CorrectionResultList))
}
