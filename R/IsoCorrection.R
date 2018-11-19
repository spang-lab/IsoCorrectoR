#' Algorithm For Natural Isotope Abundance And Tracer Purity Correction Of Data From Stable Isotope Labeling Experiments
#'
#' \code{IsoCorrection} is the main function of the IsoCorrectoR package.
#' It performs the correction of mass spectrometry data from stable isotope
#' labeling experiments with regard to natural abundance and tracer purity.
#' Data from both MS and MS/MS experiments can be corrected
#' (with any tracer isotope: 13C, 15N, 18O...), as well as high resolution data
#' from multiple-tracer experiments (e.g. 13C and 15N used simultaneously).
#'
#' @param MeasurementFile Required. The file that contains the measured data to
#'   be corrected. Only ".xls",".xlsx" and ".csv" file formats are supported.
#'
#' @param ElementFile Required. The file that contains the element information
#'   required for correction. Only ".xls",".xlsx" and ".csv" file formats are
#'   supported.
#'
#' @param MoleculeFile Required. The file that contains the information on the
#'   molecules for which data is to be corrected. Only ".xls",".xlsx" and ".csv"
#'   file formats are supported.
#'
#' @param CorrectTracerImpurity Logical. If TRUE, correction for isotopic
#'   impurity of the tracer substrate is performed.
#'
#' @param CorrectTracerElementCore Logical. If TRUE (recommended!), the tracer
#'   element atoms in the core module (usually the part of the molecule that
#'   does not come from derivatization) are considered when correcting.
#'
#' @param CalculateMeanEnrichment Logical. If TRUE, the mean isotopic enrichment
#'   is calculated for each molecule.
#'
#' @param UltraHighRes Logical. If TRUE, high resolution correction is performed
#'   on the data. Should only be set to TRUE, if you know that you have high
#'   resolution data.
#'
#' @param DirOut Character String. Defines the directory the corrected data and
#'   log-file should be written to. Default directory is set to current working
#'   directoy ('.').
#'
#' @param FileOut Character String. Defines the name of the file that contains
#'   the corrected data. The name of the fill will be
#'   IsoCorrectoR_<FileOut>.<FileFormat>. If the format is set to "csv", the
#'   name will also contain the type of the corrected data in the respective
#'   file.
#'
#' @param FileOutFormat Character String. Defines the format of the files that
#'   contain the corrected data. Can either be "csv" or "xls". If set to "csv",
#'   multiple files will be generated, one for each type of corrected data (eg.
#'   corrected data, fractions, mean enrichment...). If set to "xls", all
#'   correction results are provided in one excel file in different sheets.
#'
#' @param ReturnResultsObject Logical. If TRUE, the correction results are
#'   returned as a list in the current R_session in addition to writing the
#'   results to a file. This is useful if the corrected data has to be further
#'   processed directly in R.
#'
#' @param CorrectAlsoMonoisotopic Logical. If TRUE, monoisotopic correction
#'   results are also provided.
#'
#' @param CalculationThreshold (Advanced Option) Numeric. Defines a threshold to
#'   stop probability calculations at for making correction faster (normal
#'   resolution mode). Should be left at the default value.
#'
#' @param CalculationThreshold_UHR (Advanced Option) Numeric. Defines a
#'   threshold to stop probability calculations at for making correction faster
#'   (high resolution mode). Should be left at the default value.
#'
#' @param verbose Logical. If TRUE, status messages are sent to standard output.
#'
#' @param Testmode Logical. If TRUE, starts a testmode for development purposes.
#'   Not required for users of IsoCorrectoR.
#'
#'
#' @import quadprog utils
#'
#' @importFrom WriteXLS WriteXLS
#' @importFrom stringr str_detect
#' @importFrom utils write.csv
#' @importFrom utils packageVersion
#'
#' @export
#'
#' @examples
#'
#' # Normal resolution data
#'
#'    # 1) get path of example files
#'    path.molecule <- system.file("extdata","normal_resolution","MoleculeFile.csv",
#'       package = "IsoCorrectoR", mustWork = TRUE);
#'    path.element <- system.file("extdata","normal_resolution","ElementFile.csv",
#'       package = "IsoCorrectoR", mustWork = TRUE);
#'    path.measurement <- system.file("extdata","normal_resolution","MeasurementFile.csv",
#'       package = "IsoCorrectoR", mustWork = TRUE);
#'
#'    # 2) run correction algorithm and save results in variable
#'    correctionResults <- IsoCorrection(MeasurementFile=path.measurement,
#'       ElementFile=path.element,
#'      MoleculeFile=path.molecule)
#'
#' # High resolution data
#'
#'    # 1) get path of example files
#'    path.molecule <- system.file("extdata","high_resolution","MoleculeFile.csv",
#'       package = "IsoCorrectoR", mustWork = TRUE);
#'    path.element <- system.file("extdata","high_resolution","ElementFile.csv",
#'       package = "IsoCorrectoR", mustWork = TRUE);
#'    path.measurement <- system.file("extdata","high_resolution","MeasurementFile.csv",
#'       package = "IsoCorrectoR", mustWork = TRUE);
#'
#'    # 2) run correction algorithm and save results in variable
#'    correctionResults <- IsoCorrection(MeasurementFile=path.measurement,
#'       ElementFile=path.element,
#'      MoleculeFile=path.molecule,UltraHighRes=TRUE)
#'
#' @return The function returns a list with 4 elements
#'
#' \describe{
#'  \item{success: }{string that is "TRUE" if the correction was successful, "FALSE" if an error occured and "WARNINGS"
#'  if warnings occured}
#'  \item{results: }{a list containing a dataframe for each type of corrected data (normal, fractions, mean enrichment ...).
#'  Will be NA if ReturnResultsObject is set to FALSE}
#'  \item{log: }{list containing log information on the correction run (parameters, file names and paths, warnings and errors)}
#'  \item{error: }{contains a string with the associated error message if an error occurred, empty otherwise}
#' }
#'
#' @references See Reference 1 \url{Link to IsoCorrectoR-Paper}
#'

IsoCorrection <-
  function(MeasurementFile = NA,
             ElementFile = NA,
             MoleculeFile = NA,
             CorrectTracerImpurity = FALSE,
             CorrectTracerElementCore = TRUE,
             CalculateMeanEnrichment = TRUE,
             UltraHighRes = FALSE,
             DirOut = ".",
             FileOut = "result",
             FileOutFormat = "csv",
             ReturnResultsObject = TRUE,
             CorrectAlsoMonoisotopic = FALSE,
             CalculationThreshold = 10^-8,
             CalculationThreshold_UHR = 8,
             verbose=FALSE,
             Testmode = FALSE) {
    message("\n", rep("#", options()$width))
    message(date(), " :: Starting IsoCorrection ...\n")
    # message(rep("#", options()$width))

    version <- utils::packageVersion("IsoCorrectoR")

    timestamp <- Sys.time()
    timestampFMT <- strftime(timestamp, "%Y-%m-%d_%H%M%S")

    correctionOutput <- list()

    correctionOutput$success <- "FALSE"
    correctionOutput$results <- NA
    correctionOutput$log <- NA
    correctionOutput$error <- NA

    roundDigits <- 8

    # DEFINE LOG ENVIRONMENT

    log.env <- new.env()

    log.env$error <- character()
    log.env$warning <- character()
    log.env$general <- character()

    log.env$param <- list()

    log.env$param$Timestamp <- timestamp
    log.env$param$TimestampFMT <- timestampFMT

    log.env$param$version <- version
    log.env$param$date <- date()
    log.env$param$MeasurementFile <- MeasurementFile
    log.env$param$MoleculeFile <- MoleculeFile
    log.env$param$ElementFile <- ElementFile
    log.env$param$NamePrefix <- "IsoCorrectoR"
    log.env$param$WorkingDirectory <- getwd()
    log.env$param$LogfileName <- paste0(log.env$param$NamePrefix, ".log")
    log.env$param$FileOut <- FileOut
    log.env$param$FileOutFormat <- FileOutFormat

    log.env$param$CorrectTracerImpurity <- CorrectTracerImpurity
    log.env$param$CorrectTracerElementCore <- CorrectTracerElementCore
    log.env$param$UltraHighRes <- UltraHighRes
    log.env$param$CalculateMeanEnrichment <- CalculateMeanEnrichment

    if (!UltraHighRes) {
      log.env$param$threshold <- CalculationThreshold
    }
    else {
      log.env$param$threshold <- CalculationThreshold_UHR
    }

    tryCatch({

      # DEFINE OUTPUT DIRECTORY

      if (DirOut == "." || DirOut == "") {
        log.env$param$OutputDirectory <- log.env$param$WorkingDirectory
      }

      else {
        log.env$param$OutputDirectory <- DirOut
      }

      if (Testmode == FALSE) {
        log.env$param$OutputDirectory <-
          file.path(log.env$param$OutputDirectory, timestampFMT)
      }

      ifelse(
        !dir.exists(log.env$param$OutputDirectory),
        dir.create(log.env$param$OutputDirectory, recursive = TRUE),
        FALSE
      )

      # CHECK FUNCTION PARAMETERS

      checkIsoCorrectionParameters(
        MeasurementFile = MeasurementFile,
        ElementFile = ElementFile,
        MoleculeFile = MoleculeFile,
        CorrectTracerImpurity = CorrectTracerImpurity,
        CorrectTracerElementCore = CorrectTracerElementCore,
        CalculateMeanEnrichment = CalculateMeanEnrichment,
        UltraHighRes = UltraHighRes,
        FileOut = FileOut,
        FileOutFormat = FileOutFormat,
        ReturnResultsObject = ReturnResultsObject,
        CorrectAlsoMonoisotopic = CorrectAlsoMonoisotopic,
        CalculationThreshold = CalculationThreshold,
        CalculationThreshold_UHR = CalculationThreshold_UHR,
        logEnvironment = log.env,
        Testmode = Testmode
      )

      # DEFINE NAMES OF OUTPUT FILES

      NamesCorrectedOutput <- defineOutputFilenames(
        logEnvironment = log.env, FileOut = FileOut,
        FileOutFormat = FileOutFormat,
        CalculateMeanEnrichment = CalculateMeanEnrichment,
        CorrectAlsoMonoisotopic = CorrectAlsoMonoisotopic
      )

      # LOAD INPUT FILES, CHECK THEIR STRUCTURE, EXTRACT INFORMATION AND CHECK 
      # INPUT LOGIC

      InputFileInformation <- ExtractInputFileInformation(
        ElementFile = ElementFile, MoleculeFile = MoleculeFile,
        MeasurementFile = MeasurementFile, UltraHighRes = UltraHighRes,
        CorrectTracerImpurity = CorrectTracerImpurity,
        logEnvironment = log.env,
        verbose = verbose
      )

      MoleculeInfo <- InputFileInformation$MoleculeInfo
      ElementInfo <- InputFileInformation$ElementInfo
      MeasurementInfo <- InputFileInformation$MeasurementInfo

      MoleculesTotal <- length(MoleculeInfo)
      SamplesTotal <- ncol(MeasurementInfo)

      # CALCULATION OF PROBABILITY MATRIX

      # In this section, the probability matrix used for data correction 
      # is produced.

      if (verbose) {
        message(
          date(),
          " :: preparing ",
          MoleculesTotal,
          " molecules for correction ..."
        )
      }

      # if(verbose){message(date(), " :: UltraHighRes==", UltraHighRes)}
      if (verbose) {
        if (UltraHighRes == FALSE) {
          message(date(), " :: running in normal resolution mode.")
        } else {
          message(date(), " :: running in ultra high resolution mode,")
        }
      }


      # NORMAL RESOLUTION

      if (!UltraHighRes) {
        CombinedProbList <- ProbNormalRes(
          MoleculeInfo = MoleculeInfo,
          MoleculesTotal = MoleculesTotal,
          ElementInfo = ElementInfo,
          CorrectTracerElementCore = CorrectTracerElementCore,
          CorrectTracerImpurity = CorrectTracerImpurity,
          CalculationThreshold = CalculationThreshold,
          verbose = verbose
        )

        # HIGH RESOLUTION
      } else if (UltraHighRes) {
        CombinedProbList <- ProbUltraHighRes(
          MoleculeInfo = MoleculeInfo,
          MoleculesTotal = MoleculesTotal,
          ElementInfo = ElementInfo,
          CorrectTracerElementCore = CorrectTracerElementCore,
          CorrectTracerImpurity = CorrectTracerImpurity,
          CalculationThreshold_UHR = CalculationThreshold_UHR,
          verbose = verbose
        )
      } # (!UltraHighRes)

      # DATA CORRECTION USING THE PROBABILITY MATRIX

      CorrectedData <- RawDataCorrectionOnCompleteDataset(
        MeasurementInfo = MeasurementInfo,
        MoleculeInfo = MoleculeInfo,
        CombinedProbList = CombinedProbList,
        MoleculesTotal = MoleculesTotal,
        SamplesTotal = SamplesTotal,
        NamesCorrectedOutput = NamesCorrectedOutput,
        roundDigits = roundDigits,
        logEnvironment = log.env,
        verbose = verbose
      )

      CorrectedDataOutput <- CorrectedData$CorrectedDataOutput
      CorrectionResultList <- CorrectedData$CorrectionResultList

      # CALCULATION OF MEAN ENRICHMENT

      if (CalculateMeanEnrichment) {
        MeanEnrichmentDataOutput <-
          CalcMeanEnrichment(
            MoleculeInfo = MoleculeInfo,
            MoleculesTotal = MoleculesTotal,
            SamplesTotal = SamplesTotal,
            roundDigits = roundDigits,
            CorrectionResultList = CorrectionResultList,
            UltraHighRes = UltraHighRes,
            verbose = verbose
          )

        CorrectedDataOutput$MeanEnrichment <- MeanEnrichmentDataOutput
      }

      # WRITE CORRECTED DATA TO FILE

      writeCorrectionResultsToFile(
        CorrectedDataOutput = CorrectedDataOutput,
        NamesCorrectedOutput = NamesCorrectedOutput,
        logEnvironment = log.env,
        verbose = verbose
      )

      # END OF CORRECTION

      if (ReturnResultsObject) {
        correctionOutput$results <- CorrectedDataOutput
      }
      
      if(length(log.env$warning) == 0) {
      
        correctionOutput$success <- "TRUE"
        
      } else {
        
        correctionOutput$success <- "WARNINGS"
        
      }

      # verbose independent message
      # message("\n\n", rep("#", options()$width))
      message(date(), " :: IsoCorrection successfully completed!")
      message(rep("#", options()$width), "\n")
    },

    # ERROR HANDLING
    error = function(e) {
      commonErrorString <-
        "\n\nTHE CORRECTION PROCEDURE WAS ABORTED BECAUSE AN ERROR HAS OCCURED.
      PLEASE CHECK YOUR LOG-FILE.\n\n"

      log.env$error <- c(log.env$error, e)
      
      # verbose independent message
      
      message(paste0(commonErrorString, e))
        
    }
    )

    tryCatch({
      writeLog(log.env)
    },
    error = function(e) {
      log.env$error <- c(log.env$error, e)
      
      # verbose independent message
      
      message(paste0("THE LOG-FILE COULD NOT BE WRITTEN. ERROR: ", e[[1]]))
      
    }
    )

    correctionOutput$log <- as.list(log.env)

    correctionOutput$error <- log.env$error[1][[1]]

    return(correctionOutput)
  } # IsoCorrection
