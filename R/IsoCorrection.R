#' Algorithm For Natural Isotope Abundance And Tracer Purity Correction Of Data From Stable Isotope Labeling Experiments
#'
#' \code{IsoCorrection} is the main function of the IsoCorrectoR package. 
#' It performs the correction of mass spectrometry data from stable isotope 
#' labeling experiments with regard to natural abundance and tracer purity. 
#' Data from both MS and MS/MS experiments can be corrected 
#' (with any tracer isotope: 13C, 15N, 18O...), as well as high resolution data 
#' from multiple tracer experiments (e.g. 13C and 15N used simultaneously).
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
#' @param Testmode Logical. If TRUE, starts a testmode for development purposes.
#'   Not required for users of IsoCorrectoR.
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
#'  \item{success: }{logical that is TRUE if the correction was successful, FALSE otherwise}
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
           CalculationThreshold = 10 ^ -8,
           CalculationThreshold_UHR = 8,
           Testmode = FALSE) {
    version <- utils::packageVersion("IsoCorrectoR") #CK:TODO
    
    timestamp <- Sys.time()
    timestampFMT <- strftime(timestamp, "%Y-%m-%d_%H%M%S")
    
    correctionOutput <- list()
    
    correctionOutput$success <- FALSE
    correctionOutput$results <- NA
    correctionOutput$log <- NA
    correctionOutput$error <- NA
    
    CalcOnlyProbabilities <- FALSE
    OutputToExcel <- TRUE
    
    roundDigits <- 8
    
    #Define log environment
    
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
    
    tryCatch({
      #Define output directory
      
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
      
      #check other files and parameters
      
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
      
      #Define names of output files
      
      log.env$param$LogfileName <-
        paste0(log.env$param$NamePrefix, "_", FileOut, ".log")
      
      CorrectedDataNamesBase <-
        c("Corrected",
          "CorrectedFractions",
          "Residuals",
          "RelativeResiduals")
      
      if (CorrectAlsoMonoisotopic) {
        CorrectedDataNamesBase <-
          c(
            CorrectedDataNamesBase,
            "CorrectedMonoisotopic",
            "CorrectedMonoisotopicFractions"
          )
        
      }
      
      if (CalculateMeanEnrichment) {
        CorrectedDataNames <-
          c(CorrectedDataNamesBase, "MeanEnrichment", "RawData")
        
      }
      else {
        CorrectedDataNames <- c(CorrectedDataNamesBase, "RawData")
        
      }
      
      if (FileOutFormat == "csv") {
        CorrectedDataFileNames <-
          vector(mode = "character",
                 length = length(CorrectedDataNames))
        names(CorrectedDataFileNames) <- CorrectedDataNames
        
        for (name in CorrectedDataNames) {
          CorrectedDataFileNames[name] <-
            paste0(log.env$param$NamePrefix,
                   "_",
                   FileOut,
                   "_",
                   name,
                   ".",
                   FileOutFormat)
          
        }
        
        log.env$param$FileOut <-
          paste0(CorrectedDataFileNames, collapse = ", ")
        
      }
      else if (FileOutFormat == "xls") {
        log.env$param$FileOut <-
          paste0(log.env$param$NamePrefix,
                 "_",
                 FileOut,
                 ".",
                 FileOutFormat)
        
      }
      
      if (!UltraHighRes) {
        log.env$param$threshold <- CalculationThreshold
        
        Prob_threshold_isotopes = CalculationThreshold
        Prob_threshold_elements = CalculationThreshold
        Prob_threshold_tracer = CalculationThreshold
        Prob_threshold_impurity = CalculationThreshold
        
      }
      else {
        log.env$param$threshold <- CalculationThreshold_UHR
        
        LabelDiffConstraint = CalculationThreshold_UHR
        
      }
      
      ElementData <-
        checkElementDataStructure(data = readFileByExt(ElementFile), logEnvironment =
                                    log.env)
      MoleculeData <-
        checkMoleculeDataStructure(data = readFileByExt(MoleculeFile),
                                   logEnvironment = log.env)
      RawData <-
        checkRawData(inputFile = MeasurementFile, logEnvironment = log.env)
      
      # Extract user input from the element information file
      
      ElementInfo <-
        ElementInfoExtraction(
          ElementData = ElementData,
          UltraHighRes = UltraHighRes,
          logEnvironment = log.env
        )
      
      # Extract user input from the molecule information file
      
      MoleculeInfo <- MoleculeInfoExtraction(
        MoleculeData = MoleculeData,
        ElementArray = ElementInfo,
        UltraHighRes = UltraHighRes,
        CorrectTracerImpurity = CorrectTracerImpurity,
        MoleculesFound = RawData$MoleculesFound,
        logEnvironment = log.env
      )
      
      # Automatic determination of the transitions/measured isotopologues expected
      # according to the molecule parameters entered in the molecule information file
      
      MoleculeInfo2 <- CalculateTransitions(
        MoleculeArray = MoleculeInfo,
        ElementArray = ElementInfo,
        UltraHighRes = UltraHighRes
      )
      
      # Extraction of measurement information from the measurement file provided
      
      resultRawDataExtractionList <-
        RawDataExtraction(
          data = RawData[[1]],
          MoleculeArray = MoleculeInfo2,
          logEnvironment = log.env
        )
      
      MoleculeInfo3 <- resultRawDataExtractionList[["MoleculeInfo"]]
      
      dataRaw <- resultRawDataExtractionList[["dataRaw"]]
      
      MoleculesTotal <- length(MoleculeInfo3)
      
      SamplesTotal <- ncol(dataRaw)
      
      # CALCULATION OF PROBABILITY MATRIX
      #
      # In this section, the probability matrix used for data correction is produced.
      # To this end, for each molecule to be corrected, molecule/fragment
      # information from the molecule information input file is used
      # to calculate the probability that a certain transition or full MS ion x
      # is derived from another transition/full MS ion y due to
      # natural stable isotope abundance/tracer purity.
      
      MissingTransitionCheckList <- list()
      MaxMassShiftList <- list()
      icmResultList <- list()
      ecResultList <- list()
      msprobResultList <- list()
      mstResultList <- list()
      FragmentsCombinedProbabilityList <- list()
      
      if (!UltraHighRes) {
        message("\n",
                date(),
                " :: preparing <",
                MoleculesTotal,
                " molecules> for correction ...")
        for (MoleculeNo in seq_len(MoleculesTotal)) {
          MoleculeData <- MoleculeInfo3[[MoleculeNo]]
          NumberFragments <-
            sum(stringr::str_detect(names(MoleculeData), "Fragment"))
          Transitions <- MoleculeData[["Transitions"]]
          
          tmpIcmResultList <- list()
          tmpEcResultList <- list()
          tmpMsprobResultList <- list()
          tmpMstResultList <- list()
          tmpMaxMassShiftList <- list()
          
          # The following functions are run for both fragments of a molecule
          # in the case of MS/MS measurements or just for one molecule(-fragment) for
          # full-MS measurements.
          
          for (Fragment in seq_len(NumberFragments)) {
            # For each molecule, the expected transitions are checked and MaxMassShift
            # gives the mass shifts of the labelled molecule(-fragment(s)) with the highest mass.
            
            MaxMassShift <- apply(Transitions, 2, max)[Fragment]
            tmpMaxMassShiftList[[Fragment]] <- MaxMassShift
            
            icmResult <-
              IsoCombinationsMaster(
                MoleculeArray = MoleculeInfo3,
                MoleculeNo = MoleculeNo,
                Fragment = Fragment,
                ElementArray = ElementInfo,
                CorrectTracerElementCore = CorrectTracerElementCore,
                CorrectTracerImpurity = CorrectTracerImpurity,
                Prob_threshold_isotopes = Prob_threshold_isotopes
              )
            
            tmpIcmResultList[[Fragment]] <- icmResult
            
            # With the probabilities and mass shifts of the different
            # isotope combinations calculated, the next step is to combine
            # the isotope combination probabilities and mass shifts of
            # the different elements (and the tracer impurity, if checked) using the function
            # 'ElementCombinations'. For a given molecule(-fragment) the function combinatorically
            # multiplies all isotope combination probabilities of all non-tracer elements,
            # the tracer element and the tracer impurity with each other
            # In the same way, the isotope combination mass shifts are added to each
            # other. In the end, the probability that a molecule(-fragment) contains a combination
            # of specific element isotope combinations for its various elements
            # (and a certain number of 'impure' tracer atoms) is yielded,
            # together with the mass shift associated with
            # such a combination. Because combinatorically multiplying all
            # isotope combinations is resource intensive, the
            # calculation stops at certain probability thresholds defined by the
            # parameter CalculationThreshold. In those cases the probability of a combination is so low
            # that it does not affect the correction.
            
            ecResult <-
              ElementCombinations(
                MoleculeArray = MoleculeInfo3,
                MoleculeNo = MoleculeNo,
                Fragment = Fragment,
                IsoCombinationsMaster = icmResult,
                Prob_threshold_elements = Prob_threshold_elements,
                Prob_threshold_tracer = Prob_threshold_tracer,
                Prob_threshold_impurity = Prob_threshold_impurity,
                CorrectTracerImpurity = CorrectTracerImpurity
              )
            
            tmpEcResultList[[Fragment]] <- ecResult
            
            # The function 'MassShiftProbabilities' yields the
            # probability that a given transition of a
            # molecule(-fragment) (through naturally occuring isotopes)
            # produces a certain mass shift in relation to the same molecule(-fragment)
            # containing no isotopes of higher mass. This probability is gained
            # by summing all entries of a molecule(-fragment) of a given
            # transition in TracerElemProb_Array that have the same mass shift in
            # TracerElemMassShift_Array. If correction for tracer impurity is turned
            # on and if the molecule(-fragment) in question can contain tracer,
            # 'MassShiftProbabilities' uses TracerImpurityCombProb_Array
            # and TracerImpurityCombMassShift_Array instead.
            
            CumProbList <-
              MassShiftProbabilities(
                MoleculeArray = MoleculeInfo3,
                MoleculeNo = MoleculeNo
                ,
                Fragment = Fragment
                ,
                ElementCombinations = ecResult,
                MaxMassShift = MaxMassShift
                ,
                CorrectTracerImpurity = CorrectTracerImpurity
              )
            
            tmpMsprobResultList[[Fragment]] <- CumProbList
            
            # The purpose of the function 'MassShiftTransitions' is
            # to generate an array (FragmentMassShift) which contains
            # the mass shift of a molecule(-fragment) of a given
            # transition in relation to the same molecule(-fragment)
            # containing no isotopes of higher mass. This is required to
            # be able to yield the probability that a certain transition
            # or full MS ion x is derived from another transition/full MS ion
            # y due to natural stable isotope abundance. This probability is yielded
            # by matching the mass shift k of transition x in FragmentMassShift
            # with the mass shift k of transition y in CumProbList
            # (this is the mass shift 'produced' by y through natural abundance).
            # The probability is found associated with mass shift k in CumProbList.
            
            mstResult <-
              MassShiftTransitions(
                MoleculeArray = MoleculeInfo3,
                MoleculeNo = MoleculeNo,
                Fragment = Fragment
              )
            
            tmpMstResultList[[Fragment]] <- mstResult
            
          }#Fragment
          names(tmpIcmResultList) <-
            paste0("Fragment_", seq(seq_len(NumberFragments)))
          names(tmpEcResultList) <-
            paste0("Fragment_", seq(seq_len(NumberFragments)))
          names(tmpMsprobResultList) <-
            paste0("Fragment_", seq(seq_len(NumberFragments)))
          names(tmpMstResultList) <-
            paste0("Fragment_", seq(seq_len(NumberFragments)))
          names(tmpMaxMassShiftList) <-
            paste0("Fragment_", seq(seq_len(NumberFragments)))
          
          # In the function 'FragmentsCombinedProbability' the final
          # probability matrix CombFragmentProb is produced. In this matrix, for each
          # transition y of a molecule the probability that it produces a transition x
          # of higher mass through natural isotope abundance is given.
          # This is done by matching the mass shifts from
          # FragmentMassShift and CumProbList, as described before for the
          # function 'MassShiftTransitions'. For MS/MS data, the
          # probability that the product ion y produces the product ion
          # mass shift of transition x and and the probability that neutral
          # loss y produces the neutral loss mass shift of transition x
          # are multiplied. This way the overall probability that
          # transition y produces transition x through natural isotope abundance
          # is yielded.
          
          CombFragmentProb <-
            FragmentsCombinedProbability(
              MoleculeArray = MoleculeInfo3,
              MoleculeNo = MoleculeNo,
              CumProbList = tmpMsprobResultList,
              FragmentMassShift =
                tmpMstResultList
            )
          
          FragmentsCombinedProbabilityList[[MoleculeNo]] <-
            CombFragmentProb
          MaxMassShiftList[[MoleculeNo]] <- tmpMaxMassShiftList
          
          icmResultList[[MoleculeNo]] <- tmpIcmResultList
          ecResultList[[MoleculeNo]] <- tmpEcResultList
          msprobResultList[[MoleculeNo]] <- tmpMsprobResultList
          mstResultList[[MoleculeNo]] <- tmpMstResultList
          
          
        }#MoleculeNo
        names(MaxMassShiftList) <- names(MoleculeInfo3)
        names(icmResultList) <- names(MoleculeInfo3)
        names(ecResultList) <- names(MoleculeInfo3)
        names(msprobResultList) <- names(MoleculeInfo3) # CumProbList
        names(mstResultList) <- names(MoleculeInfo3) # FragmentMassShift
        names(FragmentsCombinedProbabilityList) <-
          names(MoleculeInfo3) # CombFragmentProb
        
      } else if (UltraHighRes) {
        message(date(), " :: UltraHighRes==", UltraHighRes)
        
        CombinedProbList <- list()
        for (MoleculeNo in seq_len(MoleculesTotal)) {
          MoleculeData <- MoleculeInfo3[[MoleculeNo]]
          MoleculeLocation <- MoleculeData[["TransitionLocation"]]
          
          CombinedProb <-
            matrix(
              ncol = nrow(MoleculeData$Transitions),
              nrow = nrow(MoleculeData$Transitions)
            )
          
          # The function ProbUltraHighRes provides the probability array
          # 'CombinedProb' which is then used for data correction of ultra high resolution /
          # multiple tracer data
          
          CombinedProb <- ProbUltraHighRes(
            MoleculeData = MoleculeData,
            ElementArray = ElementInfo,
            LabelDiffConstraint = LabelDiffConstraint,
            CorrectTracerElementCore = CorrectTracerElementCore,
            CorrectTracerImpurity = CorrectTracerImpurity
          )
          
          CombinedProbList[[MoleculeNo]] <- CombinedProb
        }#MoleculeNo
        
        names(CombinedProbList) <- names(MoleculeInfo3)
      }# (!UltraHighRes)
      
      message(date(),
              " :: preparing <",
              MoleculesTotal,
              " molecules> for correction [OK]")
      
      message("\n", rep("#", options()$width))
      #message(date()," :: all individual molecule stuff completed.")
      message(date(),
              " :: Starting data correction for <",
              SamplesTotal,
              " samples> ...")
      message(rep("#", options()$width), "\n")
      
      # DATA CORRECTION USING THE PROBABILITY MATRIX
      #
      # In the following section, for each sample and molecule from the the input raw data the
      # measured transitions will be corrected using the probability matrix generated before.
      # This is only done if the molecule to be corrected is actually present in the
      # input data and if there are no missing transitions that should be present according to
      # the molecule information file.
      
      if (!CalcOnlyProbabilities) {
        correctionResultList <- list()
        meanEnrichmentResultList <- list()
        
        for (SampleNo in seq_len(SamplesTotal)) {
          message(date(),
                  " :: processing Sample #",
                  SampleNo,
                  " [",
                  colnames(dataRaw)[SampleNo],
                  "] ")
          
          tmpCorrectionResult <- list()
          tmpMeanEnrichmentResult <- list()
          
          for (MoleculeNo in seq_len(MoleculesTotal)) {
            if (UltraHighRes) {
              CombFragmentProb <- CombinedProbList[[MoleculeNo]]
            } else{
              CombFragmentProb <- FragmentsCombinedProbabilityList[[MoleculeNo]]
            }
            
            MoleculeData <- MoleculeInfo3[[MoleculeNo]]
            MoleculeName <- names(MoleculeInfo3[MoleculeNo])
            
            # The function 'RawDataCorrection' performs the actual
            # data correction. It uses the values extracted from the measurement file
            # and corrects molecule by molecule for each sample. This is done by
            # numerically solving a linear equation system. Here,
            # the uncorrected value of a given transition is a linear combination
            # of the corrected transition values with their
            # probability of contributing to the uncorrected value
            # as their coefficients (derived from the probability matrix,
            # CombFragmentProb). The algorithm works with the
            # constraint that a solution of the linear equation
            # system cannot be < 0.
            
            correctionResult <- RawDataCorrection(
              file_num = dataRaw,
              MoleculeData = MoleculeData,
              MoleculeName = MoleculeName,
              CombFragmentProb = CombFragmentProb,
              MoleculeNo = MoleculeNo,
              SampleNo = SampleNo,
              roundDigits = roundDigits,
              logEnvironment = log.env
            )
            
            
            tmpCorrectionResult[[MoleculeNo]] <- correctionResult
            
            #}#length(MoleculeData[["MoleculeLocation"]])!=0
          }#MoleculeNo
          
          names(tmpCorrectionResult) <- names(MoleculeInfo3)
          
          correctionResultList[[SampleNo]] <- tmpCorrectionResult
          
        }#SampleNo
        
        names(correctionResultList) <- colnames(dataRaw)
        
        x <- convert2matrix(input = correctionResultList)
        
        CorrectedDataOutput <- list()
        
        CorrectedDataOutput$RawData <- as.data.frame(dataRaw)
        
        for (name in CorrectedDataNamesBase) {
          CorrectedDataOutput[[name]] <- as.data.frame(x[[name]])
          
        }
        
        if (CalculateMeanEnrichment) {
          meanEnrichmentResultList <-
            CalcMeanEnrichment(
              MoleculeInfo = MoleculeInfo3,
              MoleculesTotal = MoleculesTotal,
              SamplesTotal = SamplesTotal,
              roundDigits = roundDigits,
              correctionResultList =
                correctionResultList,
              UltraHighRes = UltraHighRes
            )
          
          names(meanEnrichmentResultList) <- colnames(dataRaw)
          
          x2 <-
            convertEnrichmentList2matrix(input = meanEnrichmentResultList)
          
          CorrectedDataOutput$MeanEnrichment <- as.data.frame(x2)
          
        }
        
      } else{
        message(
          date(),
          " CalcOnlyProbabilities==TRUE\n=> no correctionResultList generated ...\n"
        )
      }#CalcOnlyProbabilities==1
      
      if (log.env$param$FileOutFormat == "xls") {
        WriteXLS::WriteXLS(
          CorrectedDataOutput,
          file.path(log.env$param$OutputDirectory, log.env$param$FileOut),
          perl = "perl",
          row.names = TRUE
        )
        
      } else if (log.env$param$FileOutFormat == "csv") {
        for (name in CorrectedDataNames) {
          utils::write.csv(
            CorrectedDataOutput[[name]],
            file = file.path(
              log.env$param$OutputDirectory,
              CorrectedDataFileNames[name]
            )
          )
          
        }
        
      }
      
      if (ReturnResultsObject) {
        correctionOutput$results <- CorrectedDataOutput
        
      }
      
      correctionOutput$success = TRUE
      
      message("\n\n", rep("#", options()$width))
      message(date(), " :: IsoCorrection successfully completed!")
      message(rep("#", options()$width), "\n\n")
      
    },
    
    error = function(e) {
      commonErrorString <-
        "\n\nTHE CORRECTION PROCEDURE WAS ABORTED BECAUSE AN ERROR HAS OCCURED. PLEASE CHECK YOUR LOG-FILE.\n\n"
      
      log.env$error <- c(log.env$error, e)
      
      message(paste0(commonErrorString, e))
      
    })
    
    tryCatch({
      writeLog(log.env)
    },
    
    error = function(e) {
      log.env$error <- c(log.env$error, e)
      message(paste0("THE LOG-FILE COULD NOT BE WRITTEN. ERROR: ", e[[1]]))
      
    })
    
    correctionOutput$log <- as.list(log.env)
    
    correctionOutput$error <- log.env$error[1][[1]]
    
    return(correctionOutput)
    
  } #IsoCorrection