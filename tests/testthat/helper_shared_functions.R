library(readxl)

#Variables used by all/many tests

TestDataPath <- system.file('testdata', package='IsoCorrectoR')

resultTypes <- c('c', 'cf', 'me', 'cm', 'rd')
names(resultTypes) <- c('Corrected', 'CorrectedFractions', 'MeanEnrichment',
                        'CorrectedMonoisotopic', 'RawData')

#Safe loading of result files with tryCatch()

loadResults <- function(filepath, filetype='csv', sheet=NA, dataDescription=NA) {
  
  result <- tryCatch(
    {
      
      if(filetype=='csv') {
      
        result <- read.csv(filepath, row.names = 1)
        
      } else if (filetype=='xls') {
        
        result <- as.data.frame(read_excel(filepath, col_names = TRUE, sheet = sheet))
        
        rownames(result) <- result[[1]]
        result <- result[-1]
        
      } else {
        
        stop("Invalid file extension.")
        
      }
      
    },
    error = function(e) {
      
      result <- paste0("No ", dataDescription, " file")
      return(result)
      
    }
  )
  
  return(result)
  
}

#Safe call to IsoCorrection() with tryCatch

IsoCorrectionSafeCall <- function(MeasurementFile=NA,
                                  ElementFile=NA,
                                  MoleculeFile=NA,
                                  CorrectTracerImpurity=FALSE,
                                  CorrectTracerElementCore=TRUE, 
                                  CalculateMeanEnrichment=TRUE,
                                  UltraHighRes=FALSE,
                                  DirOut='.', 
                                  FileOut='result',
                                  FileOutFormat='csv',
                                  CorrectAlsoMonoisotopic=FALSE,
                                  CalculationThreshold=10^(-8),
                                  CalculationThreshold_UHR=8,
                                  ReturnResultsObject=TRUE,
                                  Testmode=TRUE) {

  correctionResults <- tryCatch(
    {
      
      correctionResults <- IsoCorrection(MeasurementFile=MeasurementFile,
                                         ElementFile=ElementFile,
                                         MoleculeFile=MoleculeFile,
                                         CorrectTracerImpurity=CorrectTracerImpurity,
                                         CorrectTracerElementCore=CorrectTracerElementCore, 
                                         CalculateMeanEnrichment=CalculateMeanEnrichment,
                                         UltraHighRes=UltraHighRes,
                                         DirOut=DirOut, 
                                         FileOut=FileOut,
                                         FileOutFormat=FileOutFormat,
                                         CorrectAlsoMonoisotopic=CorrectAlsoMonoisotopic,
                                         CalculationThreshold=CalculationThreshold,
                                         CalculationThreshold_UHR=CalculationThreshold_UHR,
                                         ReturnResultsObject=ReturnResultsObject,
                                         Testmode=Testmode)
      
    }, 
    error = function(e) {
      
      correctionResults <- list()
      
      correctionResults$success <- "UNEXPECTED_ERROR"
      correctionResults$results <- NA
      
      return(correctionResults)
      
    }
  )
}