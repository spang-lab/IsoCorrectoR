context("End to end testing of IsoCorrectoR")

fullTest <- TRUE

TestDataPath <- system.file('testdata', package='IsoCorrectoR')
ElementFile = file.path(TestDataPath, 'ElementFile.csv')

endToEndTest <- function(PathToProjects, Project, ElementFilePath,
                         ResolutionOption=FALSE, ExpectedOutcome='TRUE') {
  
  ProjectPath <- file.path(PathToProjects, Project)
  
  CalculateMeanEnrichment = TRUE
  CorrectAlsoMonoisotopic = TRUE
  CalculationThreshold = 10^(-8)
  CalculationThreshold_UHR = 8
  FileOutFormat='csv'
  
  ImpurityOptions <- c(TRUE, FALSE)
  CoreOptions <- c(TRUE, FALSE)
  names(ImpurityOptions) <- c('p', 'np')
  names(CoreOptions) <- c('co', 'nco')
  
  resultTypes <- c('c', 'cf', 'me', 're', 'cm', 'cmf', 'rd')
  names(resultTypes) <- c('Corrected', 'CorrectedFractions', 'MeanEnrichment', 'Residuals',
                          'CorrectedMonoisotopic', 'CorrectedMonoisotopicFractions', 'RawData')
  
  DirOut = file.path(ProjectPath, 'test_output')
  
  if(!dir.exists(DirOut)) {
    
    dir.create(DirOut)
    
  }
  
  for (ImpurityOption in names(ImpurityOptions)) {
    
    for (CoreOption in names(CoreOptions)) {
      
      MeasurementFile = file.path(ProjectPath, 'MeasurementFile.csv')
      MoleculeFile = file.path(ProjectPath, 'MoleculeFile.csv')
      
      FileOut = paste0(CoreOption, '_', ImpurityOption)
      
      #run IsoCorrection with given parameters
      
      correctionResults <- IsoCorrection(MeasurementFile=MeasurementFile, ElementFile=ElementFilePath, MoleculeFile=MoleculeFile,
                                         CorrectTracerImpurity=ImpurityOptions[ImpurityOption], CorrectTracerElementCore=CoreOptions[CoreOption], 
                                         CalculateMeanEnrichment=CalculateMeanEnrichment,
                                         UltraHighRes=ResolutionOption,
                                         DirOut=DirOut, 
                                         FileOut=FileOut,
                                         FileOutFormat=FileOutFormat,
                                         CorrectAlsoMonoisotopic=CorrectAlsoMonoisotopic,
                                         CalculationThreshold=CalculationThreshold,
                                         CalculationThreshold_UHR=CalculationThreshold_UHR,
                                         ReturnResultsObject=TRUE,
                                         Testmode=TRUE)
      
      #Load reference results
      
      success <- correctionResults$success
      results <- correctionResults$results
      
      testinfo <- paste0('Dataset: ', Project, ', purity: ', ImpurityOptions[ImpurityOption], ', core: ', CoreOptions[CoreOption]) 
      
      expect_equal(success, ExpectedOutcome, info = testinfo)
      
      for(resultType in names(resultTypes)) {
        
        reference <- read.csv(file.path(ProjectPath, 'reference', paste0(FileOut, '_', resultTypes[resultType], '.csv')), row.names = 1)
        testresults <- read.csv(file.path(ProjectPath, 'test_output', paste0('IsoCorrectoR_', FileOut, '_', resultType, '.csv')), row.names = 1)
        
        testinfoSpecific <- paste0(testinfo, ', result type: ', resultType)
        
        expect_equal(testresults, reference, info = testinfoSpecific)
        
      }
    }
  }
  
  #Cleanup
  
  unlink(DirOut, recursive = TRUE)
  
}

test_that("End to end tests for normal resolution", {
  
  if(!fullTest) {
  
    skip('Skipped because no full test')
    
  }
  
  generalPath <- file.path(TestDataPath, 'end_to_end_nres')

  ProjectFolders <- list.files(generalPath)
  
  for (Project in ProjectFolders) {
    
    endToEndTest(PathToProjects=generalPath, Project=Project, ElementFilePath=ElementFile,
                 ResolutionOption=FALSE)
    
  }
     
})

test_that("End to end tests for high resolution", {
  
  if(!fullTest) {
    
    skip('Skipped because no full test')
    
  }
  
  generalPath <- file.path(TestDataPath, 'end_to_end_hres')
  
  ProjectFolders <- list.files(generalPath)
  
  for (Project in ProjectFolders) {
    
    endToEndTest(PathToProjects=generalPath, Project=Project, ElementFilePath=ElementFile,
                             ResolutionOption=TRUE)
      
  }
  
})

test_that("End to end tests for normal resolution with missing values", {
  
  generalPath <- file.path(TestDataPath, 'missing_values_nres')
  
  ProjectFolders <- list.files(generalPath)
  
  for (Project in ProjectFolders) {
    
    suppressWarnings(endToEndTest(PathToProjects=generalPath, Project=Project, ElementFilePath=ElementFile,
                 ResolutionOption=FALSE, ExpectedOutcome = "WARNINGS"))
    
  }
  
})
