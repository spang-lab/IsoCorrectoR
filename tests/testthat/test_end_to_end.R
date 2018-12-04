context("End to end testing of IsoCorrectoR with csv files")

ElementFile = file.path(TestDataPath, 'ElementFile.csv')

endToEndTest <- function(PathToProjects, Project, ElementFilePath,
                         ResolutionOption=FALSE, ExpectedOutcome='TRUE', resultTypes=NA) {
  
  ProjectPath <- file.path(PathToProjects, Project)
  
  MeasurementFile = file.path(ProjectPath, 'MeasurementFile.csv')
  MoleculeFile = file.path(ProjectPath, 'MoleculeFile.csv')
  
  ImpurityOptions <- c(TRUE, FALSE)
  CoreOptions <- c(TRUE, FALSE)
  names(ImpurityOptions) <- c('p', 'np')
  names(CoreOptions) <- c('co', 'nco')
  
  DirOut = file.path(ProjectPath, 'test_output')
  
  if(!dir.exists(DirOut)) {
    
    dir.create(DirOut)
    
  }
  
  for (ImpurityOption in names(ImpurityOptions)) {
    
    for (CoreOption in names(CoreOptions)) {
      
      FileOut = paste0(CoreOption, '_', ImpurityOption)
      
      #run IsoCorrection with given parameters
      
      correctionResults <- IsoCorrectionSafeCall(MeasurementFile=MeasurementFile,
                                         ElementFile=ElementFilePath,
                                         MoleculeFile=MoleculeFile,
                                         CorrectTracerImpurity=ImpurityOptions[ImpurityOption],
                                         CorrectTracerElementCore=CoreOptions[CoreOption], 
                                         UltraHighRes=ResolutionOption,
                                         DirOut=DirOut, 
                                         FileOut=FileOut,
                                         CorrectAlsoMonoisotopic=TRUE,
                                         Testmode=TRUE)
      
      success <- correctionResults$success
      results <- correctionResults$results
      
      testinfo <- paste0('Dataset: ', Project, ', purity: ', ImpurityOptions[ImpurityOption], ', core: ', CoreOptions[CoreOption]) 
      
      expect_equal(success, ExpectedOutcome, info = testinfo)
      
      #Load reference results and test results from file and compare
      
      for(resultType in names(resultTypes)) {
        
        reference <- loadResults(filepath=file.path(ProjectPath, 'reference', paste0(FileOut, '_', resultTypes[resultType], '.csv')),
                                 filetype='csv', dataDescription='reference')
        
        testresults <- loadResults(filepath=file.path(ProjectPath, 'test_output', paste0('IsoCorrectoR_', FileOut, '_', resultType, '.csv')),
                                   filetype='csv', dataDescription='testresults')
        
        testinfoSpecific <- paste0(testinfo, ', result type: ', resultType)
        
        expect_equal(testresults, reference, info = testinfoSpecific)
        
      }
    }
  }
  
  #Cleanup
  
  unlink(DirOut, recursive = TRUE)
  
}

test_that("End to end tests for normal resolution", {
  
  generalPath <- file.path(TestDataPath, 'end_to_end_nres')

  ProjectFolders <- list.files(generalPath)
  
  for (Project in ProjectFolders) {
    
    endToEndTest(PathToProjects=generalPath, Project=Project, ElementFilePath=ElementFile,
                 ResolutionOption=FALSE, resultTypes = resultTypes)
    
  }
     
})

test_that("End to end tests for high resolution", {
  
  generalPath <- file.path(TestDataPath, 'end_to_end_hres')
  
  ProjectFolders <- list.files(generalPath)
  
  for (Project in ProjectFolders) {
    
    endToEndTest(PathToProjects=generalPath, Project=Project, ElementFilePath=ElementFile,
                             ResolutionOption=TRUE, resultTypes = resultTypes)
      
  }
  
})

test_that("End to end tests for normal resolution with missing values", {
  
  generalPath <- file.path(TestDataPath, 'missing_values_nres')
  
  ProjectFolders <- list.files(generalPath)
  
  for (Project in ProjectFolders) {
    
    suppressWarnings(endToEndTest(PathToProjects=generalPath, Project=Project, ElementFilePath=ElementFile,
                 ResolutionOption=FALSE, ExpectedOutcome = "WARNINGS", resultTypes = resultTypes))
    
  }
  
})

test_that("End to end tests for normal resolution, additional testsets", {
    
  skip('Skipped')
  
  generalPath <- file.path(TestDataPath, 'additional_testsets_nres')
  
  ProjectFolders <- list.files(generalPath)
  
  for (Project in ProjectFolders) {
    
    endToEndTest(PathToProjects=generalPath, Project=Project, ElementFilePath=ElementFile,
                 ResolutionOption=FALSE, resultTypes = resultTypes)
    
  }
  
})
