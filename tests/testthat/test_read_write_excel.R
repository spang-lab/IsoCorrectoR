context("Test reading and writing of excel files by IsoCorrectoR")

test_that("Test reading of xls files", {
  
  ExpectedOutcome = 'TRUE'
  
  Project <- 'AA1'
  
  ProjectPath <- file.path(TestDataPath, 'end_to_end_nres', Project)
  
  DirOut = file.path(ProjectPath, 'test_output')
  
  if(!dir.exists(DirOut)) {
    
    dir.create(DirOut)
    
  }
  
  MeasurementFile = file.path(ProjectPath, 'MeasurementFile.xlsx')
  MoleculeFile = file.path(ProjectPath, 'MoleculeFile.xlsx')
  ElementFile = file.path(TestDataPath, 'ElementFile.xlsx')
  
  FileOut = "co_np"
  
  #run IsoCorrection with given parameters
      
  correctionResults <- IsoCorrectionSafeCall(MeasurementFile=MeasurementFile,
                                     ElementFile=ElementFile,
                                     MoleculeFile=MoleculeFile,
                                     DirOut=DirOut, 
                                     FileOut=FileOut,
                                     FileOutFormat="csv",
                                     CorrectAlsoMonoisotopic=TRUE)
  
  success <- correctionResults$success
  results <- correctionResults$results
  
  testinfo <- paste0('Dataset: ', Project) 
  
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
  
  #Cleanup
  
  unlink(DirOut, recursive = TRUE)
     
})

test_that("Test writing of xls files", {
  
  skip('Skipped')
  
  ExpectedOutcome = 'TRUE'
  
  Project <- 'AA1'
  
  ProjectPath <- file.path(TestDataPath, 'end_to_end_nres', Project)
  
  DirOut = file.path(ProjectPath, 'test_output')
  
  if(!dir.exists(DirOut)) {
    
    dir.create(DirOut)
    
  }
  
  MeasurementFile = file.path(ProjectPath, 'MeasurementFile.csv')
  MoleculeFile = file.path(ProjectPath, 'MoleculeFile.csv')
  ElementFile = file.path(TestDataPath, 'ElementFile.csv')
  
  FileOut = "co_np"
  
  #run IsoCorrection with given parameters
  
  correctionResults <- IsoCorrectionSafeCall(MeasurementFile=MeasurementFile,
                                             ElementFile=ElementFile,
                                             MoleculeFile=MoleculeFile,
                                             DirOut=DirOut, 
                                             FileOut=FileOut,
                                             FileOutFormat="xls",
                                             CorrectAlsoMonoisotopic=TRUE)
  
  success <- correctionResults$success
  results <- correctionResults$results
  
  testinfo <- paste0('Dataset: ', Project) 
  
  expect_equal(success, ExpectedOutcome, info = testinfo)
  
  #Load reference results and test results from file and compare
  
  for(resultType in names(resultTypes)) {
    
    reference <- loadResults(filepath=file.path(ProjectPath, 'reference', paste0(FileOut, '_', resultTypes[resultType], '.csv')),
                             filetype='csv', dataDescription='reference')
    
    testresults <- loadResults(filepath=file.path(ProjectPath, 'test_output', paste0('IsoCorrectoR_', FileOut, '.xls')),
                               filetype='xls', sheet = resultType, dataDescription='testresults')
    
    testinfoSpecific <- paste0(testinfo, ', result type: ', resultType)
    
    expect_equal(testresults, reference, info = testinfoSpecific)
    
  }
  
  #Cleanup
  
  unlink(DirOut, recursive = TRUE)
  
})