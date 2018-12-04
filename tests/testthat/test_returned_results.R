context("Test correctness of results returned by the IsoCorrection() function")

test_that("Test results returned by IsoCorrection()", {
  
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
                                     CorrectAlsoMonoisotopic=TRUE)
  
  success <- correctionResults$success
  results <- correctionResults$results
  
  testinfo <- paste0('Dataset: ', Project) 
  
  expect_equal(success, ExpectedOutcome, info = testinfo)
  
  #Load reference results and test results from file and compare
  
  for(resultType in names(resultTypes)) {
    
    reference <- loadResults(filepath=file.path(ProjectPath, 'reference', paste0(FileOut, '_', resultTypes[resultType], '.csv')),
                             filetype='csv', dataDescription='reference')
    
    testresults <- results[[resultType]]
    
    testinfoSpecific <- paste0(testinfo, ', result type: ', resultType)
    
    expect_equal(testresults, reference, info = testinfoSpecific)
    
  }
  
  #Cleanup
  
  unlink(DirOut, recursive = TRUE)
     
})