# CALCULATION OF THE NON-TRACER ELEMENT COMBINATIONS 

ElementCombinationsNonTracer <- function(ProbElemList, MassShiftElemList, NumberIsoCombEff, 
                                         ElementsNonTracer, NumberElementsNonTracer, CalculationThreshold) {

  # If the number of non-tracer elements is > 1, all isotope combination probabilities of non-tracer
  # elements 1 and 2 of a molecule(-fragment) are multiplied with each other and the corresponding mass shifts are added.  Like in 'IsoCombinations' there
  # is a probability threshold which resets to the outer loop if the probability calculated is below.
  
  CombElemProbList <- list()
  tmpCombElemProb_Array <- vector()
  
  CombElemMassShiftList <- list()
  tmpCombElemMassShift_Array <- vector()
  
  NumberEleComb <- list()
  tmpNumberEleComb <- vector()
  
  if (NumberElementsNonTracer > 1) {
    
    EleComb <- 1
    
    for (IsoCombinationA in seq_len(as.numeric(NumberIsoCombEff[1]))) {
      
      for (IsoCombinationB in seq_len(as.numeric(NumberIsoCombEff[[2]]))) {
        
        tmpCombElemProb <- ProbElemList[[1]][IsoCombinationA] * ProbElemList[[2]][IsoCombinationB]
        
        if (tmpCombElemProb > CalculationThreshold) {
          
          tmpCombElemProb_Array[EleComb] <- tmpCombElemProb
          tmpCombElemMassShift_Array[EleComb] <- MassShiftElemList[[1]][IsoCombinationA] + MassShiftElemList[[2]][IsoCombinationB]
          tmpNumberEleComb <- EleComb
          
          EleComb <- EleComb + 1
        } else {
          IsoCombinationB <- NumberIsoCombEff[[2]] + 1
        }  #tmpCombElemProb > CalculationThreshold
      }  #IsoCombinationB
    }  #IsoCombinationA
    
    CombElemProbList[[2]] <- tmpCombElemProb_Array
    CombElemMassShiftList[[2]] <- tmpCombElemMassShift_Array
    NumberEleComb[[2]] <- tmpNumberEleComb
    
    rm(tmpCombElemProb_Array)
    rm(tmpCombElemMassShift_Array)
    rm(tmpNumberEleComb)
    
  } else if (NumberElementsNonTracer > 0) {
    
    # If the number of non-tracer elements is 1, the probability and mass shift arrays are equal to those generated in 'IsoCombinationsResult'.
    
    tmpCombElemProb_Array[seq_len(NumberIsoCombEff[[1]])] <- ProbElemList[[1]][seq_len(NumberIsoCombEff[[1]])]
    tmpCombElemMassShift_Array[seq_len(NumberIsoCombEff[[1]])] <- MassShiftElemList[[1]][seq_len(NumberIsoCombEff[[1]])]
    tmpNumberEleComb <- NumberIsoCombEff[[1]]
    
    CombElemProbList[[1]] <- tmpCombElemProb_Array
    CombElemMassShiftList[[1]] <- tmpCombElemMassShift_Array
    NumberEleComb[[1]] <- tmpNumberEleComb
    
    rm(tmpCombElemProb_Array)
    rm(tmpCombElemMassShift_Array)
    rm(tmpNumberEleComb)
    
  }  #NumberElementsNonTracer>1
  
  # If there are more than 2 non-tracer elements, the element combinations are calculated in iterations of the following outer loop until all elements have
  # been covered. The algorithm uses the CombElemProbList covering the previously combined element probabilities and multiplies these with the isotope
  # combination probabilities of the next element (ProbElemList).  This is fed back into CombElemProbList at the current element index of the loop). The
  # same is done additively for the mass shift array CombElemMassShiftList.
  
  if (NumberElementsNonTracer > 2) {
    
    tmpCombElemProb_Array <- vector()
    tmpCombElemMassShift_Array <- vector()
    tmpNumberEleComb <- vector()
    
    for (NonTracer in 3:NumberElementsNonTracer) {
      EleComb <- 1
      for (ElementCombination in seq_len(NumberEleComb[[NonTracer - 1]])) {
        for (IsoCombination in seq_len(NumberIsoCombEff[[NonTracer]])) {
          
          tmpCombElemProb <- CombElemProbList[[NonTracer - 1]][[ElementCombination]] * ProbElemList[[NonTracer]][[IsoCombination]]
          if (tmpCombElemProb > CalculationThreshold) {
            
            tmpCombElemProb_Array[EleComb] <- tmpCombElemProb
            tmpCombElemMassShift_Array[EleComb] <- CombElemMassShiftList[[NonTracer - 1]][[ElementCombination]] + MassShiftElemList[[NonTracer]][[IsoCombination]]
            tmpNumberEleComb <- EleComb
            EleComb <- EleComb + 1
            
          } else {
            
            IsoCombination <- NumberIsoCombEff[[NonTracer]] + 1
            
          }  #tmpCombElemProb > CalculationThreshold
        }  #IsoCombination
        
      }  #ElementCombination
      CombElemProbList[[NonTracer]] <- tmpCombElemProb_Array
      CombElemMassShiftList[[NonTracer]] <- tmpCombElemMassShift_Array
      NumberEleComb[[NonTracer]] <- tmpNumberEleComb
      
    }  #NonTracer
  }  #NumberElementsNonTracer>2
  
  return(list("CombElemProbList"=CombElemProbList,
              "CombElemMassShiftList"=CombElemMassShiftList,
              "NumberEleComb"=NumberEleComb))
  
}