# Calculates the sum of probabilities with a common mass shift

#' @importFrom stringr str_c
#' 
MassShiftProbabilities <- function(MoleculeArray, MoleculeNo, Fragment, MaxMassShift, ElementCombinations, CorrectTracerImpurity) {
  
  MoleculeData <- MoleculeArray[[MoleculeNo]][[stringr::str_c("Fragment_", Fragment)]]
  
  Transitions <- MoleculeArray[[MoleculeNo]][["Transitions"]]
  NumberTransitions <- nrow(Transitions)
  NumberTracers <- length(MoleculeData[["Tracer"]])
  
  if (CorrectTracerImpurity) {
    if (NumberTracers > 0) {
      
      Probability_Array <- ElementCombinations[["TracerImpurityCombProb_Array"]]
      MassShift_Array <- ElementCombinations[["TracerImpurityCombMassShift_Array"]]
      
    } else {
      Probability_Array <- ElementCombinations[["TracerElemProb_Array"]]
      MassShift_Array <- ElementCombinations[["TracerElemMassShift_Array"]]
      
    }
    
  } else {
    
    Probability_Array <- ElementCombinations[["TracerElemProb_Array"]]
    MassShift_Array <- ElementCombinations[["TracerElemMassShift_Array"]]
    
  }
  
  # This loop goes through all mass shifts until it reaches the maximum possible
  # mass shift of the current molecule(-fragment).  For each labelling state,
  # it finds indexes in MassShift_Array associated with the given mass shift 
  # and uses them to find the corresponding probabilites in Probability_Array.
  # Those values are eventually summed up to yield the probability that a 
  # certain labelling state produces a certain mass shift in relation to the
  # unlabelled species due to natural abundance and possibly tracer impurity.
  
  CumProbList <- list()
  
  for (MassShift in 0:MaxMassShift) {
    
    tmpCumProb.vec <- vector()
    for (TransitionNo in seq_len(NumberTransitions)) {
      
      tmpIdx <- which(MassShift_Array[[TransitionNo]] == MassShift)
      
      tmpCumProb.vec[TransitionNo] <- sum(Probability_Array[[TransitionNo]][tmpIdx])
      
    }
    CumProbList[[MassShift + 1]] <- tmpCumProb.vec
  }  #MassShift
  
  names(CumProbList) <- stringr::str_c("MassShift_", 0:MaxMassShift)
  
  return(CumProbList)
}
