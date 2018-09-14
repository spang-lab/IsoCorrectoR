# Calculates the sum of probabilities with a common mass shift

# The function 'MassShiftProbabilities' yields the
# probability that a given transition of a
# molecule(-fragment) (through naturally occuring isotopes)
# produces a certain mass shift in relation to the same molecule(-fragment)
# containing no isotopes of higher mass. This probability is gained
# by summing all entries of a molecule(-fragment) of a given
# transition in TracerElemProbList that have the same mass shift in
# TracerElemMassShiftList. If correction for tracer impurity is turned
# on and if the molecule(-fragment) in question can contain tracer,
# 'MassShiftProbabilities' uses TracerImpurityCombProbList
# and TracerImpurityCombMassShiftList instead.

#' @importFrom stringr str_c
#' 
MassShiftProbabilities <- function(MoleculeInfo, MoleculeNo, Fragment, MaxMassShift, ElementCombinations, CorrectTracerImpurity) {
  
  MoleculeData <- MoleculeInfo[[MoleculeNo]][[stringr::str_c("Fragment_", Fragment)]]
  
  Transitions <- MoleculeInfo[[MoleculeNo]][["Transitions"]]
  NumberTransitions <- nrow(Transitions)
  NumberTracers <- length(MoleculeData[["Tracer"]])
  
  if (CorrectTracerImpurity) {
    if (NumberTracers > 0) {
      
      ProbabilityList <- ElementCombinations[["TracerImpurityCombProbList"]]
      MassShiftList <- ElementCombinations[["TracerImpurityCombMassShiftList"]]
      
    } else {
      ProbabilityList <- ElementCombinations[["TracerElemProbList"]]
      MassShiftList <- ElementCombinations[["TracerElemMassShiftList"]]
      
    }
    
  } else {
    
    ProbabilityList <- ElementCombinations[["TracerElemProbList"]]
    MassShiftList <- ElementCombinations[["TracerElemMassShiftList"]]
    
  }
  
  # This loop goes through all mass shifts until it reaches the maximum possible
  # mass shift of the current molecule(-fragment).  For each labelling state,
  # it finds indexes in MassShiftList associated with the given mass shift 
  # and uses them to find the corresponding probabilites in ProbabilityList.
  # Those values are eventually summed up to yield the probability that a 
  # certain labelling state produces a certain mass shift in relation to the
  # unlabelled species due to natural abundance and possibly tracer impurity.
  
  CumProbList <- list()
  
  for (MassShift in 0:MaxMassShift) {
    
    tmpCumProb.vec <- vector()
    for (TransitionNo in seq_len(NumberTransitions)) {
      
      tmpIdx <- which(MassShiftList[[TransitionNo]] == MassShift)
      
      tmpCumProb.vec[TransitionNo] <- sum(ProbabilityList[[TransitionNo]][tmpIdx])
      
    }
    CumProbList[[MassShift + 1]] <- tmpCumProb.vec
  }  #MassShift
  
  names(CumProbList) <- stringr::str_c("MassShift_", 0:MaxMassShift)
  
  return(CumProbList)
}
